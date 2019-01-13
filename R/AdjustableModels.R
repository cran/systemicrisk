#' Calibrate empirical fitness model to a given density
#'
#' The model is an empirical fitness based model for the existence of
#' links (more details below) which contains one fixed parameter and
#' the weight of each existing link follows an exponential
#' distribution with a fixed rate parameter. This function chooses the
#' two parameters such that the density of the network (the average
#' proportion of existing links) with these given row and column sums
#' is a certain desired value.
#'
#' The empirical fitness model assumes that every node
#' \eqn{f[i]=log(l[i]+a[i])} has a fitness given by the observered row
#' and column sum and that the existence probability of a link between
#' node i and j is then given by
#' \eqn{p[i,j]=1/(1+exp(-(alpha+f[i]+f[j])))}, where alpha is an
#' additional parameter.  The resulting model uses observed quantities
#' (the row and column sums of the matrix) as input to the model and
#' is thus an empirical Bayes approach.
#'  
#' @param l row sums of matrix to be reconstructed
#' @param a column sum of matrix to be reconstructed 
#' @param targetdensity desired proportion of reconstructed entries to be positive  
#' @param L_fixed Matrix containing known values of L, where NA signifies that
#'                an element is not known. If \code{L_fixed} equates to \code{NA} (the
#'                default) then no values are assumed to be known.
#' @param nsamples_calib number of matrices to generate during calibration.
#' @param thin_calib amount of thinning to use during calibration 
#' @return Model that can be used to generate the desired samples using \code{\link{sample_HierarchicalModel}}.
#' @examples
#' ## first generate a true network
#' n <- 10 # size of network
#' ftrue <- rnorm(n) # vector of underlying fitnesses
#' p <- outer(ftrue,ftrue,FUN=function(x,y) 1/(1+exp(-(x+y))))
#' lambda <- 0.1
#' L <- matrix(nrow=n,rbinom(n*n,prob=p,size=1)*rexp(n*n,rate=lambda))
#' 
#' # then reconstruct with a target density of 0.7
#' model <- calibrate_FitnessEmp(l=rowSums(L),a=colSums(L),
#'                               targetdensity=0.7,nsamples_calib=10,thin_calib=50)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),model=model,
#'                                     nsamples=10,thin=1e2)
#' # check row sums
#' rowSums(L)
#' rowSums(Lsamp$L[[10]])
#' # check calibration
#' mean(Lsamp$L[[10]]>0)
#'
#' # now an example with some fixed entries
#' L_fixed <- L
#' L_fixed[1:(n/2),] <- NA
#' # then reconstruct with a target density of 0.9
#' model <- calibrate_FitnessEmp(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                               targetdensity=0.9,nsamples_calib=10,thin_calib=50)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                                   model=model,nsamples=10,thin=1e2)
#' mean(Lsamp$L[[10]][-(1:(n/2)),]>0) # known entries
#' mean(Lsamp$L[[10]][(1:(n/2)),]>0)  #reconstructed entries
#' @export
calibrate_FitnessEmp <- function(l,a,targetdensity,L_fixed=NA,
                                 nsamples_calib=100,thin_calib=1e2){
    n <- length(l)
    fitness <- log((l+a))
    if (is.matrix(L_fixed)) 
        n_free <- sum(is.na(L_fixed))
    else
        n_free <- n*n
    lambda <- 1/(sum(l)/(targetdensity*n_free))
    additive <- max(fitness)
    fitness <- fitness-mean(fitness)

    f <- function(alpha){
        cat("alpha=",alpha,"\n")
        pmatr <- outer(fitness,fitness,FUN=function(x,y) 1/(1+exp(-(alpha+x+y))))
        print(range(pmatr))
        model<-Model.Indep.p.lambda(model.p=Model.p.constant(n=n, p=pmatr), model.lambda=Model.lambda.constant(lambda=lambda, n=n))
        Lsamp <- sample_HierarchicalModel(l,a,L_fixed=L_fixed,model=model,
                                          nsamples=nsamples_calib,thin=thin_calib)
        mean(sapply(Lsamp$L,function(x) {
            if (is.matrix(L_fixed))
                mean(ifelse(is.na(L_fixed),x>0,NA),na.rm=TRUE)
            else sum(x>0)/n_free
        }))-targetdensity
    }
    res <- uniroot(f,interval=c(-10,10),extendInt="yes",tol=0.01,maxiter=20)

    alpha <- res$root
    cat("alpha=",alpha,"\n")
    pmatr <- outer(fitness,fitness,FUN=function(x,y) 1/(1+exp(-(alpha+x+y))))
    Model.Indep.p.lambda(model.p=Model.p.constant(n=n, p=pmatr), model.lambda=Model.lambda.constant(lambda=lambda, n=n))
}

#' Calibrate ER model to a given density
#'
#' The model is an Erdos-Renyi model for the existence of links (a
#' link exists independently of other links with a fixed probability)
#' and the weight of each existing link follows an exponential
#' distribution with a fixed rate parameter. This function chooses the
#' two parameters such that the density of the network (the average
#' proportion of existing links) is a certain desired value. Diagonal
#' elements are being set to 0.
#'  
#' @inheritParams calibrate_FitnessEmp
#' @return Model that can be used to generate the desired samples using \code{\link{sample_HierarchicalModel}}.
#' @examples
#' ## first generate a true network
#' n <- 10 # size of network
#' p <- 0.45
#' lambda <- 0.1
#' L <- matrix(nrow=n,rbinom(n*n,prob=p,size=1)*rexp(n*n,rate=lambda))
#' 
#' # then reconstruct with a target density of 0.55
#' model <- calibrate_ER(l=rowSums(L),a=colSums(L),
#'                       targetdensity=0.55,nsamples_calib=10)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),model=model,
#'                                     nsamples=10,thin=1e2)
#' # check row sums
#' rowSums(L)
#' rowSums(Lsamp$L[[10]])
#' # check calibration
#' mean(Lsamp$L[[10]]>0)
#'
#' # now an example with some fixed entries
#' L_fixed <- L
#' L_fixed[1:(n/2),] <- NA
#' # then reconstruct with a target density of 0.9
#' model <- calibrate_ER(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                               targetdensity=0.9,nsamples_calib=10)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                                   model=model,nsamples=10,thin=1e2)
#' mean(Lsamp$L[[10]][-(1:(n/2)),]>0) # known entries
#' mean(Lsamp$L[[10]][(1:(n/2)),]>0)  #reconstructed entries
#' @export
calibrate_ER <- function(l,a,targetdensity,L_fixed=NA,
                         nsamples_calib=100,thin_calib=1e2){
    n <- length(l)
    if (is.matrix(L_fixed)) 
        n_free <- sum(is.na(L_fixed))
    else
        n_free <- n*n
    lambda <- 1/(sum(l)/(targetdensity*n_free))
    f <- function(p){
        cat("p=",p,"\n")
        model<-Model.Indep.p.lambda(model.p=Model.p.constant(n=n, p=p), model.lambda=Model.lambda.constant(lambda=lambda, n=n))
        Lsamp <- sample_HierarchicalModel(l,a,L_fixed=L_fixed,model=model,
                                          nsamples=nsamples_calib,thin=1e2)
        mean(sapply(Lsamp$L,function(x) {
            if (is.matrix(L_fixed))
                mean(ifelse(is.na(L_fixed),x>0,NA),na.rm=TRUE)
            else sum(x>0)/n_free
        }))-targetdensity
    }
    res <- uniroot(f,interval=c(0,1),f.lower=-targetdensity,f.upper=1-targetdensity,tol=0.01,maxiter=20)
    
    p <- res$root
    cat("p=",p,"\n")
    Model.Indep.p.lambda(model.p=Model.p.constant(n=n, p=p), model.lambda=Model.lambda.constant(lambda=lambda, n=n))
}


#' Calibrate ER model to a given density with a nonsquare matrix
#'
#' The model is an Erdos-Renyi model for the existence of links (a
#' link exists independently of other links with a fixed probability)
#' and the weight of each existing link follows an exponential
#' distribution with a fixed rate parameter. This function chooses the
#' two parameters such that the density of the network (the average
#' proportion of existing links) is a certain desired value. This
#' function does not set diagonal values to 0.
#'
#'  
#' @inheritParams calibrate_FitnessEmp
#' @return Model that can be used to generate the desired samples using \code{\link{sample_HierarchicalModel}}.
#' @examples
#' ## first generate a true network
#' nrow <- 10 # size of network
#' ncol <- 8 # size of network
#' p <- 0.45
#' lambda <- 0.1
#' L <- matrix(nrow=nrow,rbinom(nrow*ncol,prob=p,size=1)*rexp(nrow*ncol,rate=lambda))
#' 
#' # then reconstruct with a target density of 0.55
#' model <- calibrate_ER.nonsquare(l=rowSums(L),a=colSums(L),
#'                       targetdensity=0.55,nsamples_calib=10)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),model=model,
#'                                     nsamples=10,thin=1e2)
#' # check row sums
#' rowSums(L)
#' rowSums(Lsamp$L[[10]])
#' # check calibration
#' mean(Lsamp$L[[10]]>0)
#'
#' # now an example with some fixed entries
#' L_fixed <- L
#' L_fixed[1:(nrow/2),] <- NA
#' # then reconstruct with a target density of 0.9
#' model <- calibrate_ER.nonsquare(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                               targetdensity=0.9,nsamples_calib=10)
#' Lsamp <- sample_HierarchicalModel(l=rowSums(L),a=colSums(L),L_fixed=L_fixed,
#'                                   model=model,nsamples=10,thin=1e2)
#' mean(Lsamp$L[[10]][-(1:(nrow/2)),]>0) # known entries
#' mean(Lsamp$L[[10]][(1:(nrow/2)),]>0)  #reconstructed entries
#' @export
calibrate_ER.nonsquare <- function(l,a,targetdensity,L_fixed=NA,
                         nsamples_calib=100,thin_calib=1e2){
    nrow <- length(l)
    ncol <- length(a)
    if (is.matrix(L_fixed)) 
        n_free <- sum(is.na(L_fixed))
    else
        n_free <- nrow*ncol
    lambda <- 1/(sum(l)/(targetdensity*n_free))
    f <- function(p){
        cat("p=",p,"\n")
        model<-Model.Indep.p.lambda(model.p=Model.p.constant.nonsquare(nrow=nrow,ncol=ncol, p=p),
                                    model.lambda=Model.lambda.constant.nonsquare(lambda=lambda, nrow=nrow, ncol=ncol))
        Lsamp <- sample_HierarchicalModel(l,a,L_fixed=L_fixed,model=model,
                                          nsamples=nsamples_calib,thin=1e2)
        mean(sapply(Lsamp$L,function(x) {
            if (is.matrix(L_fixed))
                mean(ifelse(is.na(L_fixed),x>0,NA),na.rm=TRUE)
            else sum(x>0)/n_free
        }))-targetdensity
    }
    res <- uniroot(f,interval=c(0,1),f.lower=-targetdensity,f.upper=1-targetdensity,tol=0.01,maxiter=20)
    
    p <- res$root
    cat("p=",p,"\n")
    Model.Indep.p.lambda(model.p=Model.p.constant.nonsquare(nrow=nrow,ncol=ncol, p=p),
                         model.lambda=Model.lambda.constant.nonsquare(lambda=lambda, nrow=nrow, ncol=ncol))
}
