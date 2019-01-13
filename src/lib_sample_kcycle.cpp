#include<Rcpp.h>
using namespace Rcpp;
using std::vector;
RNGScope scope;

//// helper functions
double min(double a,double b){
  if (a<b) return a;
  return b;
}
double max(double a,double b){
  if (a>b) return a;
  return b;
}
inline int randWrapper(const int n) { return floor(unif_rand()*n); }//for resampling with STL

//' Creates a deep copy of a matrix
//'
//' Useful  when calling \code{\link{ERE_step_cycle}} or
//' \code{\link{GibbsSteps_kcycle}} to ensure that
//' there are no side effects for the return values.
//' @param M A matrix
//' @return A deep copy of the matrix.
//' @examples
//' lambda <- matrix(0.5,nrow=2,ncol=2)
//' p <- matrix(0.7, nrow=2,ncol=2)
//' L <- matrix(rexp(4),nrow=2);
//' L
//' Lold <- L
//' Lcopy <- cloneMatrix(L)
//' ERE_step_cycle(r=c(0,1),c=c(0,1),L=L,lambda=lambda,p=p)
//'
//' L     ## new value
//' Lold  ## equal to L !!!
//' Lcopy ## still has the original value
//' @export
// [[Rcpp::export]]
NumericMatrix cloneMatrix(NumericMatrix M){
  NumericMatrix res=Rcpp::clone(M);
  return res;
}


inline double roundnum(double x, double eps){
  if (fabs(x)<eps)
    return(0.);
  else
    return(x);
}


//////////////////General kcycle
void setresDelta(NumericMatrix& L, vector<int> r, vector<int> c,  double Delta, double eps){
  int k=c.size();
  for (int i=0; i<k; i++){
    L(r[i],c[i])=roundnum(L(r[i],c[i])+Delta,eps);
    L(r[i],c[(i+1)%k])=roundnum(L(r[i],c[(i+1)%k])-Delta,eps);
  }
}

double loglDelta(NumericMatrix L, NumericMatrix lambda, NumericMatrix p, std::vector<int> r, std::vector<int> c, double Delta, double eps){
  double res=0.;
  int k=c.size();
  for (int i=0; i<k; i++){
    for (int j=0; j<2; j++){
      int ir=r[i];
      int ic=c[(i+j)%k];
      double val=L(ir,ic);
      if (j==1) val-=Delta;
      else val+=Delta;
      if (val<eps){
	res+=log(1-p(ir,ic));
      }else{
	res+=log(p(ir,ic))+log(lambda(ir,ic))-lambda(ir,ic)*val;
      }
    }
  }
  return(res);
}

//' @title Does one Gibbs Step on a cycle
//'
//' @description Execute one Gibbs step on a cycle keeping
//' row and column sums fixed
//'
//' @param r Row indies of cycle, starting at 0 (vector of length k)
//' @param c Column indices of cycle, starting at 0 (vector of length k)
//' @param L nxn matrix with nonnegative values (will be modified)
//' @param lambda nxn matrix of intensities
//' @param p nxn matrix of probabilities (must be in [0,1] and 0 on diagonal)
//' @param eps Threshold for values to be interpreted as equal to 0 (default = 1e-10)
//' @return no return value
//'
//' @examples
//' L=matrix(rexp(9),nrow=3)
//' lambda <- matrix(0.5,nrow=3,ncol=3)
//' p <- matrix(0.7, nrow=3,ncol=3)
//' ERE_step_cycle(r=c(0,1),c=c(1,2),L=L,lambda=lambda,p=p)
//' ERE_step_cycle(r=c(0,1,2),c=c(0,1,2),L=L,lambda=lambda,p=p)
//' ERE_step_cycle(r=c(0,1,2),c=c(2,1,0),L=L,lambda=lambda,p=p)
//'
//' @export
// [[Rcpp::export]]
void ERE_step_cycle(std::vector<int> r, std::vector<int> c, NumericMatrix& L, NumericMatrix lambda, NumericMatrix p, double eps=1e-10){
  int k=r.size();
  if (c.size()!=r.size())
    throw Rcpp::exception("Length of r and c do not match");
  // check if r an c are within the allowed range
  for (int i=0; i<k; i++){
    if (r[i]<0||r[i]>=L.nrow()||c[i]<0||c[i]>=L.ncol())
      throw Rcpp::exception("r and c not in the allowed range");
  }

  // Check if there is at least one enforced 0 - in this case stop.
  for (int i=0; i<k; i++){
    if (p(r[i],c[i])==0){
      if (L(r[i],c[i])>eps)
	throw Rcpp::exception("Internal error: enforced 0 not present");
      return;
    }
    if (p(r[i],c[(i+1)%k])==0){
      if (L(r[i],c[(i+1)%k])>eps)
	throw Rcpp::exception("Internal error: enforced 0 not present");
      return;
    }
  }
  // compute bound;
  double deltabound[2];
  deltabound[0]=-L(r[0],c[0]);
  deltabound[1]=L(r[0],c[1]);
  for (int i=1; i<k; i++){
    deltabound[0]=max(deltabound[0],-L(r[i],c[i]));
    deltabound[1]=min(deltabound[1],L(r[i],c[(i+1)%k]));
  }
  if (std::isnan(deltabound[0])||std::isnan(deltabound[1])){
    Rcpp::Rcout<<" " <<deltabound[0]<<" " <<deltabound[1]<<std::endl;
    throw Rcpp::exception("Internal error: upper or lower bound is nan");
  }
  double len = deltabound[1]-deltabound[0];
  if (len<eps){
    if (fabs(deltabound[0])>eps)
      throw Rcpp::exception("Unexpected state: upper bound for Delta equals lower bound for Delta should imply they are equal to 0");
    return; // nothing needs doing, Delta=0 will be used.
  }
  int nzeros[2];
  for (int j=0; j<2; j++){
    nzeros[j]=0;
    for (int i=0; i<k; i++){
      if (L(r[i],c[i])+deltabound[j]<eps)
	nzeros[j]++;
      if (L(r[i],c[(i+1)%k])-deltabound[j]<eps)
	nzeros[j]++;
    }
  }
  if (nzeros[0]<=1&&nzeros[1]<=1){
    if (nzeros[0]<1||nzeros[1]<1){
      Rcpp::Rcout<<nzeros[0]<<" "<<nzeros[1]<<" " <<deltabound[0]<<" " <<deltabound[1]<<std::endl;
      throw Rcpp::exception("Internal error: should be able to set at least one element to 0");
    }
    double lambdamarg=0;
    for (int i=0; i<k; i++){
      lambdamarg+=lambda(r[i],c[i])-lambda(r[i],c[(i+1)%k]);
    }
    lambdamarg=roundnum(lambdamarg,eps);
    double pall= loglDelta(L,lambda,p,r,c,(deltabound[0]+deltabound[1])/2.0,eps); // intermediate point to ensure no value of 0 at any position
    if (lambdamarg==0.)
      pall+=log(len);
    else
      pall+=log(-1./lambdamarg*(exp(-lambdamarg*(len/2.))-exp(-lambdamarg*(-len/2.))));
    double lh[2];
    for (int i=0; i<2; i++){
      lh[i]=loglDelta(L,lambda,p,r,c,deltabound[i],eps);
    }
    // correction for log-scale
    double maxval = pall;
    if (lh[0]>maxval) maxval=lh[0];
    if (lh[1]>maxval) maxval=lh[1];
    pall=exp(pall-maxval);
    lh[0]=exp(lh[0]-maxval);
    lh[1]=exp(lh[1]-maxval);
    double u=unif_rand()*(pall+lh[0]+lh[1]);
    if (pall+lh[0]+lh[1]==0.)
      throw Rcpp::exception("No case has positive likelihood.");
    if (pall>=u){
      double Delta;
      if (lambdamarg==0.){
	Delta = unif_rand()*len+deltabound[0];
      }else{
	Delta = -log(exp(-lambdamarg*deltabound[0])
		 -unif_rand()*(exp(-lambdamarg*deltabound[0])
			       -exp(-lambdamarg*deltabound[1]))
		 )/lambdamarg;
      }
      setresDelta(L,r,c,Delta,eps);
      return;
    }
    u-=pall;
    if (u<=lh[0]){
      setresDelta(L,r,c,deltabound[0],eps);
      return;
    }else{
      setresDelta(L,r,c,deltabound[1],eps);
      return;
    }
  }else{ // at least one case with two zeros...
    if (nzeros[0]>nzeros[1]){
      setresDelta(L,r,c,deltabound[0],eps);
      return;
    }
    if (nzeros[1]>nzeros[0]){
      setresDelta(L,r,c,deltabound[1],eps);
      return;
    }
    // tied number of large zeros, need to choose randomly according
    // to likelihood to decide
    double lh[2];
    for (int i=0; i<2; i++){
      lh[i]=loglDelta(L,lambda,p,r,c, deltabound[i],eps);
    }
    double maxval2=lh[0];
    if (lh[1]>maxval2) maxval2=lh[1];
    for (int i=0; i<2; i++){
      lh[i]=exp(lh[i]-maxval2);
    }
    double u=unif_rand()*(lh[0]+lh[1]);
    if (u<=lh[0]){
      setresDelta(L,r,c,deltabound[0],eps);
      return;
    }else{
      setresDelta(L,r,c,deltabound[1],eps);
      return;
    }
  }
}



//' Gibbs sampling step of a matrix in the ERE model
//'
//' The sampling is conditional on row and column sums and uses k-cycle steps. Then dimensions of L, lambda and p must match.
//'
//' @param L Starting matrix - will be modified to contain the results.
//' @param lambda  Matrix of intensities
//' @param p Matrix of probabilities (must be in [0,1])
//' @param it Number of iterations (default=1000)
//' @param eps Threshold for values to be interpreted as equal to 0 (default = 1e-10)
//' @param debug Should addtional debug information be printed? (0 no output, 1 output debug information)
//' @return no return value
//'
//' @examples
//' L <- matrix(c(1,2,3,4,5,6,7,8,9),nrow=3)
//' diag(L) <- 0
//' lambda <- matrix(0.5,nrow=3,ncol=3)
//' p <- matrix(0.7, nrow=3,ncol=3)
//' diag(p) <- 0
//' GibbsSteps_kcycle(L=L,lambda=lambda,p=p)
//' L
//' L <- matrix(1:16,nrow=4)
//' diag(L) <- 0
//' lambda <- matrix(0.5,nrow=4,ncol=4)
//' p <- matrix(0.25, nrow=4,ncol=4)
//' diag(p) <- 0
//' GibbsSteps_kcycle(L=L,lambda=lambda,p=p)
//' L
//' @export
// [[Rcpp::export]]
void GibbsSteps_kcycle(NumericMatrix& L, NumericMatrix lambda, NumericMatrix p,int it=1000, double eps=1e-10, int debug=0) {
  int nrow=L.nrow();
  int ncol=L.ncol();
  if (nrow<=2||ncol<=2)
    throw Rcpp::exception("L must have at least 2 rows and columns");
  if (ncol!=lambda.ncol()||nrow!=lambda.nrow())
    throw Rcpp::exception("Dimensions of lambda and L do not match.");
  if (ncol!=p.ncol()||nrow!=p.nrow())
    throw Rcpp::exception("Dimensions of p and L do not match.");

  int kmax=min(nrow,ncol);
  std::vector<double> weights(kmax);
  weights[0]=0.;
  weights[1]=1000.;// arbitrary starting number
  for (int i=2;i<kmax; i++){
    weights[i]=weights[i-1]/2.;
  }
  double weightstot=weights[1]*2.-1.;
  for (int i=kmax-1;i>=0; i--){
    weights[i]/=weightstot;
  }
  if (debug){
    Rcpp::Rcout<<"Weights for choosing cycle size: ";
    for (int i=0; i<kmax; i++)
      Rcpp::Rcout<<weights[i]<<" ";
    Rcpp::Rcout<<std::endl;
  }

  std::vector<int> rall(nrow), call(ncol);
  for (int i=0;i<nrow; i++){
    rall[i]=i;
  }
  for (int i=0;i<ncol; i++){
    call[i]=i;
  }

  std::vector<int>::iterator rit, cit;
  for (int step=0;step<it; step++){
    double ksel = unif_rand();
    int k=0;
    rit=rall.begin();
    cit=call.begin();
    //no movement if at least one zero picked in the odd and even parts of the cycle
    int nzeros[2];
    nzeros[1]=0;
    bool discard=false;
    std::swap(*cit,*(cit+randWrapper(ncol)));
    std::swap(*rit,*(rit+randWrapper(nrow)));
    nzeros[0]=(L(*rit,*cit)==0.);
    for (k=1; k<kmax; k++){
      cit++;
      std::swap(*cit,*(cit+randWrapper(ncol-k)));
      if (L(*rit,*cit)==0.){
	if (nzeros[0]>0){
	  discard=true;
	  break;
	} else {
	  nzeros[1]++;
	}
      }
      rit++;
      std::swap(*rit,*(rit+randWrapper(nrow-k)));
      if (L(*rit,*cit)==0.){
	if (nzeros[1]>0){
	  discard=true;
	  break;
	} else {
	  nzeros[0]++;
	}
      }
      ksel-=weights[k];
      if (ksel<0) break;
    }
    if (!discard){
      if (L(*rit,call[0])>0.||nzeros[0]==0){//check last point
	cit++;
	rit++;
	std::vector<int> rsub(rall.begin(),rit);
	std::vector<int> csub(call.begin(),cit);
	if (debug){
	  Rcpp::Rcout<<"Row indices:";
	  for (unsigned int i=0; i<rsub.size(); i++)
	    Rcpp::Rcout<<rsub[i]<<",";
	  Rcpp::Rcout<<" ";
	  Rcpp::Rcout<<"Column indices:";
	  for (unsigned int i=0; i<rsub.size(); i++)
	    Rcpp::Rcout<<csub[i]<<",";
	  Rcpp::Rcout<<std::endl;
	}
	ERE_step_cycle(rsub,csub,L,lambda,p,eps);
      }
    }
  }
  return;
}

