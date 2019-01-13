#include<Rcpp.h>
#include <queue>

using namespace Rcpp;
using std::vector;
typedef std::pair<int,int> PII;


//' Finds a Nonnegative Matrix Satisfying Row and Column Sums
//'
//' Given row and column sums and a matrix p which indicates which elements of the matrix can be present, this function computes a nonnegative matrix that match these row and column sums. If this is not possible then the function returns an error message.
//'
//' The function transforms the problem into a Maximum Flow
//' problem of a graph and uses the Edmonds-Karps algorithm to solve it.
//' If the error message "Could not find feasible matrix." is produced then this could be
//' due to p imposing disconnected components in the graph implied
//' by row and column sums that are not compatible with the row and column sums..
//'
//' @param r vector of row sums (nonnegative
//' @param c vector of column sums (nonnegative)
//' @param p matrix of probabilities (must be in [0,1]), matching the dimensions of r and c. Values of p=0 are interpreted that the corresponding matrix elements have to be 0. Note: p=1 does not force the corresponding matrix element to exist.
//' @param eps row and col sums can at most be different by eps.  Default 1e-9.
//' @return A feasible matrix.
//'
//' @examples
//' p=matrix(c(1,0,0,1),nrow=2)
//' findFeasibleMatrix(c(1,1),c(1,1),p=p)
//'
//' n <- 4
//' M <- matrix(nrow=n,ncol=n,rexp(n*n)*(runif(n*n)>0.6))
//' M
//' r <- rowSums(M)
//' c <- colSums(M)
//' Mnew <- findFeasibleMatrix(r=r,c=c,p=(M>0)*0.5)
//' Mnew
//' rowSums(M);rowSums(Mnew)
//' colSums(M);colSums(Mnew)
//' @export
// [[Rcpp::export]]
NumericMatrix findFeasibleMatrix(std::vector<double> r, std::vector<double> c, NumericMatrix p, double eps=1e-9){
  if (r.size()!=(unsigned)p.nrow())
    throw Rcpp::exception("Size of r and number of rows of p does not match.");
  if (c.size()!=(unsigned)p.ncol())
    throw Rcpp::exception("Size of c and number of rows of p does not match.");
  // check if r an c are nonnegative
  for (unsigned int i=0; i<r.size(); i++)
    if (r[i]<0) throw Rcpp::exception("r must have nonnegative entries.");
  for (unsigned int i=0; i<c.size(); i++)
    if (c[i]<0) throw Rcpp::exception("r must have nonnegative entries.");
  // check that sum of r and sum of c match
  double rtot=std::accumulate(r.begin(), r.end(), 0.0);
  double ctot=std::accumulate(c.begin(),c.end(),0.0);
  if (fabs(rtot-ctot)>eps)
    throw Rcpp::exception("Sums of r and c differ by more than eps.");

  // Build Graph
  vector<vector<int> > E;
  int source=0;
  int sink=r.size()+c.size()+1;
  int nV=r.size()+c.size()+2; // nubmer of vertices
  // 0=source,  then rows, 1..n cols, 2n+1=sink
  vector<vector<int> > e(nV);//adjacency list
  vector<vector<double> > capacity(nV,vector<double>(nV,0));//capacities
  vector<vector<double> > flow(nV,vector<double>(nV,0));//flow
  for (unsigned i=0; i<r.size(); i++){
    e[source].push_back(i+1);
    capacity[source].push_back(r[i]);
    e[i+1].push_back(source);
    capacity[source][i+1]=r[i];
  }
  for (unsigned i=0; i<c.size(); i++){
    e[sink].push_back(r.size()+i+1);
    e[r.size()+i+1].push_back(sink);
    capacity[r.size()+i+1][sink]=c[i];
  }
  for (unsigned i=0; i<r.size(); i++){
    for (unsigned j=0; j<c.size(); j++){
      if (p(i,j)>0){
	unsigned a=i+1;
	unsigned b=r.size()+j+1;
	e[a].push_back(b);
	e[b].push_back(a);
	capacity[a][b]=INFINITY;
      }
    }
  }

  vector<int> colour; // need this later in case of the network being not feasible
  while (true){
    //determine augmenting path.
    //breadth first search; Cormen p. 470
    int myInf=2*nV;
    colour=vector<int>(nV,0);
    vector<int> d(nV,myInf);
    vector<int> pred(nV,-1);
    colour[source]=1;
    d[source]=0;
    std::priority_queue<PII > Q;
    Q.push(PII(0,source));
    while (Q.size()>0&&colour[sink]==0){ // stop once sink has been found
      int u=Q.top().second;
      Q.pop();
      if (colour[u]==1){
	colour[u]=2;
	for (vector<int>::iterator v=e[u].begin(); v!=e[u].end(); v++){
	  if (colour[*v]==0){ // not reached yet
	    if (capacity[u][*v]>flow[u][*v]-flow[*v][u]+eps){ // check residual network
	      colour[*v]=1;
	      pred[*v]=u;
	      d[*v]=d[u]+1;
	      Q.push(PII(-d[*v],*v));
	      if (*v==sink){ // found the sink
		break;
	      }
	      if ((unsigned)*v>1+r.size()){ //if a column then check if sink can be reached
		if (capacity[*v][sink]>flow[*v][sink]-flow[sink][*v]+eps){
		  colour[sink]=1;
		  pred[sink]=*v;
		  d[sink]=d[*v]+1; // correct as graph is bipartite
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
    if (colour[sink]==0){
      break;  // no more augmenting path
    }else{  // found augmenting path
      //determine maximum flow along path
      int pos=sink;
      double flownew=capacity[pred[pos]][pos]-flow[pred[pos]][pos]+flow[pos][pred[pos]];
      pos=pred[pos];
      while(pos!=source){
	double fmax=capacity[pred[pos]][pos]-flow[pred[pos]][pos]+flow[pos][pred[pos]];
	if (fmax<flownew) flownew=fmax;
	pos=pred[pos];
      }
      // add flow along path
      pos=sink;
      while(pos!=source){
	flow[pred[pos]][pos]+=flownew;
	pos=pred[pos];
      }

    }
  }


  /// prepare output
  NumericMatrix res(r.size(),c.size());
  for (unsigned i=0; i<r.size(); i++){
    for (unsigned j=0; j<c.size(); j++){
      res(i,j)=flow[i+1][j+r.size()+1]-flow[j+r.size()+1][i+1];
    }
  }


  // now check if row and column sums are actually being met....
  vector<int> rwrong;
  vector<int> cwrong;
  for (unsigned i=0; i<r.size(); i++){
    double tot=0.;
    for (unsigned j=0; j<c.size(); j++){
      tot+=res(i,j);
    }
    if (fabs(tot-r[i])>eps)
      rwrong.push_back(i+1);  //+1 to adjust for how R indexes vectors
  }
  for (unsigned i=0; i<c.size(); i++){
    double tot=0.;
    for (unsigned j=0; j<r.size(); j++){
      tot+=res(j,i);
    }
    if (fabs(tot-c[i])>eps)
      cwrong.push_back(i+1); //+1 to adjust for how R indexes vectors
  }
  if (rwrong.size()>0||cwrong.size()>0){
    // To get the cut set need to check  the residual network graph
    // look at colour: all with colour >0 in one component, all with colour=0 in the other component.
    std::stringstream Message;
    Message<<"Could not find feasible matrix. Subrectangle in row(s)";
    for (unsigned i=0; i<r.size(); i++){
      if (colour[i+1]>0) Message <<" "<<i+1;
    }
    Message<<" and column(s)";
    for (unsigned i=0; i<c.size(); i++){
      if (colour[i+1+r.size()]==0) Message <<" "<<i+1;
    }
    Message<<" violates existence condition.";
    throw Rcpp::exception(Message.str().c_str());
  }

  return(res);
}


