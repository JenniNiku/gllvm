#define TMB_LIB_INIT R_init_gllvm
#include <TMB.hpp>
#include<math.h>
//--------------------------------------------------------
//GLLVM
//Author: Jenni Niku
//------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  //declares all data and parameters used
  DATA_MATRIX(y);
  DATA_MATRIX(x);
  DATA_MATRIX(xr);
  DATA_MATRIX(offset);

  PARAMETER_MATRIX(r0);
  PARAMETER_MATRIX(b);
  PARAMETER_MATRIX(B);
  PARAMETER_VECTOR(lambda);

  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi);
  PARAMETER(log_sigma);// log(SD for row effect)

  DATA_INTEGER(num_lv);
  DATA_INTEGER(family);
  
  PARAMETER_VECTOR(Au);
  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_VECTOR(zeta);
  
  DATA_SCALAR(extra);
  DATA_INTEGER(method);// 0=VA, 1=LA
  DATA_INTEGER(model);
  DATA_VECTOR(random);//random row

  int n = y.rows();
  int p = y.cols();
  vector<Type> iphi = exp(lg_phi);
  vector<Type> Ar = exp(lg_Ar);
  Type sigma = exp(log_sigma);

  if(random(0)<1){  r0(0,0) = 0;}

  matrix<Type> eta(n,p);
  eta.fill(0.0);
  matrix<Type> lam(n,p);
  lam.fill(0.0);
  
  matrix<Type> newlam(num_lv,p);
  if(num_lv>0){

    for (int j=0; j<p; j++){
      for (int i=0; i<num_lv; i++){
        if (j < i){
          newlam(i,j) = 0;
        } else{
          newlam(i,j) = lambda(j);
          if (i > 0){
            newlam(i,j) = lambda(i+j+i*p-(i*(i-1))/2-2*i);
          }
        }
        // set diag>0 !!!!!!!!!!!
        /*if (j == i){
          newlam(i,j) = exp(newlam(i,j));
        }*/
      }
    }

    //To create lambda as matrix upper triangle

    lam += u*newlam;
    eta = lam;
  }

  matrix<Type> mu(n,p);

  Type nll = 0.0; // initial value of log-likelihood


  if(method<1){
    eta += r0*xr + offset;

    matrix<Type> cQ(n,p);
    cQ.fill(0.0);
    
    if(random(0)>0){
      for (int i=0; i<n; i++) {
        Ar(i)=pow(Ar(i),2);
        for (int j=0; j<p;j++){
          cQ(i,j) = 0.5* Ar(i);//
        }
      }}

    if(num_lv>0){
      array<Type> A(num_lv,num_lv,n);
      for (int d=0; d<(num_lv); d++){
        for(int i=0; i<n; i++){
          A(d,d,i)=exp(Au(d*n+i));
        }
      }
      if(Au.size()>num_lv*n){
        int k=0;
        for (int c=0; c<(num_lv); c++){
          for (int r=c+1; r<(num_lv); r++){
            for(int i=0; i<n; i++){
              A(r,c,i)=Au(num_lv*n+k*n+i);
              A(c,r,i)=A(r,c,i);
            }
            k++;
          }}
      }
      /*Calculates the commonly used (1/2) theta'_j A_i theta_j
      A is a num.lv x nmu.lv x n array, theta is p x num.lv matrix*/
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          cQ(i,j) += 0.5*((newlam.col(j)).transpose()*(A.col(i).matrix()*newlam.col(j))).sum();
        }
        nll -= 0.5*(log(A.col(i).matrix().determinant()) - A.col(i).matrix().diagonal().sum());// log(det(A_i))-sum(trace(A_i))*0.5 sum.diag(A)
      }

    }

    if(model<1){
      eta += x*b;
    } else {
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)+eta1(m,0);
          m++;
        }
      }
    }
    //likelihood
    if(family==0){
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll -= dpois(y(i,j), exp(eta(i,j)+cQ(i,j)), true)-y(i,j)*cQ(i,j);
        }
        nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==1){
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll -= y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*cQ(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
        }
        nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==2) {
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
          nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j)))) - cQ(i,j);
        }
        nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==6){
      int ymax =  CppAD::Integer(y.maxCoeff());
      int K = ymax - 1;
      
      matrix <Type> zetanew(p,K);
      zetanew.fill(0.0);
      
      int idx = 0;
      for(int j=0; j<p; j++){
        int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
        int Kj = ymaxj - 1;
        if(Kj>1){
          for(int k=0; k<(Kj-1); k++){
            if(k==1){
              zetanew(j,k+1) = abs(zeta(idx+k));//second cutoffs must be positive
            }else{
              zetanew(j,k+1) = zeta(idx+k);
            }
            
          }
        }
        idx += Kj-1;
      }

      for (int i=0; i<n; i++) {
        for(int j=0; j<p; j++){
          int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
          //minimum category
          if(y(i,j)==1){
            nll -= log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)));
          }else if(y(i,j)==ymaxj){
            //maximum category
            int idx = ymaxj-2;
            nll -= log(1 - pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1)));
          }else if(ymaxj>2){
            for (int l=2; l<ymaxj; l++) {
              if(y(i,j)==l && l != ymaxj){
                nll -= log(pnorm(zetanew(j,l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,l-2)-eta(i,j), Type(0), Type(1))); 
              }
            }
          }
        
          nll += cQ(i,j);
          //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));// 
        }
        nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
     }else {
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j)));
        }
        nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    }
    nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i

  } else {
    eta += r0*xr + offset;

    if(model<1){

      eta += x*b;
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          mu(i,j) = exp(eta(i,j));
        }
      }
    } else {

      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)+eta1(m,0);
          m++;
          mu(i,j) = exp(eta(i,j));
        }
      }
    }
    //latent variable is assumed to be from N(0,1)
    if(num_lv>0){
      for (int j=0; j<u.cols();j++){
        for (int i=0; i<n; i++) {
          nll -= dnorm(u(i,j),Type(0),Type(1),true);
        }
      }}
    for(int i = 0; i < n; i++){
      nll -= dnorm(r0(i,0), Type(0), sigma, true)*random(0);
    }

    //likelihood model with the log link function
    if(family<1){
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          nll -= dpois(y(i,j), exp(eta(i,j)), true);
        }
      }
    } else if(family<2){
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          nll -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
        }
      }} else if(family<3) {
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(extra<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));
          }
        }
      } else if(family<4){
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dnorm(y(i,j), eta(i,j), iphi(j), true); //gamma family
          }
        }
      } else if(family<5){
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dtweedie(y(i,j), exp(eta(i,j)),iphi(j),extra, true); //tweedie family
          }
        }
      } else {
        iphi=iphi/(1+iphi);
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dzipois(y(i,j), exp(eta(i,j)),iphi(j), true); //zero-infl-poisson
          }
        }}
  }
  return nll;
}
