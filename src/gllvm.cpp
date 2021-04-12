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
  DATA_MATRIX(x_lv);
  DATA_MATRIX(xr);
  DATA_MATRIX(xb);
  DATA_MATRIX(offset);
  
  PARAMETER_MATRIX(r0);
  PARAMETER_MATRIX(b);
  PARAMETER_MATRIX(B);
  PARAMETER_MATRIX(Br);
  PARAMETER_MATRIX(b_lv);
  PARAMETER_VECTOR(lambda);
  PARAMETER_MATRIX(lambda2);
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi);
  PARAMETER_VECTOR(sigmaB);
  PARAMETER_VECTOR(sigmaij);
  PARAMETER_VECTOR(log_sigma);// log(SD for row effect)
  
  DATA_INTEGER(num_lv);
  DATA_INTEGER(num_lv_c);
  DATA_INTEGER(family);
  DATA_INTEGER(quadratic);
  
  PARAMETER_VECTOR(Au);
  //  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_VECTOR(Abb);
  PARAMETER_VECTOR(zeta);
  
  DATA_VECTOR(extra);
  DATA_INTEGER(method);// 0=VA, 1=LA
  DATA_INTEGER(model);
  DATA_VECTOR(random);//random row
  DATA_INTEGER(zetastruc);
  
  int n = y.rows();
  int p = y.cols();
  int l = xb.cols();
  vector<Type> iphi = exp(lg_phi);
  //vector<Type> Ar = exp(lg_Ar);
  Type sigma = exp(log_sigma(0));
  
  if(random(0)<1){  r0(0,0) = 0;}
  int nlvr = num_lv+num_lv_c;
  if(random(0)>0){  nlvr++;}
  
  matrix<Type> eta(n,p);
  eta.fill(0.0);
  matrix<Type> lam(n,p);
  lam.fill(0.0);
  matrix<Type> Cu(nlvr,nlvr); 
  Cu.fill(0.0);
  
  matrix<Type> newlam(nlvr,p);

  if(nlvr>0){
    newlam.row(0).fill(1.0);
    Cu.diagonal().fill(1.0);
    if(random(0)>0){
      Cu(0,0) = sigma*sigma;
      if(log_sigma.size()>1){
        for (int d=1; d<(nlvr); d++){
          Cu(d,0) = log_sigma(d);
          Cu(0,d) = Cu(d,0);
        }
      }
    }
    
    //To create lambda as matrix upper triangle
    if (num_lv>0){
      int tri = 0;
      if(num_lv_c>0){
        //because the lambdas for constrained and unconstrained LVs are separately identifiable
        tri += num_lv_c-1+p+(num_lv_c-1)*p-((num_lv_c-1)*(num_lv_c-1-1))/2-2*(num_lv_c-1); //number of elements for num_lv
      }
      for (int j=0; j<p; j++){
        for (int i=0; i<num_lv; i++){
          if (j < i){
            newlam(i+nlvr-num_lv,j) = 0;
          } else{
            newlam(i+nlvr-num_lv,j) = lambda(j+tri);
            if (i > 0){
              newlam(i+nlvr-num_lv,j) = lambda(i+j+i*p-(i*(i-1))/2-2*i+tri);
            }
          }
          // set diag>0 !!!!!!!!!!!
          // if (j == i){
          //   newlam(i+nlvr-num_lv,j) = exp(newlam(i+nlvr-num_lv,j));
          // }
        }
      }
    }
    if (num_lv_c>0){
      
      for (int j=0; j<p; j++){
        for (int i=0; i<num_lv_c; i++){
          if (j < i){
            newlam(i+nlvr-num_lv-num_lv_c,j) = 0;
          } else{
            newlam(i+nlvr-num_lv-num_lv_c,j) = lambda(j);
            if (i > 0){
              newlam(i+nlvr-num_lv-num_lv_c,j) = lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }
          // set diag>0 !!!!!!!!!!!
          // if (j == i){
          //   newlam(i+nlvr-num_lv,j) = exp(newlam(i+nlvr-num_lv,j));
          // }
        }
      }
    }
    
}
  
  matrix<Type> mu(n,p);
  mu.fill(0.0);
  
  using namespace density;
  
  matrix <Type> nll(n,p); // initial value of log-likelihood
  nll.fill(0.0);
  
  if(method<1){
    //quadratic coefficients
    //if random rows, add quadratic coefficients to q>0
    array<Type> D(nlvr,nlvr,p);
    D.fill(0.0);
    if(quadratic>0){
      if(nlvr>(num_lv+num_lv_c)){
        if(lambda2.cols()==1){
          for (int j=0; j<p; j++){
            for (int q=1; q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-1,0)); //common tolerances model
            }
          } 
        }else{
          for (int j=0; j<p; j++){
            for (int q=1; q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-1,j)); //full quadratic model
            }
          } 
        }
        
      }else{
        if(lambda2.cols()==1){
          for (int j=0; j<p; j++){
            for (int q=0; q<(num_lv+num_lv_c); q++){
              D(q,q,j) = fabs(lambda2(q,0)); //common tolerances model
            }
          } 
        }else{
          for (int j=0; j<p; j++){
            for (int q=0; q<(num_lv+num_lv_c); q++){
              D(q,q,j) = fabs(lambda2(q,j)); //full quadratic model
            }
          } 
        }
      }
    }
    
    eta += r0*xr + offset;
    
    matrix<Type> cQ(n,p);
    cQ.fill(0.0);
    array<Type> A(nlvr,nlvr,n);
    A.fill(0.0);
    
    if(nlvr>0){// log-Cholesky parametrization for A_i:s
      for (int d=0; d<(nlvr); d++){
        for(int i=0; i<n; i++){
          A(d,d,i)=exp(Au(d*n+i));
        }
      }
      if(Au.size()>nlvr*n){
        int k=0;
        for (int c=0; c<(nlvr); c++){
          for (int r=c+1; r<(nlvr); r++){
            for(int i=0; i<n; i++){
              A(r,c,i)=Au(nlvr*n+k*n+i);
              // A(c,r,i)=A(r,c,i);
            }
            k++;
          }}
      }
      //set VA covariances for random rows to zero for quadratic model
      if(quadratic>0&&nlvr>(num_lv+num_lv_c)){
        for(int i=0; i<n; i++){
          for (int d=0; d<(nlvr); d++){
            if(d!=0){
              A(d,0,i) = 0.0;
            }
          }
        }
      }
      for(int i=0; i<n; i++){
        if(nlvr == (num_lv+num_lv_c)) nll.row(i).array() -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()+(u.row(i)*u.row(i).transpose()).sum()))/p;
         if(nlvr>(num_lv+num_lv_c)) nll.row(i).array() -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(Cu.inverse()*(A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()-0.5*((u.row(i)*Cu.inverse())*u.row(i).transpose()).sum())/p;
        // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
      }
      nll.array() -= -0.5*log(Cu.determinant())*random(0)/p;//n*
    }
    
    
    // Include random slopes if random(1)>0
    if(random(1)>0){
      matrix<Type> sds(l,l);
      sds.fill(0.0);
      sds.diagonal() = exp(sigmaB);
      matrix<Type> S=sds*UNSTRUCTURED_CORR(sigmaij).cov()*sds;
      
      // log-Cholesky parametrization for A_bj:s
      array<Type> Ab(l,l,p);
      Ab.fill(0.0);
      for (int dl=0; dl<(l); dl++){
        for(int j=0; j<p; j++){
          Ab(dl,dl,j)=exp(Abb(dl*p+j));
        }
      }
      if(Abb.size()>l*p){
        int k=0;
        for (int c=0; c<(l); c++){
          for (int r=c+1; r<(l); r++){
            for(int j=0; j<p; j++){
              Ab(r,c,j)=Abb(l*p+k*p+j);
              // Ab(c,r,j)=Ab(r,c,j);
            }
            k++;
          }}
      }
      
      /*Calculates the commonly used (1/2) x'_i A_bj x_i
       A is a num.lv x nmu.lv x n array, theta is p x num.lv matrix*/
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          cQ(i,j) += 0.5*((xb.row(i))*((Ab.col(j).matrix()*Ab.col(j).matrix().transpose()).matrix()*xb.row(i).transpose())).sum();
        }
        nll.col(j).array() -= ((((vector <Type> (Ab.col(j).matrix().diagonal())).log()).sum() - 0.5*(S.inverse()*(Ab.col(j).matrix()*Ab.col(j).matrix().transpose()).matrix()).trace()-0.5*(Br.col(j).transpose()*(S.inverse()*Br.col(j))).sum()))/n;// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
      }
      eta += xb*Br;
      nll.array() -= -0.5*log(S.determinant())*random(1)/n;//n*
    }
    
    if(model<1){
      eta += x*b;
    } else {
      // Fourth corner model
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)*extra(1)+eta1(m,0); //extra(1)=0 if beta0comm=TRUE
          m++;
        }
      }
    }
 if(num_lv_c>0){
      if(random(0)>0){
        //b_lv has rows K, colummns q
        //so, first column is random intercept,
        //second to num_lv_c are constrained
        //rest lvs are unconstrained
        //this needs to be done after the integration of (just) the std. normal
        for (int q=0; q<num_lv_c; q++){
          u.col(q+1) += x_lv*b_lv.col(q);
        }
      }else{
        for (int q=0; q<num_lv_c; q++){
          u.col(q) += x_lv*b_lv.col(q);
        }
        }
    }
      if(nlvr>0){
        lam += u*newlam;
      }
        
    if(quadratic < 1 && nlvr > 0){
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          cQ(i,j) += 0.5*((newlam.col(j)).transpose()*((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()*newlam.col(j))).sum();
        }
      }
      eta += lam;
    }
    matrix <Type> e_eta(n,p);
    e_eta.fill(0.0);
    if(quadratic>0 && nlvr > 0){
      matrix <Type> Acov(nlvr,nlvr);
      //quadratic model approximation
      //Poisson
      if(family==0){
        matrix <Type> B(nlvr,nlvr);
        matrix <Type> v(nlvr,1);
        for (int i=0; i<n; i++) {
          Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
          matrix <Type> Q = atomic::matinv(Acov);
          for (int j=0; j<p;j++){
            B = (2*D.col(j).matrix()+Q);
            v = (newlam.col(j)+Q*u.row(i).transpose());
            Type detB = atomic::logdet(B);
            Type detA = ((vector <Type> (A.col(i).matrix().diagonal())).log()).sum(); //log-determinant of cholesky
            e_eta(i,j) += exp(cQ(i,j) + eta(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()-detB)-detA); //add all the other stuff to the quadratic approximation
            eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
          }
        }
      }
      // //NB, gamma, exponential
      if(family==1|family==4|family==8){
        matrix <Type> B(nlvr,nlvr);
        matrix <Type> v(nlvr,1);
        for (int i=0; i<n; i++) {
          Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
          matrix <Type> Q = atomic::matinv(Acov);
          for (int j=0; j<p;j++){
            B = (-2*D.col(j).matrix()+Q);
            v = (-newlam.col(j)+Q*u.row(i).transpose());
            Type detB = log((B.llt().matrixL()).determinant());//required like this due to potential negative semi-definiteness
            Type detA = ((vector <Type> (A.col(i).matrix().diagonal())).log()).sum();
            e_eta(i,j) += exp(-eta(i,j) - cQ(i,j)+0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value())-detA-detB);
            eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
            
          }
        }
      }
      //Binomial, Gaussian, Ordinal
      if(family==2|family==3|family==7){
        for (int i=0; i<n; i++) {
          Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
          for (int j=0; j<p;j++){
            cQ(i,j) += 0.5*(newlam.col(j)*newlam.col(j).transpose()*Acov).trace() + (D.col(j).matrix()*Acov*D.col(j).matrix()*Acov).trace() +2*(u.row(i)*D.col(j).matrix()*Acov*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*Acov*newlam.col(j)).value();
            eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
          }
        }
      }
    }
    
    if(family==0){//poisson
      if(quadratic<1){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= dpois(y(i,j), exp(eta(i,j)+cQ(i,j)), true)-y(i,j)*cQ(i,j);
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= y(i,j)*eta(i,j) - e_eta(i,j) - lfactorial(y(i,j));
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }
    } else if(family==1){//NB
      if(quadratic<1){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*cQ(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }  
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= -iphi(j)*eta(i,j) -(y(i,j)+iphi(j))*log(1+iphi(j)*e_eta(i,j))+ lgamma(y(i,j)+iphi(j))+ iphi(j)*log(iphi(j)) -lgamma(iphi(j)) -lfactorial(y(i,j));
            //log(1+phi*e_eta) = log(phi+1/e_eta)+log(e_eta)
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }
      
    } else if(family==2) {//binomial probit
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
          nll(i,j) -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j)))) - cQ(i,j);
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==3) {//gaussian
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll(i,j) -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j)));
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==4) {//gamma
      if(quadratic<1){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));          
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -=  ( -eta(i,j) - e_eta(i,j)*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }
      
    } else if(family==7 && zetastruc == 1){//ordinal
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
              zetanew(j,k+1) = fabs(zeta(idx+k));//second cutoffs must be positive
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
            nll(i,j) -= log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)));
          }else if(y(i,j)==ymaxj){
            //maximum category
            int idx = ymaxj-2;
            nll(i,j) -= log(1 - pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1)));
          }else if(ymaxj>2){
            for (int l=2; l<ymaxj; l++) {
              if(y(i,j)==l && l != ymaxj){
                nll(i,j) -= log(pnorm(zetanew(j,l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,l-2)-eta(i,j), Type(0), Type(1))); 
              }
            }
          }
          
          nll(i,j) += cQ(i,j);
          //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));// 
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    }else if(family==7 && zetastruc==0){
      int ymax =  CppAD::Integer(y.maxCoeff());
      int K = ymax - 1;
      
      vector <Type> zetanew(K);
      zetanew.fill(0.0);
      for(int k=0; k<(K-1); k++){
        if(k==1){
          zetanew(k+1) = fabs(zeta(k));//second cutoffs must be positive
        }else{
          zetanew(k+1) = zeta(k);
        }
      }
      for (int i=0; i<n; i++) {
        for(int j=0; j<p; j++){
          //minimum category
          if(y(i,j)==1){
            nll(i,j) -= log(pnorm(zetanew(0) - eta(i,j), Type(0), Type(1)));
          }else if(y(i,j)==ymax){
            //maximum category
            int idx = ymax-2;
            nll(i,j) -= log(1 - pnorm(zetanew(idx) - eta(i,j), Type(0), Type(1)));
          }else if(ymax>2){
            for (int l=2; l<ymax; l++) {
              if(y(i,j)==l && l != ymax){
                nll(i,j) -= log(pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1)));
              }
            }
          }
          nll(i,j) += cQ(i,j);
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==8) {// exp dist
      if(quadratic<1){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) );
          }
        }  
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= ( -eta(i,j) - e_eta(i,j)*y(i,j) );
          }
        }
      }
      
    }
    // nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
    
  } else {

    // Include random slopes if random(1)>0
    if(random(1)>0){
      vector<Type> sdsv = exp(sigmaB);
      density::UNSTRUCTURED_CORR_t<Type> neg_log_MVN(sigmaij);
      for (int j=0; j<p;j++){
        nll.col(j).array() += VECSCALE(neg_log_MVN,sdsv)(vector<Type>(Br.col(j)))/n;
      }
      eta += xb*Br;
    }
    
    //latent variables and random site effects (r_i,u_i) from N(0,Cu)
    if(nlvr>0){
      
      MVNORM_t<Type> mvnorm(Cu);
      for (int i=0; i<n; i++) {
        nll.row(i).array() += mvnorm(u.row(i))/p;
      }
    }
    if(num_lv_c>0){
      if(random(0)>0){
        //b_lv has rows K, colummns d
        //so, first column is random intercept,
        //second to num_lv_c are constrained
        //rest lvs are unconstrained
        //this needs to be done after the integration of (just) the std. normal
        for (int q=0; q<num_lv_c; q++){
          u.col(q+1) += x_lv*b_lv.col(q);
        }
      }else{
        for (int q=0; q<num_lv_c; q++){
          u.col(q) += x_lv*b_lv.col(q);
        }
      }
    }
    if(nlvr>0){
      lam += u*newlam;
    }
    eta += r0*xr + offset;
    if(nlvr>0){
      eta += lam;
    }
    
    
    if(model<1){
      
      eta += x*b;
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          mu(i,j) = exp(eta(i,j));
        }
      }
      
    } else {
      // Fourth corner model
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)*extra(1)+eta1(m,0);
          m++;
          mu(i,j) = exp(eta(i,j));
        }
      }
    }
    
    
    //likelihood model with the log link function
    if(family==0){//poisson family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          nll(i,j) -= dpois(y(i,j), exp(eta(i,j)), true);
        }
      }
    } else if(family==1){//negative.binomial family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          nll(i,j) -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
        }
      }} else if(family==2) {//binomial family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            nll(i,j) -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));
          }
        }
      } else if(family==3){//gaussian family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll(i,j) -= dnorm(y(i,j), eta(i,j), iphi(j), true); 
          }
        }
      } else if(family==4){//gamma family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll(i,j) -= dgamma(y(i,j), iphi(j), exp(eta(i,j))/iphi(j), true); 
          }
        }
      } else if(family==5){//tweedie family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll(i,j) -= dtweedie(y(i,j), exp(eta(i,j)),iphi(j),extra(0), true); 
          }
        }
      } else if(family==6) {//zero-infl-poisson
        iphi=iphi/(1+iphi);
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll(i,j) -= dzipois(y(i,j), exp(eta(i,j)),iphi(j), true); 
          }
        }
      } else if(family==8) {// exponential family
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll(i,j) -= dexp(y(i,j), exp(-eta(i,j)), true);  // (-eta(i,j) - exp(-eta(i,j))*y(i,j) );
          }
        }
      } else if(family==9) {// beta family
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            nll(i,j) -= dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
          }
        }
      }
  }
  REPORT(nll);//only works for VA!!
  REPORT(newlam);
  REPORT(b_lv);
  return nll.sum();
}
