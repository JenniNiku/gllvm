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
  DATA_MATRIX(y); // matrix of responses
  DATA_MATRIX(x); // matrix of covariates
  DATA_MATRIX(xr); 
  DATA_MATRIX(xb); // envs with random slopes
  DATA_ARRAY(dr0); // design matrix for rows, (times, n, nr)
  DATA_MATRIX(offset); //offset matrix
  
  PARAMETER_MATRIX(r0); // site/row effects
  PARAMETER_MATRIX(b); // matrix of species specific intercepts and coefs
  PARAMETER_MATRIX(B); // coefs of 4th corner model
  PARAMETER_MATRIX(Br); // random slopes for envs
  PARAMETER_VECTOR(lambda); // lv loadings
  PARAMETER_MATRIX(lambda2);
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi); // dispersion params/extra zero probs for ZIP
  PARAMETER_VECTOR(sigmaB); // sds for random slopes
  PARAMETER_VECTOR(sigmaij);// cov terms for random slopes covariance
  PARAMETER_VECTOR(log_sigma);// log(SD for row effect) and 
  
  DATA_INTEGER(num_lv); // number of lvs
  DATA_INTEGER(family); // family index
  DATA_INTEGER(quadratic); // quadratic model, 0=no, 1=yes
  
  PARAMETER_VECTOR(Au); // variational covariances for u
   PARAMETER_VECTOR(lg_Ar); // variational covariances for r0
  PARAMETER_VECTOR(Abb);  // variational covariances for Br
  PARAMETER_VECTOR(zeta); // ordinal family param
  
  DATA_VECTOR(extra); // extra values, power of 
  DATA_INTEGER(method);// 0=VA, 1=LA
  DATA_INTEGER(model);// which model, basic or 4th corner
  DATA_VECTOR(random);//1=random, 0=fixed row params
  DATA_INTEGER(zetastruc); //zeta param structure for ordinal model
  DATA_INTEGER(rstruc); //Type for random rows. default = 0, when same as u:s. If 1, dr0 defines the structure. If 2, Points within groups has covariance struct defined by cstruc
  DATA_INTEGER(times); //number of time points
  DATA_INTEGER(cstruc); //correlation structure for row.params, 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm
  DATA_VECTOR(dc); //coordinates for sites, used for exponentially decaying cov. struc
  
  matrix<Type> dr = dr0.matrix();
  // REPORT(dr);
  
  int n = y.rows();
  int p = y.cols();
  // int nt =n;
  int nr =n;
  if(rstruc>0){
    nr = dr.cols();
  }
  
  
  int l = xb.cols();
  vector<Type> iphi = exp(lg_phi);
  vector<Type> Ar = exp(lg_Ar);
  Type sigma = exp(log_sigma(0));
  
  if(random(0)<1){  r0(0,0) = 0;}
  int nlvr = num_lv;
  if((random(0)>0) & (n == nr)){
    nlvr++;
    
    if(num_lv>0){
      u.conservativeResize(u.rows(), u.cols()+1);
      for (int i=0; i<n; i++){
        for (int q=num_lv; q>0; q--){
          u(i,q) = u(i,q-1);
        }
      }
      
      for (int i=0; i<n; i++){
        u(i,0) = r0(i,0); //can easily be extended to multiple rn
      }
    } else {
      u = r0;
    }
  }
  
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
    if((nlvr>num_lv)){
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
      
      for (int j=0; j<p; j++){
        for (int i=0; i<num_lv; i++){
          if (j < i){
            newlam(i+nlvr-num_lv,j) = 0;
          } else{
            newlam(i+nlvr-num_lv,j) = lambda(j);
            if (i > 0){
              newlam(i+nlvr-num_lv,j) = lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }
          // set diag>0 !!!!!!!!!!!
          // if (j == i){
          //   newlam(i+nlvr-num_lv,j) = exp(newlam(i+nlvr-num_lv,j));
          // }
        }
      }
    }
    
    lam += u*newlam;
    REPORT(u);
    REPORT(Cu);
    REPORT(newlam);
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
      if(nlvr>num_lv){
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
            for (int q=0; q<num_lv; q++){
              D(q,q,j) = fabs(lambda2(q,0)); //common tolerances model
            }
          } 
        }else{
          for (int j=0; j<p; j++){
            for (int q=0; q<num_lv; q++){
              D(q,q,j) = fabs(lambda2(q,j)); //full quadratic model
            }
          } 
        }
      }
    }
    REPORT(D);
    
    eta += offset;
    if((random(0)==0)){
      eta += r0*xr;
    }
    
    matrix<Type> cQ(n,p);
    cQ.fill(0.0);
    array<Type> A(nlvr,nlvr,n);
    A.fill(0.0);
    
    if(nlvr>0){
      
      if(nlvr>num_lv){
        for(int i=0; i<n; i++){
          A(0,0,i)=exp(lg_Ar(i));
        }
        if(lg_Ar.size()>n){
          for (int r=1; r<(nlvr); r++){
            for(int i=0; i<n; i++){
              A(r,0,i)=lg_Ar(r*n+i);
            }}
        }
      }
      
      
    if(num_lv>0){
      // log-Cholesky parametrization for A_i:s
      for (int d=0; d<(num_lv); d++){
        for(int i=0; i<n; i++){
          A(d+(nlvr-num_lv),d+(nlvr-num_lv),i)=exp(Au(d*n+i));
          // A(d,d,i)=exp(Au(d*n+i));
        }
      }
      if(Au.size()>(num_lv*n)){
        int k=0;
        for (int c=0; c<(num_lv); c++){
          for (int r=c+1; r<(num_lv); r++){
            for(int i=0; i<n; i++){
              A(r+(nlvr-num_lv),c+(nlvr-num_lv),i)=Au(num_lv*n+k*n+i);
              // A(r,c,i)=Au(nlvr*n+k*n+i);
              // A(c,r,i)=A(r,c,i);
            }
            k++;
          }}
      }
    }
      //set VA covariances for random rows to zero for quadratic model
      if((quadratic>0) && (nlvr>num_lv)){
        for(int i=0; i<n; i++){
          for (int d=0; d<(nlvr); d++){
            if(d!=0){
              A(d,0,i) = 0.0;
            }
          }
        }
      }
      for(int i=0; i<n; i++){
        if(nlvr == num_lv) nll.row(i).array() -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()+(u.row(i)*u.row(i).transpose()).sum()))/p;
         if(nlvr>num_lv) nll.row(i).array() -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(Cu.inverse()*(A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()-0.5*((u.row(i)*Cu.inverse())*u.row(i).transpose()).sum())/p;
        // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
      }
      nll.array() -= 0.5*(nlvr - log(Cu.determinant())*random(0))/p; //n*
      
      REPORT(nlvr);
      REPORT(A);
    }
    
    // Row/Site effects
    if(((random(0)>0) & (nlvr==num_lv)) & (rstruc>0)){
      if(rstruc == 1){
        if(cstruc==0){
          for (int j=0; j<p;j++){
            cQ.col(j) = cQ.col(j) + 0.5*(dr*Ar.matrix());
            eta.col(j) = eta.col(j) + dr*r0;
          }
          for (int i=0; i<nr; i++) {//i<n //!!!
            nll.array() -= 0.5*(1 + log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2) - log(sigma))/(n*p)*random(0);
          }
        } else {
          int j,d,r;
          
          matrix<Type> Sr(nr,nr);
          if(cstruc==1){// AR1 covariance
            Type rho = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*pow(rho,(d-j))*sigma;       // ar1 correlation
                Sr(j,d)=Sr(d,j);
              }
            }
          } else if(cstruc==2){// exp decaying
            Type alf = exp(log_sigma(1));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*exp(-(dc(d)-dc(j))/alf)*sigma;
                Sr(j,d)=Sr(d,j);
              }
            }
          } else {// Compound Symm  if(cstruc==3)
            Type rhob = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*rhob*sigma;
                Sr(j,d)=Sr(d,j);
              }
            }
          }
          
          matrix<Type> Arm(nr,nr);
          for (d=0; d<(nr); d++){
            Arm(d,d)=Ar(d);
          }
          
          if(lg_Ar.size()>nr){
            int k=0;
            for (d=0; d<(nr); d++){
              for (r=d+1; r<(nr); r++){
                Arm(r,d)=lg_Ar(nr+k);
                k++;
              }}
          }
          
          for (j=0; j<p;j++){
            cQ.col(j) = cQ.col(j) + 0.5*(dr*(Arm*Arm.transpose()).diagonal().matrix());
            eta.col(j) = eta.col(j) + dr*r0;
          }
          
          nll.array() -= 0.5*(log((Arm*Arm.transpose()).determinant()) - (Sr.inverse()*(Arm*Arm.transpose())).diagonal().sum()-(r0.transpose()*(Sr.inverse()*r0)).sum())/(n*p);// log(det(Ar_i))-sum(trace(Sr^(-1)Ar_i))*0.5 + ar_i*(Sr^(-1))*ar_i
          
          nll.array() -= 0.5*(nr-log(Sr.determinant()))/(n*p);
            // REPORT(Arm);
            // REPORT(Sr);
        }
        
      } else if(rstruc == 2){
        int i,j,d,r;
        matrix<Type> Sr(times,times);
        
        if(cstruc==1){// AR1 covariance
          Type rho = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*pow(rho,(d-j))*sigma;       // ar1 correlation
              Sr(j,d)=Sr(d,j);
            }
          }
        } else if(cstruc==2){// exp decaying
          Type alf = exp(log_sigma(1));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*exp(-(dc(d)-dc(j))/alf)*sigma;
              Sr(j,d)=Sr(d,j);
            }
          }
        } else {// Compound Symm  if(cstruc==3)
          Type rhob = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*rhob*sigma;
              Sr(j,d)=Sr(d,j);
            }
          }
        }
        
        array<Type> Arm(times,times,nr);
        for(i=0; i<nr; i++){
          for (d=0; d<(times); d++){
            Arm(d,d,i)=Ar(i*times+d);
          }
        }
        if(lg_Ar.size()>(nr*times)){
          int k=0;
          for (d=0; d<(times); d++){
            for (r=d+1; r<(times); r++){
              for(int i=0; i<nr; i++){//i<nr
                Arm(r,d,i)=lg_Ar(nr*times+k*nr+i);
                // Arm(d,r,i)=Arm(r,d,i);
              }
              k++;
            }}
        }
        
        for (j=0; j<p;j++){
          for (i=0; i<nr; i++) {
            for (d=0; d<(times); d++){
              cQ(i*times + d,j) += 0.5*(Arm.col(i).matrix().row(d)*Arm.col(i).matrix().row(d).transpose()).sum(); //Arm(d,d,i); 
            }
          }
          eta.col(j).array() += r0.array();
        }
        r0.resize(times, nr);
        for (i=0; i<nr; i++) {
          nll.array() -= 0.5*(log((Arm.col(i).matrix()*Arm.col(i).matrix().transpose()).determinant()) - (Sr.inverse()*(Arm.col(i).matrix()*Arm.col(i).matrix().transpose())).diagonal().sum()-((r0.col(i).matrix()).transpose()*(Sr.inverse()*(r0.col(i).matrix()))).sum())/(n*p);
          // log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
        }
        nll.array() -= 0.5*nr*(times - log(Sr.determinant()))/(n*p);
        // REPORT(Arm);
        // REPORT(Sr);
      }
      REPORT(nr);
      REPORT(r0);
      // eta += dr*r0;
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
      if(Abb.size()>(l*p)){
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
      nll.array() -= 0.5*(l - log(S.determinant())*random(1))/n;//n*
      // REPORT(S);  
      // REPORT(sds);  
      // REPORT(xb);  
      // REPORT(Ab);  
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
    if((quadratic < 1) && (nlvr > 0)){
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          cQ(i,j) += 0.5*((newlam.col(j)).transpose()*((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()*newlam.col(j))).sum();
        }
      }
      eta += lam;
    }
    matrix <Type> e_eta(n,p);
    e_eta.fill(0.0);
    if((quadratic>0) && (nlvr > 0)){
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
      if(((family==1)|(family==4))|(family==8)){
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
        // REPORT(B);
        // REPORT(V);
      }
      //Binomial, Gaussian, Ordinal
      if(((family==2)|(family==3))|(family==7)){
        for (int i=0; i<n; i++) {
          Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
          for (int j=0; j<p;j++){
            cQ(i,j) += 0.5*(newlam.col(j)*newlam.col(j).transpose()*Acov).trace() + (D.col(j).matrix()*Acov*D.col(j).matrix()*Acov).trace() +2*(u.row(i)*D.col(j).matrix()*Acov*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*Acov*newlam.col(j)).value();
            eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
          }
        }
      }
    }
    // REPORT(eta);
    // REPORT(cQ);
    
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
          nll(i,j) -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j))) - log(M_PI)/2;
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
      
    } else if(family==5){ // Tweedie EVA
      Type v = extra(0);
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          // Tweedie log-likelihood:
          nll(i,j) -= dtweedie(y(i,j), exp(eta(i,j)), iphi(j), v, true);
          if (y(i,j) == 0) {
            // Hessian-trace part:
            nll(i,j) += (1/iphi(j)) * (2-v)*exp(2*eta(i,j))*exp(-v*eta(i,j)) * cQ(i,j);
          } else if (y(i,j) > 0) {
            nll(i,j) -= (1/iphi(j)) * (y(i,j)*(1-v)*exp((1-v)*eta(i,j)) - (2-v)*exp((2-v)*eta(i,j))) * cQ(i,j);
          }
        }
      }
    } else if((family==7) && (zetastruc == 1)){//ordinal
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
    } else if((family==7) && (zetastruc==0)){
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
              if((y(i,j)==l) && (l != ymax)){
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
      
    } else if(family==9) { // Beta EVA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          // define mu, mu' and mu''
          Type mu = 0.0;
          Type mu_prime = 0.0;
          Type mu_prime2 = 0.0;
          if (extra(0) == 0) { // logit

            CppAD::vector<Type> z(4);
            z[0] = eta(i,j);
            z[1] = 0;
            z[2] = 1/(1+exp(-z[0]));
            z[3] = exp(z[0])/(exp(z[0])+1);

            mu = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
            mu_prime = mu * (1-mu);
            mu_prime2 = mu_prime * (1-2*mu);

          } else if (extra(0) == 1) { // probit
            mu = pnorm(eta(i,j), Type(0), Type(1));
            mu_prime = dnorm(eta(i,j), Type(0), Type(1));
            mu_prime2 = (-eta(i,j))*mu_prime;
          }
          CppAD::vector<Type> a(2);
          CppAD::vector<Type> b(2);
          a[0] = mu*iphi(j);
          a[1] = 1;
          b[0] = (1-mu)*iphi(j);
          b[1] = 1;
          CppAD::vector<Type> aa = a;
          CppAD::vector<Type> bb = b;
          aa[1] = 2;
          bb[1] = 2;
          Type dig_a = Type(atomic::D_lgamma(a)[0]);
          Type dig_b = Type(atomic::D_lgamma(b)[0]);
          Type trig_a = Type(atomic::D_lgamma(aa)[0]);
          Type trig_b = Type(atomic::D_lgamma(bb)[0]);
    //       
          nll(i,j) -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
          nll(i,j) -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
          nll(i,j) -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
    //       
        }
      }
    }
    // nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
    
  } else {
    eta += offset;
    if(random(0)==0){
      eta += r0*xr;
    }
    
    if(nlvr>0){
      eta += lam;
    }
    
    // Include random slopes if random(1)>0
    if(random(1)>0){
      vector<Type> sdsv = exp(sigmaB);
      density::UNSTRUCTURED_CORR_t<Type> neg_log_MVN(sigmaij);
      for (int j=0; j<p;j++){
        nll.col(j).array() += VECSCALE(neg_log_MVN,sdsv)(vector<Type>(Br.col(j)))/n;
      }
      eta += xb*Br;
    }
    
    // Row/Site effects
    if(((random(0)>0) & (nlvr==num_lv)) & (rstruc>0)){
      int i,j,d;
      
      if(rstruc == 1){
        if(cstruc ==0){
          matrix<Type> Sr(1,1);
          Sr(0,0) = sigma*sigma;
          MVNORM_t<Type> mvnorm(Sr);
          for (int i=0; i<nr; i++) {
            nll.array() += mvnorm(r0.row(i))/(n*p);
            // nll += dnorm(r0(i), Type(0), sigma);
          }
        } else {
          matrix<Type> Sr(nr,nr);
          if(cstruc==1){// AR1 covariance
            Type rho = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*pow(rho,(d-j))*sigma;       // ar1 correlation
                Sr(j,d)=Sr(d,j);
              }
            }
          } else if(cstruc==2){// exp decaying
            Type alf = exp(log_sigma(1));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*exp(-(dc(d)-dc(j))/alf)*sigma;
                Sr(j,d)=Sr(d,j);
              }
            }
          } else {// Compound Symm  if(cstruc==3)
            Type rhob = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
            for (d=0;d<nr;d++) {
              Sr(d,d)=sigma*sigma;
              for (j=0;j<d;j++){
                Sr(d,j)=sigma*rhob*sigma;
                Sr(j,d)=Sr(d,j);
              }
            }
          }
            MVNORM_t<Type> mvnorm(Sr);
            nll.array() += mvnorm(r0.col(0))/(n*p);
        }
        
        for (int j=0; j<p;j++){
          eta.col(j) = eta.col(j) + dr*r0;
        }
      } else {
        matrix<Type> Sr(times,times);
        
        if(cstruc==1){// AR1 covariance
          Type rho = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*pow(rho,(d-j))*sigma;       // ar1 correlation
              Sr(j,d)=Sr(d,j);
            }
          }
        } else if(cstruc==2){// exp decaying
          Type alf = exp(log_sigma(1));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*exp(-(dc(d)-dc(j))/alf)*sigma;
              Sr(j,d)=Sr(d,j);
            }
          }
        } else {// Compound Symm  if(cstruc==3)
          Type rhob = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          for (d=0;d<times;d++) {
            Sr(d,d)=sigma*sigma;
            for (j=0;j<d;j++){
              Sr(d,j)=sigma*rhob*sigma;
              Sr(j,d)=Sr(d,j);
            }
          }
        }
        
        
        for (j=0; j<p;j++){
          // cQ.col(j).array() += 0.5*Ar.array();
          eta.col(j).array() += r0.array();
        }
        MVNORM_t<Type> mvnorm(Sr);
        r0.resize(times, nr);
        for (i=0; i<nr; i++) {
          nll.array() += mvnorm(vector <Type> (r0.col(i)))/(n*p);
          // nll -= 0.5*(log((Arm.col(i).matrix()*Arm.col(i).matrix().transpose()).determinant()) - (Sr.inverse()*(Arm.col(i).matrix()*Arm.col(i).matrix().transpose())).diagonal().sum()-((r0.col(i).matrix()).transpose()*(Sr.inverse()*(r0.col(i).matrix()))).sum());// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
        }
        // REPORT(Sr);
        // REPORT(rho);
      }
      // REPORT(r0);
      // REPORT(eta);
      // REPORT(nlvr);
      // REPORT(nr);
      // REPORT(sigma);
      // REPORT(log_sigma);
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
    //latent variables and random site effects (r_i,u_i) from N(0,Cu)
    if(nlvr>0){
      
      MVNORM_t<Type> mvnorm(Cu);
      for (int i=0; i<n; i++) {
        nll.row(i).array() += mvnorm(u.row(i))/p;
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

    return nll.sum();
}
