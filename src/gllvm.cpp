#define TMB_LIB_INIT R_init_gllvm
#include <TMB.hpp>
#include<math.h>
#include "distrib.h"

//--------------------------------------------------------
//GLLVM
//Authors: Jenni Niku, Bert van der Veen, Pekka Korhonen
//------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  //declares all data and parameters used
  DATA_MATRIX(y); // matrix of responses
  DATA_MATRIX(x); // matrix of covariates
  DATA_MATRIX(x_lv); // matrix of covariates for Reduced Rank and/or constrained ord
  DATA_MATRIX(xr); 
  DATA_MATRIX(xb); // envs with random slopes
  DATA_ARRAY(dr0); // design matrix for rows, (times, n, nr)
  DATA_MATRIX(offset); //offset matrix
  
  PARAMETER_MATRIX(r0); // site/row effects
  PARAMETER_MATRIX(b); // matrix of species specific intercepts and coefs
  PARAMETER_MATRIX(bH); // matrix of species specific intercepts and coefs for beta hurdle model
  PARAMETER_MATRIX(B); // coefs of 4th corner model
  PARAMETER_MATRIX(Br); // random slopes for envs
  PARAMETER_MATRIX(b_lv); //slopes for RRR and constrained ord, VA means for random slopes
  //Left columns are for constrained ordination, Right for RRR
  PARAMETER_VECTOR(sigmaLV);//SD for LV
  PARAMETER_VECTOR(lambda); // lv loadings
  PARAMETER_MATRIX(lambda2);// quadratic lv loadings
  PARAMETER_MATRIX(thetaH);// hurdle model lv loadings
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi); // dispersion params/extra zero probs for ZIP
  PARAMETER_VECTOR(sigmaB); // sds for random slopes
  PARAMETER_VECTOR(sigmab_lv); // sds for random slopes constr. ord.
  PARAMETER_VECTOR(sigmaij);// cov terms for random slopes covariance
  PARAMETER_VECTOR(log_sigma);// log(SD for row effect) and 
  PARAMETER_MATRIX(rho_lvc);// correlation parameters for correlated LVs, matrix of q x 1 for corExp/corCS, qx2 for Matern
  
  DATA_INTEGER(num_lv); // number of lvs
  DATA_INTEGER(num_lv_c); //number of constrained lvs
  DATA_INTEGER(num_RR); //number of RRR dimensions
  DATA_INTEGER(num_corlv); //number of correlated lvs
  DATA_INTEGER(family); // family index
  DATA_INTEGER(quadratic); // quadratic model, 0=no, 1=yes
  
  PARAMETER_VECTOR(Au); // variational covariances for u
  PARAMETER_VECTOR(lg_Ar); // variational covariances for r0
  PARAMETER_VECTOR(Abb);  // variational covariances for Br
  // PARAMETER_VECTOR(scaledc);// scale parameters for dc, of length of dc.cols()
  PARAMETER_VECTOR(Ab_lv); //variational covariances for b_lv
  PARAMETER_VECTOR(zeta); // ordinal family param
  
  DATA_VECTOR(extra); // extra values, power of 
  DATA_INTEGER(method);// 0=VA, 1=LA, 2=EVA
  DATA_INTEGER(model);// which model, basic or 4th corner
  DATA_IVECTOR(random);//(0)1=random, (0)0=fixed row params, for Br: (1)1 = random slopes, (1)0 = fixed, for b_lv: (2)1 = random slopes, (2)0 = fixed slopes
  DATA_INTEGER(zetastruc); //zeta param structure for ordinal model
  DATA_INTEGER(rstruc); //Type for random rows. default = 0, when same as u:s. If 1, dr0 defines the structure. If 2, Points within groups has covariance struct defined by cstruc
  DATA_INTEGER(times); //number of time points
  DATA_INTEGER(cstruc); //correlation structure for row.params, 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm
  DATA_MATRIX(dc); //coordinates for sites, used for exponentially decaying cov. struc
  DATA_INTEGER(Astruc); //Structure of the variational covariance, 0=diagonal, 1=RR, (2=sparse cholesky not implemented yet)
  DATA_IMATRIX(NN); //nearest neighbours,
  
  
  matrix<Type> dr = dr0.matrix();
  // REPORT(dr);
  int Klv = x_lv.cols();
  int n = y.rows();
  int p = y.cols();
  // int nt =n;
  int nr =n;
  int nu =n; //CorLV
  if(rstruc>0){
    nr = dr.cols();
  }
  if(num_corlv>0){ //CorLV
    nu = dr.cols();
  }
  
  int l = xb.cols();
  vector<Type> iphi = exp(lg_phi);
  vector<Type> Ar = exp(lg_Ar);
  Type sigma = exp(log_sigma(0));
  
  // Set first row param to zero, if row effects are fixed
  if(random(0)<1){  r0(0,0) = 0;}
  int nlvr = num_lv+num_lv_c;//treating constr. ord random slopes as a LV, to use existing infrastructure for integration
  
  matrix<Type> ucopy = u;
  if(num_corlv>0){
    nlvr=0; num_lv=0; num_lv_c=0;
    quadratic=0;
    num_RR=0;
  }
  
  // Distance matrix calculated from the coordinates
  matrix<Type> DiSc(dc.cols(),dc.cols());
  matrix<Type> dc_scaled(dc.rows(),dc.cols());
  // matrix<Type> DistM(dc.rows(),dc.rows());
  // if(((num_corlv>0) | (((random(0)>0) & (nlvr==(num_lv+num_lv_c))) & (rstruc>0))) & ((cstruc==2) | (cstruc>3))){
  //   matrix<Type> DiSc(dc.cols(),dc.cols());
  //   DiSc.fill(0.0);
  // 
  //   for(int j=0; j<dc.cols(); j++){
  //     DiSc(j,j) += 1/exp(2*scaledc(j));
  //     // dc.col(j) *= 1/exp(scaledc(j));
  //   }
  //   // sigma_lvc(0,0) = 0;
  // 
  //   DistM.fill(0.0);
  //   for (int d=0;d<dc.rows();d++) {
  //     for (int j=0;j<d;j++){
  //       DistM(d,j)=sqrt( ((dc.row(d)-dc.row(j))*DiSc*(dc.row(d)-dc.row(j)).transpose()).sum() ); // + extra(2);
  //   //     DistM(j,d)=DistM(d,j);
  //     }
  //   }
  // }
  
  // if row params are in the same form as LVs, let's put them together
  if((random(0)>0) & (n == nr)){
    nlvr++;
    
    if((num_lv+num_lv_c)>0){
      u.conservativeResize(u.rows(), u.cols()+1);
      // for (int i=0; i<n; i++){
      for (int q=(num_lv+num_lv_c); q>0; q--){
        // u(i,q) = u(i,q-1);
        u.col(q) = u.col(q-1);
      }
      // }
      
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
  
  Type nll = 0; // initial value of log-likelihood
  
  matrix<Type> b_lv2(x_lv.cols(),nlvr);
  matrix <Type> Delta(nlvr,nlvr);
  b_lv2.fill(0.0);
  Delta.fill(0.0);
  
  array<Type> D;
  if( (quadratic>0) && (method!=1)){
    D = array<Type> (nlvr+num_RR,nlvr+num_RR,p);
    D.fill(0.0);
  }else if(quadratic>0){
    D = array<Type> (nlvr,nlvr,p);
    D.fill(0.0);
  }
  
  matrix<Type> newlam(nlvr,p);
  matrix<Type> RRgamma(num_RR,p);
  //K*K*d or d*d*K
  int sbl12 = 0;
  int sbl3 = 0;
  if((sigmab_lv.size()==Klv)|(sigmab_lv.size()==Type(1))){
    sbl12 = num_lv_c + num_RR;
    sbl3 = Klv;
  }else if(sigmab_lv.size()==(num_lv_c+num_RR)){
    sbl12 = Klv;
    sbl3 = num_lv_c + num_RR;
  }
  array<Type> Sigmab_lv(sbl12,sbl12,sbl3);
  
  if(random(2)>0){
    sigmab_lv = exp(sigmab_lv);
    sigmab_lv *= sigmab_lv;
    
    Sigmab_lv.fill(0.0);
    if(sigmab_lv.size()==(num_lv_c+num_RR)){//randomB=="LV", Sigma_q = sigma_q I_klv
      for (int q=0; q<(num_RR+num_lv_c); q++){
        matrix <Type> temp = Sigmab_lv.col(q).matrix();
        temp.diagonal().array() = sigmab_lv(q);
        Sigmab_lv.col(q) = temp.array();
      }
    }else if(sigmab_lv.size()==Klv){//randomB=="P", Sigma_klv = sigma_klv I_d
      for (int klv=0; klv<Klv; klv++){
        matrix <Type> temp = Sigmab_lv.col(klv).matrix();
        temp.diagonal().array() = sigmab_lv(klv);
        Sigmab_lv.col(klv) = temp.array();
      }
    }else if(sigmab_lv.size()==Type(1)){
      for (int klv=0; klv<Klv; klv++){
        matrix <Type> temp = Sigmab_lv.col(klv).matrix();
        temp.diagonal().array() = sigmab_lv(0);
        Sigmab_lv.col(klv) = temp.array();
      }
    }
    // else if(sigmab_lv.size()==((num_RR+num_lv_c)*Klv)){//probably a bad idea
    //   for (int klv=0; klv<Klv; klv++){
    //     matrix <Type> temp = Sigmab_lv.col(klv).matrix();
    //     for (int d=0; d<num_RR; d++){
    //       temp(d,d) = sigmab_lv(d+klv*(num_RR+num_lv_c));
    //       //        temp.diagoanl() = sigmab_lv(Eigen::seq(klv*(num_RR+num_lv_c),klv*(num_RR+num_lv_c)+(num_RR+num_lv_c)-1));//doesnt work..
    //       
    //     }
    //     Sigmab_lv.col(klv) = temp.array();
    //   }
    // }
  }
  if((nlvr>0)|(num_RR>0)){
    
    if(nlvr>0){
      newlam.row(0).fill(1.0);
      Cu.diagonal().fill(1.0);
      
      if((random(0)>0) && (n == nr)){
        for (int d=1; d<nlvr; d++){
          // Delta(d,d) = exp(sigmaLV(d-1)); //!!!
          Delta(d,d) = fabs(sigmaLV(d-1));
        }
        Delta(0,0) = 1;
        Cu(0,0) = sigma*sigma;
        if(log_sigma.size()>1){
          for (int d=1; d<nlvr; d++){
            Cu(d,0) = log_sigma(d);
            Cu(0,d) = Cu(d,0);
          }
        }
      }else if((num_lv+num_lv_c)>0){
        for (int d=0; d<nlvr; d++){
          // Delta(d,d) = exp(sigmaLV(d));
          Delta(d,d) = fabs(sigmaLV(d));
        }
      }
    }
    //To create lambda as matrix upper triangle
    // put LV loadings into a matrix
    if (num_lv>0){
      int tri = 0;
      if((num_lv_c+num_RR)>0){
        //because the lambdas for constrained and unconstrained LVs are separately identifiable and in the same vector
        tri += (num_lv_c+num_RR)*p-((num_lv_c+num_RR)*(num_lv_c+num_RR)-(num_lv_c+num_RR))/2-(num_lv_c+num_RR); //num_lv_c-1+p+(num_lv_c-1)*p-((num_lv_c-1)*(num_lv_c-1-1))/2-2*(num_lv_c-1); //number of elements for num_lv
      }
      for (int j=0; j<p; j++){
        for (int i=0; i<num_lv; i++){
          if(j<i){
            newlam(i+nlvr-num_lv-num_lv_c,j) = 0;
          }else if (j == i){
            newlam(i+nlvr-num_lv,j) = 1;
          }else if(j>i){
            newlam(i+nlvr-num_lv,j) = lambda(j+i*p-(i*(i-1))/2-2*i+tri-1);
          }
        }
      }
    }
    //species scores for constrained ordination and RRR
    if ((num_lv_c+num_RR)>0){
      for (int j=0; j<p; j++){
        for (int i=0; i<(num_lv_c+num_RR); i++){
          if(i<num_lv_c){
            if (j < i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = 0;
            } else if (j == i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = 1;
            }else if (j > i){
              newlam(i+nlvr-num_lv-num_lv_c,j) = lambda(j+i*p-(i*(i-1))/2-2*i-1);//lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }else{
            if (j < i){
              RRgamma(i-num_lv_c,j) = 0;
            } else if (j == i){
              RRgamma(i-num_lv_c,j) = 1;
            }else if (j > i){
              RRgamma(i-num_lv_c,j) = lambda(j+i*p-(i*(i-1))/2-2*i-1);//lambda(i+j+i*p-(i*(i-1))/2-2*i);
            }
          }
          
        }
      }
    }
  }
  
  // Loadings for correlated latent variables //CorLV
  matrix<Type> newlamCor(num_corlv,p);
  // matrix <Type> Delta_clv(num_corlv,num_corlv);
  if((num_corlv)>0){
    //To create lambda as matrix upper triangle
    // put LV loadings into a matrix
    for (int j=0; j<p; j++){
      for (int i=0; i<num_corlv; i++){
        if(j<i){
          newlamCor(i,j) = 0;
        }else if (j == i){
          newlamCor(i,j) = 1;
          // newlamCor(i,j) = exp(sigmaLV(i));
        }else if(j>i){
          newlamCor(i,j) = lambda(j+i*p-(i*(i-1))/2-2*i-1);
        }
      }
    }
    for (int d=0; d<num_corlv; d++){
      // Delta_clv(d,d) = fabs(sigmaLV(d));
      newlamCor.row(d)*=fabs(sigmaLV(d));
    }
  }
  
  matrix<Type> mu(n,p);
  mu.fill(0.0);
  
  using namespace density;
  using namespace gllvm;
  
  
  // Variational approximation
  if((method<1) | (method>1)){
    
    //quadratic coefficients for ordination
    //if random rows, add quadratic coefficients for num_RR to D otherwise
    //they go into D_RR below
    //The ordering here is num_lv_c-num_lv-num_RR so that the code works for
    //fixed-effects B and random effects B
    //The order we need to pick them from lambda2 is 
    //num_lv_c-num_RR-num_lv however, to ensure everything on the R-side works
    //like it should. So, things are messy but they work.
    if((quadratic>0) && ((num_lv+num_lv_c+num_RR*random(2))>0)){
      if(nlvr>(num_lv+num_lv_c)){
        if(num_lv_c>0){
          if(lambda2.cols()==1){
            for (int j=0; j<p; j++){
              for (int q=1; q<(num_lv_c+1); q++){
                D(q,q,j) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=1; q<(num_lv_c+1); q++){
                D(q,q,j) = fabs(lambda2(q-1,j)); //full quadratic model
              }
            } 
          }
        }
        if((num_RR*random(2))>0){
          if(lambda2.cols()==1){
            //make sure that num_RR comes at the end..has to be
            //like this due to the difference between fixed and random Bs
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1); q<(num_lv_c+1+num_RR); q++){
                D(q+num_lv,q+num_lv,j) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1); q<(num_lv_c+1+num_RR); q++){
                D(q+num_lv,q+num_lv,j) = fabs(lambda2(q-1,j)); //full quadratic model
              }
            } 
          }
        }
        if(num_lv>0){
          if(lambda2.cols()==1){
            //make sure that num_lv is taken from the middle even with num_RR
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1+num_RR); q<(num_lv_c+1+num_RR+num_lv); q++){
                D(q-num_RR,q-num_RR,j) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1+num_RR); q<(num_lv_c+1+num_RR+num_lv); q++){
                D(q-num_RR,q-num_RR,j) = fabs(lambda2(q-1,j)); //full quadratic model
              }
            } 
          }
        }
      }else{
        if(num_lv_c>0){
          if(lambda2.cols()==1){
            for (int j=0; j<p; j++){
              for (int q=0; q<num_lv_c; q++){
                D(q,q,j) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=0; q<num_lv_c; q++){
                D(q,q,j) = fabs(lambda2(q,j)); //full quadratic model
              }
            } 
          }
        }
        if((num_RR*random(2))>0){
          if(lambda2.cols()==1){
            //make sure that num_RR comes at the end..has to be
            //like this due to the difference between fixed and random Bs
            for (int j=0; j<p; j++){
              for (int q=num_lv_c; q<(num_lv_c+num_RR); q++){
                D(q+num_lv,q+num_lv,j) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=num_lv_c; q<(num_lv_c+num_RR); q++){
                D(q+num_lv,q+num_lv,j) = fabs(lambda2(q,j)); //full quadratic model
              }
            } 
          }
        }
        if(num_lv>0){
          if(lambda2.cols()==1){
            //make sure that num_lv is taken from the middle even with num_RR
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+num_RR); q<(num_lv_c+num_RR+num_lv); q++){
                D(q-num_RR,q-num_RR,j) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+num_RR); q<(num_lv_c+num_RR+num_lv); q++){
                D(q-num_RR,q-num_RR,j) = fabs(lambda2(q,j)); //full quadratic model
              }
            } 
          }
        }
      }
    }
    
    // add offset
    eta += offset;
    // add fixed row effects
    if((random(0)==0)){
      eta += r0*xr;
    }
    
    matrix<Type> cQ(n,p);
    cQ.fill(0.0);
    
    array<Type> A;
    if( (random(2)>0) && (num_RR!=0)){
      A = array<Type> (nlvr+num_RR,nlvr+num_RR,n);
      A.fill(0.0);
    }else{
      A = array<Type> (nlvr,nlvr,n);//plus one because just trying
      A.fill(0.0);
    }
    
    // Set up variational covariance matrix for LVs 
    if(nlvr>0){
      // Include variational covs of row effects, if structure is same for both
      if(nlvr>(num_lv+num_lv_c)){
        for(int i=0; i<n; i++){
          A(0,0,i)=exp(lg_Ar(i));
        }
        if(lg_Ar.size()>n){
          for (int r=1; r<nlvr; r++){
            for(int i=0; i<n; i++){
              A(r,0,i)=lg_Ar(r*n+i);
            }}
        }
      }
      
      
      if((num_lv+num_lv_c)>0){
        // log-Cholesky parametrization for A_i:s
        // don't include num_RR for random slopes, comes in later
        for (int d=0; d<(num_lv+num_lv_c); d++){
          for(int i=0; i<n; i++){
            A(d+(nlvr-num_lv-num_lv_c),d+(nlvr-num_lv-num_lv_c),i)=exp(Au(d*n+i));
            // A(d,d,i)=exp(Au(d*n+i));
          }
        }
        if(Au.size()>((num_lv+num_lv_c)*n)){
          int k=0;
          for (int c=0; c<(num_lv+num_lv_c); c++){
            for (int r=c+1; r<(num_lv+num_lv_c); r++){
              for(int i=0; i<n; i++){
                A(r+(nlvr-num_lv-num_lv_c),c+(nlvr-num_lv-num_lv_c),i)=Au((num_lv+num_lv_c)*n+k*n+i);
                // A(r,c,i)=Au(nlvr*n+k*n+i);
                // A(c,r,i)=A(r,c,i);
              }
              k++;
            }}
        }
      }
      
      //set VA covariances for random rows to zero for quadratic model
      //but not with quadratic model. constrained LVs, and row-eff.
      if((quadratic>0)&&(nlvr>(num_lv+num_lv_c))&&((num_lv+num_lv_c+num_RR*random(2))>0)){
        for(int i=0; i<n; i++){
          for (int d=0; d<nlvr; d++){
            if(d!=0){
              A(d,0,i) = 0.0;
            }
          }
        }
      }
      
      // // Add VA terms to logL
      if(random(2)<1){
        //Go this route if no random Bs
        for(int i=0; i<n; i++){
          if(nlvr == (num_lv+num_lv_c)) nll -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()+(u.row(i)*u.row(i).transpose()).sum()));
          if(nlvr>(num_lv+num_lv_c)) nll -= (((vector <Type> (A.col(i).matrix().diagonal())).log()).sum() - 0.5*(Cu.inverse()*(A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()).diagonal().sum()-0.5*((u.row(i)*Cu.inverse())*u.row(i).transpose()).sum());
          // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
          nll -= 0.5*(nlvr - log(Cu.determinant())*random(0));
        }
        //scale LVs with standard deviations, as well as the VA covariance matrices
        u *= Delta;
        
        for (int i=0; i<n; i++) {
          A.col(i) = (Delta*A.col(i).matrix()).array(); 
        }
      }else{
        //Go this route with random Bs (since size of Cu and A.col(i) are then not the same)
        matrix <Type>Atemp(nlvr,nlvr);
        Atemp.fill(0.0);
        for(int i=0; i<n; i++){
          Atemp = A.col(i).matrix().topLeftCorner(nlvr,nlvr);//to exlcude the 0 rows & columns for num_RR
          if(nlvr == (num_lv+num_lv_c)) nll -= (((vector <Type> (Atemp.diagonal())).log()).sum() - 0.5*((Atemp*Atemp.transpose()).diagonal().sum()+(u.row(i)*u.row(i).transpose()).sum()));
          if(nlvr>(num_lv+num_lv_c)) nll -= (((vector <Type> (Atemp.diagonal())).log()).sum() - 0.5*(Cu.inverse()*(Atemp*Atemp.transpose()).matrix()).diagonal().sum()-0.5*((u.row(i)*Cu.inverse())*u.row(i).transpose()).sum());
          // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
          nll -= 0.5*(nlvr - log(Cu.determinant())*random(0));
        }
        
        
        //scale LVs with standard deviations, as well as the VA covariance matrices
        u *= Delta;
        if(num_RR>0){
          Delta.conservativeResize(nlvr+num_RR,nlvr+num_RR);
          for(int d=nlvr; d<(nlvr+num_RR); d++){
            Delta.col(d).fill(0.0);
            Delta.row(d).fill(0.0);
          } 
        }
        
        for (int i=0; i<n; i++) {
          A.col(i) = (Delta*A.col(i).matrix()).array(); 
        }
      }
      
    }
    
    //random slopes for constr. ord.
    if((random(2)>0) && ((num_RR+num_lv_c)>0)){
      //resize A, u, D, and add RRGamma to newlam.
      //add columns to u on the right for num_RR with random slopes
      //resize u
      if(num_RR>0){
        u.conservativeResize(n, nlvr + num_RR); 
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          u.col(d).fill(0.0);
        }
        
        //resize and fill newlam, we don't use RRgamma further with random Bs
        //easiest to do is slap RRgamma at the end of newlam
        //this makes the order of newlam, A, u, and D inconsistent with the R-side of things
        //nicer would be to have to same order as in R, but that isn't possible since
        //it requires going down the same route for fixed and random B
        //which would only work with diagonal of 0s in A
        //And that needs to be invertible for the quadratic case, so that is not possible
        //potentially adjust in the future if I find out an alternative route.
        newlam.conservativeResize(nlvr+num_RR,p);
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          newlam.col(d).fill(0.0);
        }
        nlvr += num_RR;
        newlam.bottomRows(num_RR) = RRgamma;
      }
      
      // Variational covariance for random slopes
      // log-Cholesky parametrization for Ab_k:s
      array<Type> AB_lv(sbl12,sbl12,sbl3);
      
      AB_lv.fill(0.0);
      for (int q=0; q<(sbl12); q++){
        for(int d=0; d<sbl3; d++){
          AB_lv(q,q,d)=exp(Ab_lv(q*sbl3+d));
        }
      }
      if(Ab_lv.size()>((sbl12)*sbl3)){
        int k=0;
        for (int c=0; c<(sbl12); c++){
          for (int r=c+1; r<(sbl12); r++){
            for(int d=0; d<sbl3; d++){
              AB_lv(r,c,d)=Ab_lv((sbl12)*sbl3+k*sbl3+d);
              // Ab(c,r,j)=Ab(r,c,j);
            }
            k++;
          }}
      }
      //VA likelihood parts for random slopes
      if(sbl3 == Klv){//randomB="P" & randomB=="single"
        for(int klv=0; klv<Klv; klv++){
          nll -= ((((vector <Type> (AB_lv.col(klv).matrix().diagonal())).log()).sum() - 0.5*((Sigmab_lv.col(klv).matrix()).inverse()*(AB_lv.col(klv).matrix()*AB_lv.col(klv).matrix().transpose())).trace()-0.5*(b_lv.row(klv)*((Sigmab_lv.col(klv).matrix()).inverse()*b_lv.row(klv).transpose())).sum()));// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
          nll -= 0.5*(num_lv_c+num_RR-((vector <Type> (Sigmab_lv.col(klv).matrix().diagonal())).log()).sum());
        }
      }else{//randomB="LV"
        for(int q=0; q<(num_lv_c+num_RR); q++){
          nll -= ((((vector <Type> (AB_lv.col(q).matrix().diagonal())).log()).sum() - 0.5*((Sigmab_lv.col(q).matrix()).inverse()*(AB_lv.col(q).matrix()*AB_lv.col(q).matrix().transpose())).trace()-0.5*(b_lv.col(q).transpose()*((Sigmab_lv.col(q).matrix()).inverse()*b_lv.col(q))).sum()));// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
          nll -= 0.5*(Klv-((vector <Type> (Sigmab_lv.col(q).matrix().diagonal())).log()).sum());
        }
      }
      
      
      //now rebuild A and u with covariances for random slopes so that existing infrastructure below can be used
      //in essence, q(XBsigmab_lv + eDelta) ~ N(uDelta + \sum \limits^K X_ik b_lv_k , Delta A Delta + \sum \limits^K X_ik^2 AB_lv_k )
      //so build u and A accordingly (and note covariance due to Bs if num_lv_c and num_RR > 0)
      
      if(sbl3 == Klv){//variance per predictor
        if((num_lv_c>0) && (num_RR == 0)){
          // matrix <Type> b_lv2 =  b_lv;//.leftCols(num_lv_c);
          if((random(0)>0) && (n == nr)){
            u.middleCols(1, num_lv_c) += x_lv*b_lv;
          }else{
            u.leftCols(num_lv_c) += x_lv*b_lv;
          }
          
          matrix<Type> temp(nlvr,nlvr);
          temp.fill(0.0);
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0);
          if((random(0) == 0) || (n != nr)){
            for(int i=0; i<n; i++){
              temp = A.col(i).matrix()*A.col(i).matrix().transpose();
              for(int klv=0; klv<Klv; klv++){
                temp.topLeftCorner(num_lv_c,num_lv_c) += x_lv(i,klv)*x_lv(i,klv)*AB_lv.col(klv).matrix()*AB_lv.col(klv).matrix().transpose();//cholesky of variance block for num_lv_c
              }
              L =  temp.llt().matrixL();//can't do only a part due to potential covariance with num_lv
              A.col(i) =  L.array();//have to recompute cholesky of covariance due to summation
            }
          }else if((n == nr) && (random(0)>0)){//if row effects are included in u and A
            for(int i=0; i<n; i++){
              temp = A.col(i).matrix()*A.col(i).matrix().transpose();
              for(int klv=0; klv<Klv; klv++){
                temp.block(1,1,num_lv_c,num_lv_c) += x_lv(i,klv)*x_lv(i,klv)*AB_lv.col(klv).matrix()*AB_lv.col(klv).matrix().transpose();//cholesky of variance block for num_lv_c
              }
              L =  temp.llt().matrixL();//can't do only a part due to potential covariance with num_lv
              A.col(i) =  L.array();//have to recompute cholesky of covariance due to summation
            }
          }
        }
        
        if((num_RR>0) && (num_lv_c == 0)){
          // matrix <Type> b_lv3 =  b_lv;//.rightCols(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv;
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0);
          for(int i=0; i<n; i++){
            L = A.col(i).matrix();
            for(int klv=0; klv<Klv; klv++){
              L.bottomRightCorner(num_RR,num_RR) += x_lv(i,klv)*x_lv(i,klv)*AB_lv.col(klv).matrix()*AB_lv.col(klv).matrix().transpose();//cholesky of variance block for num_lv_c
            }
            
            L.bottomRightCorner(num_RR,num_RR) =  (L.bottomRightCorner(num_RR,num_RR)).llt().matrixL();//block diagonal structure so only need to re-do part of this matrix, the bottom right
            A.col(i) =  L.array();//have to recompute cholesky of covariance due to summation
          }
        }
        
        //separate case because now we have both, so that we need to build the covariances between these two as well, and put them back in the right place in A
        if((num_RR>0) && (num_lv_c>0)){
          matrix <Type> b_lv2 =  b_lv.leftCols(num_lv_c);
          matrix <Type> b_lv3 =  b_lv.rightCols(num_RR);
          
          matrix<Type> temp(num_RR+num_lv_c,num_RR+num_lv_c);
          temp.fill(0.0);
          
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0);
          
          if((random(0)>0) && (n == nr)){
            u.middleCols(1, num_lv_c) += x_lv*b_lv2;
          }else{
            u.leftCols(num_lv_c) += x_lv*b_lv2;
          }
          u.rightCols(num_RR) += x_lv*b_lv3;
          for(int i=0; i<n; i++){
            L = A.col(i).matrix()*A.col(i).matrix().transpose();
            temp.setZero();
            for(int klv=0; klv<Klv; klv++){
              temp +=  x_lv(i,klv)*x_lv(i,klv)*AB_lv.col(klv).matrix()*AB_lv.col(klv).matrix().transpose();//num_lv_c variance block
            }
            
            if((random(0)==0) || (n != nr)){
              L.topLeftCorner(num_lv_c,num_lv_c) += temp.topLeftCorner(num_lv_c,num_lv_c);
              L.bottomRightCorner(num_RR,num_RR) += temp.bottomRightCorner(num_RR,num_RR);
              
              L.bottomLeftCorner(num_RR,num_lv_c) += temp.bottomLeftCorner(num_RR,num_lv_c);
              L.topRightCorner(num_lv_c,num_RR) += temp.topRightCorner(num_lv_c,num_RR);
              
            }else if ((random(0) > 0) && (n == nr)){//if row effects are included in u and A
              L.block(1,1,num_lv_c,num_lv_c) += temp.topLeftCorner(num_lv_c,num_lv_c);
              L.bottomRightCorner(num_RR,num_RR) += temp.bottomRightCorner(num_RR,num_RR);
              
              L.block(nlvr-num_RR,1,num_RR,num_lv_c) += temp.bottomLeftCorner(num_RR,num_lv_c);//should be bottom left corner
              L.block(1,nlvr-num_RR,num_lv_c,num_RR) += temp.topRightCorner(num_lv_c,num_RR);//should be top right corner
              
            }
            L = L.llt().matrixL();
            A.col(i) = L.array();
          }
          
        }
      }else if(sbl3 == (num_lv_c+num_RR)){//variance per LV
        if(num_lv_c>0){
          matrix <Type> b_lv2 =  b_lv.leftCols(num_lv_c);
          u.leftCols(num_lv_c) += x_lv*b_lv2;
        }
        if(num_RR>0){
          matrix <Type> b_lv3 =  b_lv.rightCols(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv3;
        }
        //much easier, since we assume independence between LVs
        array <Type> temp(nlvr,nlvr,n);
        temp.fill(0.0);
        matrix <Type> L(nlvr,nlvr);
        L.fill(0.0);
        
        if((random(0)<1) || (n != nr)){
          if(num_lv_c>0){
            for(int i=0; i<n; i++){
              for(int q=0; q<num_lv_c; q++){
                temp(q,q,i) = (x_lv.row(i)*(AB_lv.col(q).matrix()*AB_lv.col(q).matrix().transpose())*x_lv.row(i).transpose()).sum();
              }
            }
          }
          if(num_RR>0){
            for(int i=0; i<n; i++){
              for(int q=(num_lv_c+num_lv); q<(num_lv_c+num_lv+num_RR); q++){
                temp(q,q,i) = (x_lv.row(i)*(AB_lv.col(q-num_lv).matrix()*AB_lv.col(q-num_lv).matrix().transpose())*x_lv.row(i).transpose()).sum();
              }
            }
          }
          for(int i=0; i<n; i++){
            L = A.col(i).matrix()*A.col(i).matrix().transpose() + temp.col(i).matrix();
            L = L.llt().matrixL();
            A.col(i) = L.array();
          }
          REPORT(temp);
        }else if((random(0)>1) && (n == nr)){//if row effects are included in u and A
          array <Type> temp(nlvr,nlvr,n);
          temp.fill(0.0);
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0);
          
          if(num_lv_c>0){
            matrix <Type> b_lv2 =  b_lv.leftCols(num_lv_c);
            u.middleCols(1, num_lv_c) += x_lv*b_lv2;
          }
          if(num_RR>0){
            matrix <Type> b_lv3 =  b_lv.rightCols(num_RR);
            u.rightCols(num_RR) += x_lv*b_lv3;
          }
          
          if(num_lv_c>0){
            for(int i=0; i<n; i++){
              for(int q=1; q<(num_lv_c+1); q++){
                temp(q,q,i) = (x_lv.row(i)*AB_lv.col(q-1).matrix()*AB_lv.col(q-1).matrix().transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          if(num_RR>0){
            for(int i=0; i<n; i++){
              for(int q=(num_lv_c+num_lv+1); q<(num_lv_c+num_lv+num_RR+1); q++){
                temp(q,q,i) = (x_lv.row(i)*AB_lv.col(q-num_lv-1).matrix()*AB_lv.col(q-num_lv-1).matrix().transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          for(int i=0; i<n; i++){
            L = A.col(i).matrix()*A.col(i).matrix().transpose() + temp.col(i).matrix();
            L = L.llt().matrixL();
            A.col(i) = L.array();
          }
        }
      }
      // REPORT(u);
      REPORT(A);
      REPORT(AB_lv);
    }
    
    
    // Include random slopes if random(1)>0
    if(random(1)>0){
      matrix<Type> sds(l,l);
      sds.fill(0.0);
      sds.diagonal() = exp(sigmaB);
      matrix<Type> S=sds*UNSTRUCTURED_CORR(sigmaij).cov()*sds;
      // REPORT(S);
      
      // Variational covariance for random slopes
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
       A is a num.lv x num.lv x n array, theta is p x num.lv matrix*/
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          cQ(i,j) += 0.5*((xb.row(i))*((Ab.col(j).matrix()*Ab.col(j).matrix().transpose()).matrix()*xb.row(i).transpose())).sum();
        }
        nll -= ((((vector <Type> (Ab.col(j).matrix().diagonal())).log()).sum() - 0.5*(S.inverse()*(Ab.col(j).matrix()*Ab.col(j).matrix().transpose()).matrix()).trace()-0.5*(Br.col(j).transpose()*(S.inverse()*Br.col(j))).sum()));// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
      }
      eta += xb*Br;
      nll -= 0.5*(l - log(S.determinant())*random(1))*p;//n*
    }
    
    
    if(model<1){
      // basic gllvm, gllvm.TMB.R
      eta += x*b;
    } else {
      // Fourth corner model TMB.trait.R
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          eta(i,j)+=b(0,j)*extra(1)+eta1(m,0); //extra(1)=0 if beta0comm=TRUE
          m++;
        }
      }
    }
    
    matrix <Type> e_eta(n,p);
    e_eta.fill(0.0);
    //components for reduced rank regression terms
    if((num_RR>0) && (random(2)<1)){
      //predictor coefficients RRR.  num_RR comes after num_lv_c
      //Since later dimensions are more likely to have less residual variance
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      
      //quadratic terms for fixed-effects only RRR
      //-num_lv to ensure that we pick num_RR from the middle
      if(quadratic>0){
        matrix <Type> D_RR(num_RR,num_RR);
        D_RR.fill(0.0);
        //quadratic coefficients for RRR
        if(lambda2.cols()==1){
          for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
            D_RR(d-num_lv_c,d-num_lv_c) = fabs(lambda2(d,0));
          }
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
          
        }else{
          for (int j=0; j<p;j++){
            for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
              D_RR(d-num_lv_c,d-num_lv_c) = fabs(lambda2(d,j));
            }
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
            
          }
        }
        
      }
    }else if((quadratic>0) && (random(2)>0)){
      //slap D's at end for num_RR and random slopes
      //-num_lv to ensure that we pick num_RR from the middle
      if(nlvr>(num_lv+num_lv_c+(num_RR*random(2)))){
        if(lambda2.cols()==1){
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c+1); q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-1-num_lv,0)); //common tolerances model
            }
          }
        }else{
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c+1); q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-1-num_lv,j)); //full quadratic model
            }
          }
        }
        
      }else{
        if(lambda2.cols()==1){
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c); q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-num_lv,0)); //common tolerances model
            }
          }
        }else{
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c); q<nlvr; q++){
              D(q,q,j) = fabs(lambda2(q-num_lv,j)); //full quadratic model
            }
          }
        }
      }
    }
    
    
    // Structured Row/Site effects
    if(((random(0)>0) & (nlvr==(num_lv+num_lv_c))) & (rstruc>0)){
      // Group specific random row effects:
      if(rstruc == 1){
        if(cstruc==0){
          for (int j=0; j<p;j++){
            cQ.col(j) = cQ.col(j) + 0.5*(dr*Ar.matrix());
            eta.col(j) = eta.col(j) + dr*r0;
          }
          for (int i=0; i<nr; i++) {//i<n //!!!
            nll -= 0.5*(1 + log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2) - 2*log(sigma))*random(0); ///(n*p)
          }
        } else {
          // group specific random row effects, which are correlated between groups
          int j,d,r;
          
          matrix<Type> Sr(nr,nr);
          if(cstruc==1){// AR1 covariance
            Sr = gllvm::corAR1(sigma, log_sigma(1), nr);
          } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
            Sr = gllvm::corCS(sigma, log_sigma(1), nr);
          } else {
            DiSc.fill(0.0);
            for(int j=0; j<dc.cols(); j++){
              DiSc(j,j) += 1/exp(log_sigma(1+j));
            }
            dc_scaled = dc*DiSc;
            if(cstruc==2){// exp decaying
              Sr = gllvm::corExp(sigma, Type(0), nr, dc_scaled);
              // Sr = gllvm::corExp(sigma, (log_sigma(1)), nr, DistM);
            } else if(cstruc==4) {// Matern
              Sr = gllvm::corMatern(sigma, Type(0), log_sigma(dc.cols()+1), nr, dc_scaled);
              //   Sr = gllvm::corMatern(sigma, log_sigma(1), log_sigma(2), nr, DistM);
            }
          }
          
          // Variational covariance for row effects
          matrix<Type> Arm(nr,nr);
          for (d=0; d<(nr); d++){
            Arm(d,d)=Ar(d);
          }
          
          if((lg_Ar.size()>nr) & (Astruc>0)){ // unstructured Var.cov
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
          
          nll -= 0.5*(log((Arm*Arm.transpose()).determinant()) - (atomic::matinv(Sr)*(Arm*Arm.transpose())).diagonal().sum()-(r0.transpose()*(atomic::matinv(Sr)*r0)).sum());// /(n*p)log(det(Ar_i))-sum(trace(Sr^(-1)Ar_i))*0.5 + ar_i*(Sr^(-1))*ar_i
          
          nll -= 0.5*(nr-log(Sr.determinant()));
          // REPORT(Arm);
          // REPORT(Sr);
        }
        
      } else if(rstruc == 2){
        // site specific random row effects, which are correlated within groups
        int i,j,d,r;
        matrix<Type> Sr(times,times);
        
        // Define covariance matrix
        if(cstruc==1){// AR1 covariance
          Sr = gllvm::corAR1(sigma, log_sigma(1), times);
        } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
          Sr = gllvm::corCS(sigma, log_sigma(1), times);
        } else{ 
          DiSc.fill(0.0);
          for(int j=0; j<dc.cols(); j++){
            DiSc(j,j) += 1/exp(log_sigma(1+j));
          }
          dc_scaled = dc*DiSc;
          if(cstruc==2){// exp decaying
            Sr = gllvm::corExp(sigma, Type(0), times, dc_scaled);
            // Sr = gllvm::corExp(sigma, (log_sigma(1)), times, DistM);
          } else if(cstruc==4) {// Matern
            Sr = gllvm::corMatern(sigma, Type(0), log_sigma(dc.cols()+1), times, dc_scaled);
            // Sr = gllvm::corMatern(sigma, log_sigma(1), log_sigma(2), times, DistM);
          }
        }
        
        // Variational covariance for row effects
        array<Type> Arm(times,times,nr);
        for(i=0; i<nr; i++){
          for (d=0; d<(times); d++){
            Arm(d,d,i)=Ar(i*times+d);
          }
        }
        if((lg_Ar.size()>(nr*times)) & (Astruc>0)){ // unstructured Var.cov
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
          nll -= log(Arm.col(i).matrix().determinant()) + 0.5*( - (atomic::matinv(Sr)*(Arm.col(i).matrix()*Arm.col(i).matrix().transpose())).diagonal().sum()-((r0.col(i).matrix()).transpose()*(atomic::matinv(Sr)*(r0.col(i).matrix()))).sum());
          // log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
        }
        nll -= 0.5*nr*(times - log(Sr.determinant()));
        // Report
        // REPORT(Arm);
        // REPORT(Sr);
      }
      // REPORT(nr);
      // REPORT(r0);
      // eta += dr*r0;
    }
    
    // Correlated LVs
    if(num_corlv>0) { //CorLV
      int i,j,d;
      int arank = 2;
      matrix<Type> AQ(num_corlv,num_corlv);
      AQ.fill(0.0); AQ.diagonal().fill(1.0);
      
      if(ucopy.rows() == nu){
        REPORT(ucopy);
        REPORT(nu);
        
        if(cstruc==0){
          array<Type> Alvm(num_corlv,num_corlv,nu);
          Alvm.fill(0.0);
          eta += (dr*ucopy)*newlamCor;
          
          // Variational covariance for row effects
          for (int q=0; q<(num_corlv); q++){
            for (d=0; d<(nu); d++){
              Alvm(q,q,d)=exp(Au(q*nu+d));
            }
          }
          if((Astruc>0) & (Au.size()>((num_corlv)*nu))){//unstructured cov
            int k=0;
            for (int c=0; c<(num_corlv); c++){
              for (int r=c+1; r<(num_corlv); r++){
                for(d=0; d<nu; d++){
                  Alvm(r,c,d)=Au(nu*num_corlv+k*nu+d);
                }
                k++;
              }}
          }
          
          
          for (d=0; d<nu; d++) {
            nll -= log(Alvm.col(d).matrix().determinant()) + 0.5*( - (Alvm.col(d).matrix()*Alvm.col(d).matrix().transpose()).diagonal().sum() - (ucopy.row(d).matrix()*ucopy.row(d).matrix().transpose()).sum());
            // for (d=0; d<nu; d++)
            for (j=0; j<p;j++){
              cQ.col(j) += 0.5*dr.col(d)*((newlamCor.col(j).transpose()*(Alvm.col(d).matrix()*Alvm.col(d).matrix().transpose()))*newlamCor.col(j));
            }
          }
          nll -= 0.5*(nu*num_corlv);
          
        } else {
          vector<matrix<Type> > Slv(num_corlv);
          matrix<Type> Slvinv(nu,nu);
          Slvinv.fill(0.0);
          // matrix<Type> Slv(nu,nu);
          matrix<Type> uq(nu,1);
          eta += (dr*ucopy)*newlamCor;
          
          if(Astruc<3){
            array<Type> Alvm(nu,nu,num_corlv);
            Alvm.fill(0.0);
            // matrix<Type> Alvm(nu,nu);
            for(int q=0; q<num_corlv; q++){
              // site specific LVs, which are correlated between groups
              // Slv.fill(0.0);
              uq = ucopy.col(q);
              
              // group specific lvs
              if(cstruc==1){// AR1 covariance
                Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
              } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
                Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
              } else {
                DiSc.fill(0.0);
                for(int j=0; j<dc.cols(); j++){
                  DiSc(j,j) += 1/exp(rho_lvc(q,j));
                }
                dc_scaled = dc*DiSc;
                if(cstruc==2){// exp decaying
                  Slv(q) = gllvm::corExp(Type(1), Type(0), nu, dc_scaled);
                  // Slv(q) = gllvm::corExp(Type(1), (rho_lvc(q,0)), nu, DistM);
                } else if(cstruc==4) {// Compound Symm  if(cstruc==3)
                  Slv(q) = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), nu, dc_scaled);
                  // Slv(q) = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), nu, DistM);
                }
              }
              
              // Variational covariance for row effects
              for (d=0; d<(nu); d++){
                Alvm(d,d,q)=exp(Au(q*nu+d));
              }
              
              if((Astruc>0) & (Au.size() > nu*num_corlv)){//reduced rank cov
                // if(Au.size()>(times*nu*num_corlv)){}
                int k=0;
                if(Astruc==1){
                  for (d=0; d<nu; d++){
                    for (int r=d+1; r<(nu); r++){
                      Alvm(r,d,q)=Au(nu*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                } else if(Astruc==2) {
                  arank = NN.rows();
                  // arank = NN.cols();
                  // for (d=0; (d<nu); d++){
                  for (int r=0; r<(arank); r++){
                    Alvm(NN(r,0)-1,NN(r,1)-1,q)=Au(nu*num_corlv+k*num_corlv+q);
                    // int d2 = NN(d,r)-1;
                    // if(d2<d){
                    //   Alvm(d,d2)=Au(nu*num_corlv+k*num_corlv+q);
                    // } else {
                    // Alvm(d2,d)=Au(nu*num_corlv+k*num_corlv+q);
                    // }
                    k++;
                  }
                  // }
                }
                REPORT(k);
              }
              
              for (j=0; j<p;j++){
                cQ.col(j) = cQ.col(j) + 0.5*pow(newlamCor(q,j),2)*(dr*(Alvm.col(q).matrix()*Alvm.col(q).matrix().transpose()).diagonal().matrix());
              }
              Slvinv = atomic::matinv(Slv(q));
              nll -= log(Alvm.col(q).matrix().determinant()) + 0.5*(- (Slvinv*(Alvm.col(q).matrix()*Alvm.col(q).matrix().transpose())).diagonal().sum()-( uq.transpose()*(Slvinv*uq) ).sum());
              // nll -= log(Alvm.col(q).matrix().determinant()) + 0.5*(- (atomic::matinv(Slv(q))*(Alvm.col(q).matrix()*Alvm.col(q).matrix().transpose())).diagonal().sum()-( uq.transpose()*(atomic::matinv(Slv(q))*uq) ).sum());
              
              nll -= 0.5*(nu-log(Slv(q).determinant()));
              // REPORT(Arm);
              // REPORT(Sr);
              
            }
            REPORT(Alvm);
          } else if(num_corlv>1){
            matrix<Type> Alvm(nu,nu);
            Alvm.fill(0.0);
            //   // Kronecker Variational covariance
            for (d=0; d<(nu); d++){
              Alvm(d,d)=exp(Au(d));
            }
            //   //reduced rank cov
            //     // if(Au.size()>(times*nu*num_corlv)){}
            int k=0;
            arank = NN.rows();
            if(Au.size()>(nu+num_corlv*(num_corlv+1)/2)) {
              if(Astruc == 4) {
                for (int r=0; r<(arank); r++){
                  Alvm(NN(r,0)-1,NN(r,1)-1)=Au(nu+k);
                  k++;
                }
              } else if(Astruc == 3) {
                for (d=0; d<nu; d++){
                  for (int r=d+1; r<(nu); r++){
                    Alvm(r,d)=Au(nu+k);
                    k++;
                  }
                }
              }
            }
            
            for (d=0; d<num_corlv; d++){
              AQ(d,d)=exp(Au(nu+k));
              k++;
              for (int r=d+1; r<(num_corlv); r++){
                AQ(r,d)=Au(nu+k);
                k++;
              }
            }
            Alvm *= Alvm.transpose();
            AQ *= AQ.transpose();
            
            for (j=0; j<p;j++){
              cQ.col(j) = cQ.col(j) + 0.5*(dr*Alvm.diagonal())*((newlamCor.col(j).transpose()*AQ)*newlamCor.col(j));
            }
            nll -= 0.5*num_corlv*log(Alvm.determinant()) + 0.5*nu*log(AQ.determinant()) + 0.5*num_corlv*nu;
            //
            for(int q=0; q<num_corlv; q++){
              // site specific LVs, which are correlated between groups
              // Slv.fill(0.0);
              uq = ucopy.col(q);
              
              // group specific lvs
              if(cstruc==1){// AR1 covariance
                Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
              } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
                Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
              } else {
                DiSc.fill(0.0);
                for(int j=0; j<dc.cols(); j++){
                  DiSc(j,j) += 1/exp(rho_lvc(q,j));
                }
                dc_scaled = dc*DiSc;
                if(cstruc==2){// exp decaying
                  Slv(q) = gllvm::corExp(Type(1), Type(0), nu, dc_scaled);
                  // Slv(q) = gllvm::corExp(Type(1), (rho_lvc(q,0)), nu, DistM);
                } else if(cstruc==4) {// Compound Symm  if(cstruc==3)
                  Slv(q) = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), nu, dc_scaled);
                  // Slv(q) = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), nu, DistM);
                }
              }
              
              Slvinv = atomic::matinv(Slv(q));
              nll -= 0.5*(- AQ(q,q)*(Slvinv*Alvm).diagonal().sum()-( uq.transpose()*(Slvinv*uq) ).sum());
              // nll -= 0.5*(- AQ(q,q)*(atomic::matinv(Slv(q))*Alvm).diagonal().sum()-( uq.transpose()*(atomic::matinv(Slv(q))*uq) ).sum());
              nll -= -0.5*log(Slv(q).determinant());
            }
            //
            REPORT(Alvm);
          }
          REPORT(Slv);
        }
      } else {
        
        eta += ucopy*newlamCor;
        vector<matrix<Type> > Slv(num_corlv);
        // matrix<Type> Slv(times,times);
        matrix<Type> Slvinv(times,times);
        Slvinv.fill(0.0);
        
        // int acol = times;
        // if(Astruc>0){
        //   acol = arank;
        // }
        //     
        if(Astruc<3){
          
          array<Type> Alvm(times*nu,times*nu,num_corlv);
          Alvm.fill(0.0);
          // array<Type> Alvm(times,times,nu);
          // matrix<Type> uq(times*nu,1);
          
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            
            // Variational covariance for row effects
            //diagonal
            for(i=0; i<nu; i++){
              for (d=0; d<(times); d++){
                Alvm(i*times+d,i*times+d,q)=exp(Au(q*n+i*times+d));
              }
            }
            
            if((Astruc>0) & (Au.size() > nu*times*num_corlv)){//reduced rank cov
              // if(Au.size()>(times*nu*num_corlv)){}
              int k=0;
              if(Astruc==1){
                for(i=0; i<nu; i++){
                  for (d=0; ((d<arank) & (d<times)); d++){
                    // for (d=0; (d<times); d++){//(num_lv+num_lv_c)*n+k*n+q
                    // for (int r=d; r<(times); r++){
                    // if(r==d){
                    // Alvm(r,d,i)=exp(Au(k*num_corlv+q));
                    // } else {
                    // Alvm(r,d,i)=Au(k*num_corlv+q);
                    // }
                    for (int r=d+1; r<(times); r++){
                      Alvm(i*times+r,i*times+d,q)=Au(nu*times*num_corlv+k*num_corlv+q);
                      // Alvm(r,d,i)=Au(nu*times*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                }
              } else if(Astruc==2) {
                arank = NN.rows();
                for(i=0; i<nu; i++){
                  for (int r=0; r<(arank); r++){
                    Alvm(i*times+NN(r,0)-1,i*times+NN(r,1)-1,q)=Au(nu*times*num_corlv+k*num_corlv+q);
                    k++;
                  }
                }
                // arank = NN.cols();
                // for(i=0; i<nu; i++){
                //   for (d=0; (d<times); d++){
                //     for (int r=0; r<(arank); r++){
                //       int d2 = NN(d,r)-1;
                //       if(d2<d){
                //         Alvm(i*times+d,i*times+d2,q)=Au(nu*times*num_corlv+k*num_corlv+q);
                //         // k++;
                //       } else {
                //         Alvm(i*times+d2,i*times+d,q)=Au(nu*times*num_corlv+k*num_corlv+q);
                //       }
                //       k++;
                //     }
                //   }
                // }
              }
              // REPORT(k);
            }
            
            
            for (j=0; j<p;j++){
              for (i=0; i<(times*nu); i++) {
                cQ(i,j) += 0.5*pow(newlamCor(q,j),2)*(Alvm.col(q).matrix().row(i)*Alvm.col(q).matrix().row(i).transpose()).sum();
              }
            }
            nll -= log(Alvm.col(q).matrix().determinant());
            
            // Slv.fill(0.0);
            // Alvm.fill(0.0);
            // uq = ucopy.col(q);
            
            // Define covariance matrix
            if(cstruc==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else { 
              DiSc.fill(0.0);
              for(int j=0; j<dc.cols(); j++){
                DiSc(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled = dc*DiSc;
              if(cstruc==2){// exp decaying
                Slv(q) = gllvm::corExp(Type(1), Type(0), times, dc_scaled);
                // Slv(q) = gllvm::corExp(Type(1), (rho_lvc(q,0)), times, DistM);
              } else if(cstruc==4) {// matern
                Slv(q) = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), times, dc_scaled);
                // Slv(q) = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), times, DistM);
              }
            }
            
            nll -= 0.5*nu*(times - log(Slv(q).determinant()));
            
            Slvinv = atomic::matinv(Slv(q));
            for (i=0; i<nu; i++) {
              nll -=  0.5*(- (Slvinv*(Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix()*Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix().transpose())).diagonal().sum()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(Slvinv*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // nll -= log(Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix().determinant()) - 0.5*((atomic::matinv(Slv(q))*(Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix()*Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix().transpose())).diagonal().sum()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(atomic::matinv(Slv(q))*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // nll -= 0.5*(log((Alvm.col(i).matrix()*Alvm.col(i).matrix().transpose()).determinant()) - (Slv(q).inverse()*(Alvm.col(i).matrix()*Alvm.col(i).matrix().transpose())).diagonal().sum()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(Slv(q).inverse()*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
            }
            
            // Report
            // REPORT(Alvm);
          }
          
          REPORT(Alvm);
        } else if(num_corlv>1){
          //   // Kron A=AQ*Alvm
          matrix<Type> Alvm(times*nu,times*nu);
          Alvm.fill(0.0);
          //
          //   // Variational covariance
          //   //diagonal
          for(i=0; i<nu; i++){
            for (d=0; d<(times); d++){
              Alvm(i*times+d,i*times+d)=exp(Au(i*times+d));
            }
          }
          
          //   //reduced rank cov
          int k=0;
          arank = NN.rows();
          if(Au.size()>(nu*times+num_corlv*(num_corlv+1)/2)) {
            if(Astruc == 4) {
              for(i=0; i<nu; i++){
                for (int r=0; r<(arank); r++){
                  Alvm(i*times+NN(r,0)-1,i*times+NN(r,1)-1)=Au(nu*times+k);
                  k++;
                }
              }
            } else if(Astruc == 3){
              for(i=0; i<nu; i++){
                for (d=0; (d<times); d++){
                  for (int r=d+1; r<(times); r++){
                    Alvm(i*times+r,i*times+d)=Au(nu*times+k);
                    k++;
                  }
                }
              }
            }
          }
          //
          for (d=0; d<num_corlv; d++){
            AQ(d,d)=exp(Au(nu*times+k));
            k++;
            for (int r=d+1; r<(num_corlv); r++){
              AQ(r,d)=Au(nu*times+k);
              k++;
            }
          }
          // Alvm *= Alvm.transpose();
          // AQ *= AQ.transpose();
          
          for (j=0; j<p;j++){
            // for (i=0; i<(nu*times);i++){
            //   cQ(i,j) += 0.5*(Alvm.row(i)*Alvm.row(i).transpose()).sum()*((newlamCor.col(j).transpose()*(AQ*AQ.transpose()))*newlamCor.col(j)).sum();
            // }
            cQ.col(j) = cQ.col(j) + 0.5*(Alvm*Alvm.transpose()).diagonal().matrix()*((newlamCor.col(j).transpose()*(AQ*AQ.transpose()))*newlamCor.col(j));
            // // cQ.col(j) = cQ.col(j) + 0.5*Alvm.diagonal().matrix()*((newlamCor.col(j).transpose()*AQ)*newlamCor.col(j));
          }
          nll -= num_corlv*log(Alvm.determinant()) + times*nu*log(AQ.determinant()) + 0.5*num_corlv*times*nu;
          // nll -= 0.5*num_corlv*log(Alvm.determinant()) + 0.5*times*nu*log(AQ.determinant()) + 0.5*num_corlv*times*nu;
          //
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            // Slv.fill(0.0);
            
            // Define covariance matrix
            if(cstruc==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else {
              DiSc.fill(0.0);
              for(int j=0; j<dc.cols(); j++){
                DiSc(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled = dc*DiSc;
              if(cstruc==2){// exp decaying
                Slv(q) = gllvm::corExp(Type(1), Type(0), times, dc_scaled);
                // Slv(q) = gllvm::corExp(Type(1), (rho_lvc(q,0)), times, DistM);
              } else if(cstruc==4) {// matern
                Slv(q) = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), times, dc_scaled);
                // Slv(q) = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), times, DistM);
              }
            }
            
            nll -= - 0.5*nu*log(Slv(q).determinant());
            Slvinv = atomic::matinv(Slv(q));
            
            for (i=0; i<nu; i++) {
              nll -=  0.5*(- (AQ.row(q)*AQ.row(q).transpose()).sum()*(Slvinv*(Alvm.block(i*times,i*times,times,times).matrix()*Alvm.block(i*times,i*times,times,times).matrix().transpose())).diagonal().sum() - ((ucopy.block(i*times,q,times,1).matrix()).transpose()*(Slvinv*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // nll -=  0.5*(- AQ(q,q)*(atomic::matinv(Slv(q))*Alvm.block(i*times,i*times,times,times).matrix()).diagonal().sum()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(atomic::matinv(Slv(q))*(ucopy.block(i*times,q,times,1).matrix()))).sum());
            }
            
          }
          REPORT(Alvm);
        }
        
        REPORT(Slv);
        // REPORT(Alvm);
      }
      REPORT(AQ);
      REPORT(nu);
    }
    
    
    if(nlvr>0){
      //constrained ordination terms
      if((num_lv_c>0) && (random(2)<1)){
        //predictor coefficients for constrained ordination
        if((random(0)>0) && (n == nr)){
          //first column are zeros in case of random intercept
          b_lv2.middleCols(1,num_lv_c) = b_lv.leftCols(num_lv_c);
          
        }else{
          b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        }
        
        eta += x_lv*b_lv2*newlam;
        //quadratic term for constrained ordination
        if(quadratic>0){
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose();
            }
          }
        }
      }
      lam += u*newlam;
      
      // Update cQ for non quadratic latent variable model and
      // also takes this route if there are quadratic constrained LVs with random row-effect
      if((quadratic < 1) || ( (nlvr==1 && random(2)<0 && num_RR>0))){
        
        //Binomial, Gaussian, Ordinal
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            cQ(i,j) += 0.5*((newlam.col(j)).transpose()*((A.col(i).matrix()*A.col(i).matrix().transpose()).matrix()*newlam.col(j))).sum();
          }
        }
        eta += lam;
      }
      // do not take this route not with quadratic model, constrained LVs and random row-effects.
      if(( (quadratic>0 && nlvr > 0) && (num_lv+num_lv_c)>0 ) || ( quadratic>0 && random(2) > 0 )){
        matrix <Type> Acov(nlvr,nlvr);
        //quadratic model approximation
        //Poisson
        if(family==0){
          matrix <Type> B(nlvr,nlvr);
          matrix <Type> v(nlvr,1);
          for (int i=0; i<n; i++) {
            Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
            matrix <Type> Ql = (A.col(i).matrix()).inverse();
            matrix <Type> Q = Ql*Ql.transpose();
            for (int j=0; j<p;j++){
              B = (2*D.col(j).matrix()+Q);
              if((random(2)>0) || ((num_lv>0)&(num_lv_c==0))){
                v = (newlam.col(j)+Q*u.row(i).transpose());
              }else if(random(2)<1){
                //extra term for concurrent ordination
                v = (newlam.col(j)+Q*u.row(i).transpose() - 2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose());
              }
              
              Type detB = atomic::logdet(B);
              Type detA = ((vector <Type> (A.col(i).matrix().diagonal())).log()).sum(); //log-determinant of cholesky
              e_eta(i,j) += exp(cQ(i,j) + eta(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()-detB)-detA); //add all the other stuff to the quadratic approximation
              eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
              
              if((num_lv_c>0) && (random(2)<1)){
                eta(i,j) -= 2*u.row(i)*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
        // //NB, gamma, exponential
        if((family==1)|(family==4)|(family==8)){
          matrix <Type> B(nlvr,nlvr);
          matrix <Type> v(nlvr,1);
          for (int i=0; i<n; i++) {
            Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
            matrix <Type> Ql = (A.col(i).matrix()).inverse();
            matrix <Type> Q = Ql*Ql.transpose();
            for (int j=0; j<p;j++){
              B = (-2*D.col(j).matrix()+Q);
              if((random(2)>0) || ((num_lv>0)&(num_lv_c==0))){
                v = (-newlam.col(j)+Q*u.row(i).transpose());
              }else if(random(2)<1){
                //extra term for concurrent ordination
                v = (-newlam.col(j)+Q*u.row(i).transpose() + 2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose());
              }
              
              Type detB = log((B.llt().matrixL()).determinant());//required like this due to potential negative semi-definiteness
              Type detA = ((vector <Type> (A.col(i).matrix().diagonal())).log()).sum(); //log-determinant of cholesky
              e_eta(i,j) += exp(-eta(i,j) - cQ(i,j)+0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value())-detA-detB);
              eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
              if((num_lv_c>0) && (random(2)<1)){
                eta(i,j) -= 2*u.row(i)*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
        //Binomial, Gaussian, Ordinal
        if((family==2)|(family==3)|(family==7)){
          for (int i=0; i<n; i++) {
            Acov = (A.col(i).matrix()*A.col(i).matrix().transpose()).matrix();
            for (int j=0; j<p;j++){
              if((random(2)>0) || ((num_lv>0)&(num_lv_c==0))){
                cQ(i,j) += 0.5*(newlam.col(j)*newlam.col(j).transpose()*Acov).trace() + (D.col(j).matrix()*Acov*D.col(j).matrix()*Acov).trace() +2*(u.row(i)*D.col(j).matrix()*Acov*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*Acov*newlam.col(j)).value();
              }else if(random(2)<1){
                //extra terms for concurrent ordination
                cQ(i,j) += 0.5*((newlam.col(j)-2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose())*(newlam.col(j)-2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose()).transpose()*Acov).trace() + (D.col(j).matrix()*Acov*D.col(j).matrix()*Acov).trace() +2*(u.row(i)*D.col(j).matrix()*Acov*D.col(j).matrix()*(u.row(i)).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*Acov*(newlam.col(j)-2*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose())).value();
              }
              eta(i,j) += lam(i,j) - (u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value() - (D.col(j).matrix()*Acov).trace();
              if((num_lv_c>0) && (random(2)<1)){
                eta(i,j) -= 2*u.row(i)*D.col(j).matrix()*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
      }
    }
    // REPORT(eta);
    // REPORT(cQ);
    
    if(family==0){//poisson
      if((quadratic < 1) || ( (quadratic > 0 && (num_lv+num_lv_c)<1 && nlvr >0) )){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= dpois(y(i,j), exp(eta(i,j)+cQ(i,j)), true)-y(i,j)*cQ(i,j);
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= y(i,j)*eta(i,j) - e_eta(i,j) - lfactorial(y(i,j));
          }
        }
      }
    } else if((family == 1) & (method<1)){//NB VA
      if((quadratic < 1) || ( (quadratic > 0 && (num_lv+num_lv_c)<1 && nlvr >0) )){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            // nll -= Type(gllvm::dnegbinva(y(i,j), eta(i,j), iphi(j), cQ(i,j)));
            nll -= y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*cQ(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
          }
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= -iphi(j)*eta(i,j) -(y(i,j)+iphi(j))*log(1+iphi(j)*e_eta(i,j))+ lgamma(y(i,j)+iphi(j))+ iphi(j)*log(iphi(j)) -lgamma(iphi(j)) -lfactorial(y(i,j));
            //log(1+phi*e_eta) = log(phi+1/e_eta)+log(e_eta)
          }
        }
      }
      
    } else if ((family == 1) & (method>1)) { // NB EVA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
          nll += (((iphi(j)+y(i,j)) / (iphi(j)+exp(eta(i,j)))) * exp(eta(i,j)) - ((iphi(j)+y(i,j))*pow(iphi(j)+exp(eta(i,j)),-2))*pow(exp(eta(i,j)),2)) * cQ(i,j);
          
          // nll += gllvm::nb_Hess(y(i,j), eta(i,j), iphi(j)) * cQ(i,j);
          // nll -= lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) - lgamma(y(i,j)+1) + y(i,j)*eta(i,j) + iphi(j)*log(iphi(j))-(y(i,j)+iphi(j))*log(exp(eta(i,j))+iphi(j));
          // nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
          // nll += (((iphi(j)+y(i,j)) / (iphi(j)+exp(eta(i,j)))) * exp(eta(i,j)) - ((iphi(j)+y(i,j))*pow(iphi(j)+exp(eta(i,j)),-2))*pow(exp(eta(i,j)),2)) * cQ(i,j);
        }
      }
    } else if((family == 2) & (method<1)) {//binomial probit VA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
          nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j)))) - cQ(i,j);
        }
      }
    } else if ((family == 2) & (method>1)) { // Binomial EVA
      if (extra(0) == 0) { // logit
        for (int i=0; i<n; i++) {
          for (int j=0; j<p; j++) {
            // nll -= gllvm::dbinom_logit_eva(y(i,j), eta(i,j), cQ(i,j));
            
            Type mu = 0.0;
            Type mu_prime = 0.0;
            
            CppAD::vector<Type> z(4);
            z[0] = eta(i,j);
            z[1] = 0;
            z[2] = 1/(1+exp(-z[0]));
            z[3] = exp(z[0])/(exp(z[0])+1);
            
            mu = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
            mu_prime = mu * (1-mu);
            nll -= y(i,j) * eta(i,j) + log(1-mu);
            nll += mu_prime*cQ(i,j);
          }
        }
      } else if (extra(0) == 1) { // probit
        for (int i=0; i<n; i++) {
          for (int j=0; j<p; j++) {
            Type etaP = pnorm_approx(Type(eta(i,j)));   //pnorm funktion approksimaatio
            nll -= y(i,j)*log(etaP) + (1-y(i,j))*log(1-etaP); //
            Type etaD =  dnorm(Type(eta(i,j)), Type(0), Type(1), true);   // log normal density evaluated at eta(i,j)
            nll -= ((y(i,j)*(etaP*exp(etaD)*(-eta(i,j))-pow(exp(etaD),2))*pow(1-etaP,2) + (1-y(i,j))*((1-etaP)*exp(etaD)*eta(i,j)-pow(exp(etaD),2))*pow(etaP,2) )/(etaP*etaP*(etaP*etaP-2*etaP+1)))*cQ(i,j); //Tää toimii ok tähän etaD = (log=true)
          }
        }
      }
    } else if(family==3) {//gaussian
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          nll -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j))) - log(M_PI)/2;
        }
      }
    } else if(family==4) {//gamma
      if((quadratic < 1) || ( (quadratic > 0 && (num_lv+num_lv_c)<1 && nlvr >0) )){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));
          }
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -=  ( -eta(i,j) - e_eta(i,j)*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));
          }
        }
      }
      
    } else if(family==5){ // Tweedie EVA
      Type v = extra(0);
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          // Tweedie log-likelihood:
          nll -= dtweedie(y(i,j), exp(eta(i,j)), iphi(j), v, true);
          if (y(i,j) == 0) {
            // Hessian-trace part:
            nll += (1/iphi(j)) * (2-v)*exp(2*eta(i,j))*exp(-v*eta(i,j)) * cQ(i,j);
          } else if (y(i,j) > 0) {
            nll -= (1/iphi(j)) * (y(i,j)*(1-v)*exp((1-v)*eta(i,j)) - (2-v)*exp((2-v)*eta(i,j))) * cQ(i,j);
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
            nll -= log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)));
          }else if(y(i,j)==ymaxj){
            //maximum category
            int idx = ymaxj-2;
            nll -= log(1 - pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1)));
          }else if(ymaxj>2){
            for (int l=2; l<ymaxj; l++) {
              if((y(i,j)==l) && (l != ymaxj)){
                nll -= log(pnorm(zetanew(j,l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,l-2)-eta(i,j), Type(0), Type(1)));
              }
            }
          }
          
          nll += cQ(i,j);
          //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));//
        }
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
            nll -= log(pnorm(zetanew(0) - eta(i,j), Type(0), Type(1)));
          }else if(y(i,j)==ymax){
            //maximum category
            int idx = ymax-2;
            nll -= log(1 - pnorm(zetanew(idx) - eta(i,j), Type(0), Type(1)));
          }else if(ymax>2){
            for (int l=2; l<ymax; l++) {
              if((y(i,j)==l) && (l != ymax)){
                nll -= log(pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1)));
              }
            }
          }
          nll += cQ(i,j);
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      }
    } else if(family==8) {// exp dist
      if((quadratic < 1) || ( (quadratic > 0 && (num_lv+num_lv_c)<1 && nlvr >0) )){
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) );
          }
        }
      }else{
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= ( -eta(i,j) - e_eta(i,j)*y(i,j) );
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
          nll -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
          nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
          nll -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
          //
        }
      }
    }
    // nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
    
  } else {
    if(random(2)>0){
      // REPORT(Sigmab_lv); //!!!!
      //MVNORM_t<Type> mvnorm(Sigmab_lv);
      if(sbl3==Klv){
        for (int klv=0; klv<Klv; klv++) {
          nll += MVNORM(Sigmab_lv.col(klv).matrix())(b_lv.row(klv));
        }
      }else{
        for (int q=0; q<(num_lv_c+num_RR); q++) {
          nll += MVNORM(Sigmab_lv.col(q).matrix())(b_lv.col(q));
        }
      }
      
    }
    
    // REPORT(ucopy);
    // REPORT(num_corlv);
    // REPORT(nu);
    // REPORT(dr);
    // REPORT(cstruc);
    
    matrix<Type> etaH(n,p); 
    etaH.fill(0.0);
    
    //For fixed-effects RRR with and without quadratic term
    if(num_RR>0){
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      if(quadratic>0){
        matrix <Type> D_RR(num_RR,num_RR);
        D_RR.fill(0.0);
        if(lambda2.cols()==1){
          for (int d=0; d<num_RR;d++){
            D_RR(d,d) = fabs(lambda2(d,0));
          }
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
          
        }else{
          for (int j=0; j<p;j++){
            for (int d=0; d<num_RR;d++){
              D_RR(d,d) = fabs(lambda2(d,j));
            }
            
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
        }
        
      }
    }
    
    // Laplace approximation
    
    // add offset to lin. predictor 
    eta += offset;
    if(random(0)==0){
      eta += r0*xr;
    }
    
    
    // Include random slopes if random(1)>0
    if(random(1)>0){
      vector<Type> sdsv = exp(sigmaB);
      density::UNSTRUCTURED_CORR_t<Type> neg_log_MVN(sigmaij);
      for (int j=0; j<p;j++){
        nll += VECSCALE(neg_log_MVN,sdsv)(vector<Type>(Br.col(j)));
      }
      eta += xb*Br;
    }
    
    //latent variables and random site effects (r_i,u_i) from N(0,Cu)
    if(nlvr>0){
      MVNORM_t<Type> mvnorm(Cu);
      for (int i=0; i<n; i++) {
        nll += mvnorm(u.row(i));
      }
      //variances of LVs
      u *= Delta;
      
      if(num_lv_c>0){
        if((random(0)>0) && (n == nr)){
          //first column are zeros in case of random intercept
          b_lv2.middleCols(1,num_lv_c) = b_lv.leftCols(num_lv_c);
          
        }else{
          b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        }  
        eta += x_lv*b_lv2*newlam;
      }
      // add LV term to lin. predictor 
      lam += u*newlam;
      eta += lam;
      if(family==10){
        // etaH += lam;
        etaH += ucopy*thetaH;
      }
    }
    
    
    // Structured Row/Site effects
    if(((random(0)>0) & (nlvr==(num_lv+num_lv_c))) & (rstruc>0)){
      int i,j;
      // Group specific random row effects:
      if(rstruc == 1){
        if(cstruc ==0){
          matrix<Type> Sr(1,1);
          Sr(0,0) = sigma*sigma;
          MVNORM_t<Type> mvnorm(Sr);
          for (int i=0; i<nr; i++) {
            nll += mvnorm(r0.row(i));
          }
        } else {
          // group specific random row effects, which are correlated between groups
          matrix<Type> Sr(nr,nr);
          // if(cstruc==1){// AR1 covariance
          //   Type rho = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          //   for (d=0;d<nr;d++) {
          //     Sr(d,d)=sigma*sigma;
          //     for (j=0;j<d;j++){
          //       Sr(d,j)=sigma*pow(rho,(d-j))*sigma;       // ar1 correlation
          //       Sr(j,d)=Sr(d,j);
          //     }
          //   }
          // } else if(cstruc==2){// exp decaying
          //   Type alf = exp(log_sigma(1));
          //   for (d=0;d<nr;d++) {
          //     Sr(d,d)=sigma*sigma;
          //     for (j=0;j<d;j++){
          //       Sr(d,j)=sigma*exp(-sqrt(((dc.row(d)-dc.row(j))*(dc.row(d)-dc.row(j)).transpose()).sum())/alf)*sigma;
          //       Sr(j,d)=Sr(d,j);
          //     }
          //   }
          // } else {// Compound Symm  if(cstruc==3)
          //   Type rhob = log_sigma(1) / sqrt(1.0 + pow(log_sigma(1), 2));
          //   for (d=0;d<nr;d++) {
          //     Sr(d,d)=sigma*sigma;
          //     for (j=0;j<d;j++){
          //       Sr(d,j)=sigma*rhob*sigma;
          //       Sr(j,d)=Sr(d,j);
          //     }
          //   }
          // }
          // Define covariance matrix
          if(cstruc==1){// AR1 covariance
            Sr = gllvm::corAR1(sigma, log_sigma(1), nr);
          } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
            Sr = gllvm::corCS(sigma, log_sigma(1), nr);
          } else {
            DiSc.fill(0.0);
            for(int j=0; j<dc.cols(); j++){
              DiSc(j,j) += 1/exp(log_sigma(1+j));
            }
            dc_scaled = dc*DiSc;
            if(cstruc==2){// exp decaying
              Sr = gllvm::corExp(sigma, Type(0), nr, dc_scaled);
              // Sr = gllvm::corExp(sigma, (log_sigma(1)), nr, DistM);
            } else if(cstruc==4) {// matern
              Sr = gllvm::corMatern(sigma, Type(0), log_sigma(dc.cols()+1), nr, dc_scaled);
              // Sr = gllvm::corMatern(sigma, log_sigma(1), log_sigma(2), nr, DistM);
            }
          }
          MVNORM_t<Type> mvnorm(Sr);
          nll += mvnorm(r0.col(0));
        }
        
        for (int j=0; j<p;j++){
          eta.col(j) = eta.col(j) + dr*r0;
        }
      } else {
        // site specific random row effects, which are correlated within groups
        
        // Define covariance matrix
        matrix<Type> Sr(times,times);
        // Define covariance matrix
        if(cstruc==1){// AR1 covariance
          Sr = gllvm::corAR1(sigma, log_sigma(1), times);
        } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
          Sr = gllvm::corCS(sigma, log_sigma(1), times);
        } else {
          DiSc.fill(0.0);
          for(int j=0; j<dc.cols(); j++){
            DiSc(j,j) += 1/exp(log_sigma(1+j));
          }
          dc_scaled = dc*DiSc;
          if(cstruc==2){// exp decaying
            Sr = gllvm::corExp(sigma, (Type(0)), times, dc_scaled);
            // Sr = gllvm::corExp(sigma, (log_sigma(1)), times, DistM);
          } else if(cstruc==4) {// matern
            Sr = gllvm::corMatern(sigma, Type(0), log_sigma(dc.cols()+1), times, dc_scaled);
            // Sr = gllvm::corMatern(sigma, log_sigma(1), log_sigma(2), times, DistM);
          }
        }
        
        for (j=0; j<p;j++){
          eta.col(j).array() += r0.array();
        }
        MVNORM_t<Type> mvnorm(Sr);
        r0.resize(times, nr);
        for (i=0; i<nr; i++) {
          nll += mvnorm(vector <Type> (r0.col(i)));
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
    
    // Correlated LVs
    if(num_corlv>0) {
      int i;
      if(ucopy.rows() == nu){
        eta += (dr*ucopy)*newlamCor;
        if(family==10){ // betaH
          etaH += (dr*ucopy)*thetaH;
        }
        
        // group specific lvs
        if(cstruc==0){// no covariance
          matrix<Type> Slv(num_corlv,num_corlv);
          Slv.fill(0.0);
          Slv.diagonal().fill(1.0);
          MVNORM_t<Type> mvnorm(Slv);
          for (int i=0; i<nu; i++) {
            nll += mvnorm(ucopy.row(i));
          }
          REPORT(Slv);
        } else {
          
          matrix<Type> Slv(nu,nu);
          
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated between groups
            Slv.fill(0.0);
            
            if(cstruc==1){// AR1 covariance
              Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
            } else {
              DiSc.fill(0.0);
              for(int j=0; j<dc.cols(); j++){
                DiSc(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled = dc*DiSc;
              if(cstruc==2){// exp decaying
                Slv = gllvm::corExp(Type(1), Type(0), nu, dc_scaled);
                // Slv = gllvm::corExp(Type(1), (rho_lvc(q,0)), nu, DistM);
              } else if(cstruc==4) {// matern
                Slv = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), nu, dc_scaled);
                // Slv = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), nu, DistM);
              }
            }
            
            MVNORM_t<Type> mvnormS1(Slv);
            nll += mvnormS1(ucopy.col(q));
            
          }
          REPORT(Slv);
        }
      } else {
        
        matrix<Type> Slv(times,times);
        eta += ucopy*newlamCor;
        if(family==10){// betaH
          etaH += ucopy*thetaH;
        }
        for(int q=0; q<num_corlv; q++){
          // site specific LVs, which are correlated within groups
          Slv.fill(0.0);
          // Define covariance matrix
          if(cstruc==1){// AR1 covariance
            Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
          } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
            Slv = gllvm::corCS(Type(1), rho_lvc(q,0), times);
          } else {
            DiSc.fill(0.0);
            for(int j=0; j<dc.cols(); j++){
              DiSc(j,j) += 1/exp(rho_lvc(q,j));
            }
            dc_scaled = dc*DiSc;
            if(cstruc==2){// exp decaying
              Slv = gllvm::corExp(Type(1), Type(0), times, dc_scaled);
              // Slv = gllvm::corExp(Type(1), (rho_lvc(q,0)), times, DistM);
            } else if(cstruc==4) {// matern
              Slv = gllvm::corMatern(Type(1), Type(0), rho_lvc(q,dc.cols()), times, dc_scaled);
              // Slv = gllvm::corMatern(Type(1), rho_lvc(q,0), rho_lvc(q,1), times, DistM);
            }
          }
          
          MVNORM_t<Type> mvnormS2(Slv);
          
          for (i=0; i<nu; i++) {
            nll += mvnormS2(ucopy.block(i*times,q,times,1));
            
          }
          
        }
        REPORT(Slv);
      }
      // REPORT(nu);
    }
    
    
    if(model<1){
      // gllvm.TMB.R
      if(family==10){
        etaH += x*bH;
      }
      eta += x*b;
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          mu(i,j) = exp(eta(i,j));
        }
      }
      
    } else {
      // Fourth corner model, TMBtrait.R
      if(family==10){
        matrix<Type> eta1h=x*bH;
        eta1h.resize(n, p);
        etaH += eta1h;
      }
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
          nll -= dpois(y(i,j), exp(eta(i,j)), true);
        }
      }
    } else if(family==1){//negative.binomial family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          nll -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
        }
      }} else if(family==2) {//binomial family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));
          }
        }
      } else if(family==3){//gaussian family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dnorm(y(i,j), eta(i,j), iphi(j), true); 
          }
        }
      } else if(family==4){//gamma family
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dgamma(y(i,j), iphi(j), exp(eta(i,j))/iphi(j), true); 
          }
        }
      } else if(family==5){//tweedie familyF
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dtweedie(y(i,j), exp(eta(i,j)),iphi(j),extra(0), true); 
          }
        }
      } else if(family==6) {//zero-infl-poisson
        iphi=iphi/(1+iphi);
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= dzipois(y(i,j), exp(eta(i,j)),iphi(j), true); 
          }
        }
      } else if(family==8) {// exponential family
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            nll -= dexp(y(i,j), exp(-eta(i,j)), true);  // (-eta(i,j) - exp(-eta(i,j))*y(i,j) );
          }
        }
      } else if(family==9) {// beta family
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {mu(i,j) = pnorm(eta(i,j));}
            nll -= dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
          }
        }
      } else if(family==10) {// beta family
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            if(extra(0)<1) {
              etaH(i,j) = exp(etaH(i,j))/(exp(etaH(i,j))+1);
              mu(i,j) = mu(i,j)/(mu(i,j)+1);
            } else {
              etaH(i,j) = pnorm(etaH(i,j));
              mu(i,j) = pnorm(eta(i,j));
            }
            if (y(i,j) == 0) {
              nll -= log(1-mu(i,j));
            } else{
              nll -= log(mu(i,j)) + dbeta(y(i,j), Type(etaH(i,j)*iphi(j)), Type((1-etaH(i,j))*iphi(j)), 1);
            }
          }
        }
        // REPORT(mu);
        // REPORT(etaH);
      }
  }
  
  return nll;
}
