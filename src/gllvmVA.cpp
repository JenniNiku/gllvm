#include <TMB.hpp>
#include "distrib.h"
#define TMB_LIB_INIT R_init_gllvmVA
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
  // if(((num_corlv>0) || (((random(0)>0) && (nlvr==(num_lv+num_lv_c))) && (rstruc>0))) && ((cstruc==2) || (cstruc>3))){
  //   matrix<Type> DiSc(dc.cols(),dc.cols());
  //   DiSc.setZero();
  // 
  //   for(int j=0; j<dc.cols(); j++){
  //     DiSc(j,j) += 1/exp(2*scaledc(j));
  //     // dc.col(j) *= 1/exp(scaledc(j));
  //   }
  //   // sigma_lvc(0,0) = 0;
  // 
  //   DistM.setZero();
  //   for (int d=0;d<dc.rows();d++) {
  //     for (int j=0;j<d;j++){
  //       DistM(d,j)=sqrt( ((dc.row(d)-dc.row(j))*DiSc*(dc.row(d)-dc.row(j)).transpose()).sum() ); // + extra(2);
  //   //     DistM(j,d)=DistM(d,j);
  //     }
  //   }
  // }
  
  // if row params are in the same form as LVs, let's put them together
  if((random(0)>0) && (n == nr)){
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
  eta.setZero();
  matrix<Type> lam(n,p);
  matrix<Type> Cu(nlvr,nlvr); 
  Cu.setZero();  
  
  Type nll = 0; // initial value of log-likelihood
  
  matrix<Type> b_lv2(x_lv.cols(),nlvr);
  matrix<Type> RRgamma(num_RR,p);
  b_lv2.setZero();
  RRgamma.setZero();
  
  matrix <Type> Delta(nlvr,nlvr);
  Delta.setZero();
  
  
  vector<matrix<Type>> D(p);
  
  if( (quadratic>0) && (method!=1)){
    Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Dmat(nlvr+num_RR);
    Dmat.setZero();
    for (int j=0; j<p; j++){
      D(j) = Dmat;
    }
  }else if(quadratic>0){
    Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Dmat(nlvr);
    Dmat.setZero();
    for (int j=0; j<p; j++){
      D(j) = Dmat;
    }
  }
  
  matrix <Type> newlam(nlvr,p);
  newlam.setZero();  
  
  //K*K*d or d*d*K
  int sbl12 = 0;
  int sbl3 = 0;
  if((sigmab_lv.size()==Klv)||(sigmab_lv.size()==Type(1))){
    sbl12 = num_lv_c + num_RR;
    sbl3 = Klv;
  }else if(sigmab_lv.size()==(num_lv_c+num_RR)){
    sbl12 = Klv;
    sbl3 = num_lv_c + num_RR;
  }
  
  vector<matrix<Type>> Sigmab_lv(sbl3);
  if(random(2)>0){
    sigmab_lv = exp(sigmab_lv);
    sigmab_lv *= sigmab_lv;
    
    if(sigmab_lv.size()>Type(1)){//randomB=="LV", Sigma_q = sigma_q I_klv
      Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Sigmab_lvtemp(sbl12);
      for (int q=0; q<sbl3; q++){
        Sigmab_lv(q) = Sigmab_lvtemp;
        Sigmab_lv(q).diagonal().array() = sigmab_lv(q);
      }
    }else if(sigmab_lv.size()==Type(1)){
      Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Sigmab_lvtemp(sbl12);
      for (int klv=0; klv<Klv; klv++){
        Sigmab_lv(klv) = Sigmab_lvtemp;
        Sigmab_lv(klv).diagonal().array() = sigmab_lv(0);
      }
    }
  }
  
  if((nlvr>0)||(num_RR>0)){
    
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
    //To create lambda as matrix Upper triangle
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
  matrix<Type> newlamCor;
  // matrix <Type> Delta_clv(num_corlv,num_corlv);
  if((num_corlv)>0){
    newlamCor = matrix <Type> (num_corlv,p);
    //To create lambda as matrix Upper triangle
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
  
    // Variational approximation
    //quadratic coefficients for ordination
    //if random rows, add quadratic coefficients for num_RR to D otherwise
    //they go into D_RR below
    //The ordering here is num_lv_c-num_lv-num_RR so that the code works for
    //fixed-effects B and random effects B
    //The order we need to pick them from lambda2 is 
    //num_lv_c-num_RR-num_lv however, to ensure everything on the R-side works
    if((quadratic>0) && ((num_lv+num_lv_c+num_RR*random(2))>0)){
      if(nlvr>(num_lv+num_lv_c)){
        if(num_lv_c>0){
          if(lambda2.cols()==1){
            for (int j=0; j<p; j++){
              for (int q=1; q<(num_lv_c+1); q++){
                D(j).diagonal()(q) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=1; q<(num_lv_c+1); q++){
                D(j).diagonal()(q) = fabs(lambda2(q-1,j)); //full quadratic model
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
                D(j).diagonal()(q+num_lv) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1); q<(num_lv_c+1+num_RR); q++){
                D(j).diagonal()(q+num_lv) = fabs(lambda2(q-1,j)); //full quadratic model
              }
            } 
          }
        }
        if(num_lv>0){
          if(lambda2.cols()==1){
            //make sure that num_lv is taken from the middle even with num_RR
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1+num_RR); q<(num_lv_c+1+num_RR+num_lv); q++){
                D(j).diagonal()(q-num_RR) = fabs(lambda2(q-1,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+1+num_RR); q<(num_lv_c+1+num_RR+num_lv); q++){
                D(j).diagonal()(q-num_RR) = fabs(lambda2(q-1,j)); //full quadratic model
              }
            } 
          }
        }
      }else{
        if(num_lv_c>0){
          if(lambda2.cols()==1){
            for (int j=0; j<p; j++){
              for (int q=0; q<num_lv_c; q++){
                D(j).diagonal()(q) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=0; q<num_lv_c; q++){
                D(j).diagonal()(q) = fabs(lambda2(q,j)); //full quadratic model
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
                D(j).diagonal()(q+num_lv) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=num_lv_c; q<(num_lv_c+num_RR); q++){
                D(j).diagonal()(q+num_lv) = fabs(lambda2(q,j)); //full quadratic model
              }
            } 
          }
        }
        if(num_lv>0){
          if(lambda2.cols()==1){
            //make sure that num_lv is taken from the middle even with num_RR
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+num_RR); q<(num_lv_c+num_RR+num_lv); q++){
                D(j).diagonal()(q-num_RR) = fabs(lambda2(q,0)); //common tolerances model
              }
            } 
          }else{
            for (int j=0; j<p; j++){
              for (int q=(num_lv_c+num_RR); q<(num_lv_c+num_RR+num_lv); q++){
                D(j).diagonal()(q-num_RR) = fabs(lambda2(q,j)); //full quadratic model
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
    cQ.setZero();
    vector<matrix<Type>> A(n);
    
    if( (random(2)>0) && (num_RR!=0)){
      for(int i=0; i<n; i++){
        A(i).resize(nlvr+num_RR,nlvr+num_RR);
        A(i).setZero();
      }
    }else{
      for(int i=0; i<n; i++){
        A(i).resize(nlvr,nlvr);
        A(i).setZero();
      }
    }
    
    // lltOfB.matrixL() = A(0).template triangularView<Lower>;//wouuld be great if we could store A(i) each as a triangular matrix where the upper zeros are ignored
    // Set up variational covariance matrix for LVs 
    if(nlvr>0){
      // Include variational covs of row effects, if structure is same for both
      if(nlvr>(num_lv+num_lv_c)){
        for(int i=0; i<n; i++){
          A(i)(0,0)=exp(lg_Ar(i));
        }
        if(lg_Ar.size()>n){
          for (int r=1; r<nlvr; r++){
            for(int i=0; i<n; i++){
              A(i)(r,0)=lg_Ar(r*n+i);
            }}
        }
      }
      
      
      if((num_lv+num_lv_c)>0){
        // log-Cholesky parametrization for A_i:s
        // don't include num_RR for random slopes, comes in later
        for (int d=0; d<(num_lv+num_lv_c); d++){
          for(int i=0; i<n; i++){
            A(i)(d+(nlvr-num_lv-num_lv_c),d+(nlvr-num_lv-num_lv_c))=exp(Au(d*n+i));
            // A(d,d,i)=exp(Au(d*n+i));
          }
        }
        if(Au.size()>((num_lv+num_lv_c)*n)){
          int k=0;
          for (int c=0; c<(num_lv+num_lv_c); c++){
            for (int r=c+1; r<(num_lv+num_lv_c); r++){
              for(int i=0; i<n; i++){
                A(i)(r+(nlvr-num_lv-num_lv_c),c+(nlvr-num_lv-num_lv_c))=Au((num_lv+num_lv_c)*n+k*n+i);
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
              A(i)(d,0) = 0.0;
            }
          }
        }
      }
      
      // // Add VA terms to logL
      if(random(2)<1){
        //Go this route if no random Bs
        vector <Type> Adiag(nlvr);
        matrix <Type> CuI;
        for(int i=0; i<n; i++){
          Adiag = A(i).diagonal();
          if(nlvr == (num_lv+num_lv_c)) nll -= (Adiag.log()).sum() - 0.5*((A(i)*A(i).transpose()).trace()+(u.row(i)*u.row(i).transpose()).sum());
          if(nlvr>(num_lv+num_lv_c)) {
            CuI = Cu.inverse();//for small matrices use .inverse rather than atomic::matinv
            nll -= (Adiag.log()).sum() - 0.5*(CuI*A(i)*A(i).transpose()).trace()-0.5*((u.row(i)*CuI)*u.row(i).transpose()).sum();
          }
          
          // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
          nll -= 0.5*(nlvr - log(Cu.determinant())*random(0));
        }
        //scale LVs with standard deviations, as well as the VA covariance matrices
        u *= Delta;
        
        for (int i=0; i<n; i++) {
          A(i) *= Delta; 
        }
      }else{
        //Go this route with random Bs (since size of Cu and A.col(i) are then not the same)
        matrix <Type>Atemp(nlvr,nlvr);
        vector <Type> AtempDiag(nlvr);
        matrix <Type> CuI;
        for(int i=0; i<n; i++){
          Atemp = A(i).topLeftCorner(nlvr,nlvr);//to exlcude the 0 rows && columns for num_RR
          AtempDiag = A(i).diagonal();
          if(nlvr == (num_lv+num_lv_c)) nll -= (AtempDiag.log()).sum() - 0.5*((Atemp*Atemp.transpose()).trace()+(u.row(i)*u.row(i).transpose()).sum());
          if(nlvr>(num_lv+num_lv_c)) {
            CuI = Cu.inverse();//for small matrices use .inverse rather than atomic::matinv
            nll -= (AtempDiag.log()).sum() - 0.5*(CuI*Atemp*Atemp.transpose()).trace()-0.5*(u.row(i)*CuI*u.row(i).transpose()).sum(); 
          }
          // log(det(A_i))-sum(trace(Cu^(-1)*A_i))*0.5 sum.diag(A)
          nll -= 0.5*(nlvr - log(Cu.determinant())*random(0));
        }
        
        
        //scale LVs with standard deviations, as well as the VA covariance matrices
        u *= Delta;
        if(num_RR>0){
          Delta.conservativeResize(nlvr+num_RR,nlvr+num_RR);
          for(int d=nlvr; d<(nlvr+num_RR); d++){
            Delta.col(d).setZero();
            Delta.row(d).setZero();
          } 
        }
        
        for (int i=0; i<n; i++) {
          A(i) *= Delta; 
        }
      }
      
    }
    
    //random slopes for constr. ord.
    if((random(2)>0) && ((num_RR+num_lv_c)>0)){
      //resize A, u, D, and add RRGamma to newlam.
      //add columns to u on the right for num_RR with random slopes
      if(num_RR>0){
        u.conservativeResize(n, nlvr + num_RR);  
        //resize and fill newlam, we don't use RRgamma further with random Bs
        //easiest to do is slap RRgamma at the end of newlam
        //this makes the order of newlam, A, u, and D inconsistent with the R-side of things
        //nicer would be to have to same order as in R, but that isn't possible since
        //it requires going down the same route for fixed and random B
        //which would only work with diagonal of 0s in A
        //And that needs to be invertible for the quadratic case, so that is not possible
        newlam.conservativeResize(nlvr+num_RR,p);
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          u.col(d).fill(0.0);
          newlam.row(d).fill(0.0);
        }
        nlvr += num_RR;
        newlam.bottomRows(num_RR) = RRgamma;
      }
      
      // Variational covariance for random slopes
      vector<matrix<Type>> AB_lv(sbl3);
      for(int d=0; d<sbl3; d++){
       AB_lv(d).resize(sbl12,sbl12);
       AB_lv(d).setZero();
      }
      
      for (int q=0; q<(sbl12); q++){
        for(int d=0; d<sbl3; d++){
          AB_lv(d)(q,q)=exp(Ab_lv(q*sbl3+d));
        }
      }
      if(Ab_lv.size()>((sbl12)*sbl3)){
        int k=0;
        for (int c=0; c<(sbl12); c++){
          for (int r=c+1; r<(sbl12); r++){
            for(int d=0; d<sbl3; d++){
              AB_lv(d)(r,c)=Ab_lv((sbl12)*sbl3+k*sbl3+d);
              // Ab(c,r,j)=Ab(r,c,j);
            }
            k++;
          }}
      }
      //VA likelihood parts for random slope
        vector <Type> AB_lvDiag(sbl12);
        vector <Type> Sigmab_lvDiag(sbl12);
        for(int klv=0; klv<sbl3; klv++){
          AB_lvDiag = AB_lv(klv).diagonal();
          Sigmab_lvDiag = Sigmab_lv(klv).diagonal();
          nll -= ((AB_lvDiag.log()).sum() - 0.5*(Sigmab_lv(klv).inverse()*AB_lv(klv)*AB_lv(klv).transpose()).trace()-0.5*(b_lv.row(klv)*Sigmab_lv(klv).inverse()*b_lv.row(klv).transpose()).sum());// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
          nll -= 0.5*(sbl12-(Sigmab_lvDiag.log()).sum());
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
              temp = A(i)*A(i).transpose();
              for(int klv=0; klv<Klv; klv++){
                temp.topLeftCorner(num_lv_c,num_lv_c) += x_lv(i,klv)*x_lv(i,klv)*AB_lv(klv)*AB_lv(klv).transpose();//cholesky of variance block for num_lv_c
              }
              L =  temp.llt().matrixL();//can't do only a part due to potential covariance with num_lv
              A(i) =  L;//have to recompute cholesky of covariance due to summation
            }
          }else if((n == nr) && (random(0)>0)){//if row effects are included in u and A
            for(int i=0; i<n; i++){
              temp = A(i)*A(i).transpose();
              for(int klv=0; klv<Klv; klv++){
                temp.block(1,1,num_lv_c,num_lv_c) += x_lv(i,klv)*x_lv(i,klv)*AB_lv(klv)*AB_lv(klv).transpose();//cholesky of variance block for num_lv_c
              }
              L =  temp.llt().matrixL();//can't do only a part due to potential covariance with num_lv
              A(i) =  L;//have to recompute cholesky of covariance due to summation
            }
          }
        }
        
        if((num_RR>0) && (num_lv_c == 0)){
          // matrix <Type> b_lv3 =  b_lv;//.rightCols(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv;
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0); 
          for(int i=0; i<n; i++){
            L = A(i)*A(i).transpose();
            for(int klv=0; klv<Klv; klv++){
              L.bottomRightCorner(num_RR,num_RR) += x_lv(i,klv)*x_lv(i,klv)*AB_lv(klv)*AB_lv(klv).transpose();//cholesky of variance block for num_lv_c
            }
            
            A(i) =  L.llt().matrixL();//have to recompute cholesky of covariance due to summation
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
            L = A(i)*A(i).transpose();
            temp.setZero();
            for(int klv=0; klv<Klv; klv++){
              temp +=  x_lv(i,klv)*x_lv(i,klv)*AB_lv(klv)*AB_lv(klv).transpose();//num_lv_c variance block
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
            A(i) = L;
          }
          
        }
      }else if(sbl3 == (num_lv_c+num_RR)){//variance per LV
        if(num_RR>0){
          matrix <Type> b_lv3 =  b_lv.rightCols(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv3;
        }
        //much easier, since we assume independence between LVs
        matrix<Type> temp(nlvr,nlvr);
        matrix <Type> L(nlvr,nlvr);
        L.fill(0.0);
        
        if((random(0)<1) || (n != nr)){
          if(num_lv_c>0){
            matrix <Type> b_lv2 =  b_lv.leftCols(num_lv_c);
            u.leftCols(num_lv_c) += x_lv*b_lv2;
          }
          
          if(num_lv_c>0){
            for(int i=0; i<n; i++){
              for(int q=0; q<num_lv_c; q++){
                temp(q,q) = (x_lv.row(i)*AB_lv(q)*AB_lv(q).transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          if(num_RR>0){
            for(int i=0; i<n; i++){
              for(int q=(num_lv_c+num_lv); q<(num_lv_c+num_lv+num_RR); q++){
                temp(q,q) = (x_lv.row(i)*AB_lv(q-num_lv)*AB_lv(q-num_lv).transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          for(int i=0; i<n; i++){
            L = A(i)*A(i).transpose() + temp;
            L = L.llt().matrixL();
            A(i) = L;
          }
          
        }else if((random(0)>1) && (n == nr)){//if row effects are included in u and A
          matrix<Type> temp(nlvr,nlvr);
          // for(int i=0; i<n; i++){
          // temp(i).resize(nlvr,nlvr);
          // temp(i).setZero();
          // }
          
          matrix <Type> L(nlvr,nlvr);
          L.fill(0.0);
          
          if(num_lv_c>0){
            matrix <Type> b_lv2 =  b_lv.leftCols(num_lv_c);
            u.middleCols(1, num_lv_c) += x_lv*b_lv2;
          }
          
          if(num_lv_c>0){
            for(int i=0; i<n; i++){
              for(int q=1; q<(num_lv_c+1); q++){
                temp(q,q) = (x_lv.row(i)*AB_lv(q-1)*AB_lv(q-1).transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          if(num_RR>0){
            for(int i=0; i<n; i++){
              for(int q=(num_lv_c+num_lv+1); q<(num_lv_c+num_lv+num_RR+1); q++){
                temp(q,q) = (x_lv.row(i)*AB_lv(q-num_lv-1)*AB_lv(q-num_lv-1).transpose()*x_lv.row(i).transpose()).sum();
              }
            }
          }
          for(int i=0; i<n; i++){
            L = A(i)*A(i).transpose() + temp;
            L = L.llt().matrixL();
            A(i) = L;
          }
        }
      }
    }
    
    
    // Include random slopes if random(1)>0
    if(random(1)>0){
      matrix<Type> sds(l,l);
      sds.setZero();
      sds.diagonal() = exp(sigmaB);
      matrix<Type> S=sds*density::UNSTRUCTURED_CORR(sigmaij).cov()*sds;

      // Variational covariance for random slopes
      // log-Cholesky parametrization for A_bj:s
      vector<matrix<Type>> Ab(p);
      for(int j=0; j<p; j++){
        Ab(j).resize(l,l);
        Ab(j).setZero();
      }
      
      for (int dl=0; dl<(l); dl++){
        for(int j=0; j<p; j++){
          Ab(j)(dl,dl)=exp(Abb(dl*p+j));
        }
      }
      if(Abb.size()>(l*p)){
        int k=0;
        for (int c=0; c<(l); c++){
          for (int r=c+1; r<(l); r++){
            for(int j=0; j<p; j++){
              Ab(j)(r,c)=Abb(l*p+k*p+j);
              // Ab(c,r,j)=Ab(r,c,j);
            }
            k++;
          }}
      }
      
      /*Calculates the commonly used (1/2) x'_i A_bj x_i
       A is a num.lv x num.lv x n array, theta is p x num.lv matrix*/
      matrix <Type> SI(sigmaij.size(),sigmaij.size());
      vector <Type> AbDiag(l);
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          cQ(i,j) += 0.5*((xb.row(i))*Ab(j)*Ab(j).transpose()*xb.row(i).transpose()).sum();
        }
        AbDiag = Ab(j).diagonal();
        SI = atomic::matinv(S);
        nll -= ((AbDiag.log()).sum() - 0.5*(SI*Ab(j)*Ab(j).transpose()).trace()-0.5*(Br.col(j).transpose()*SI*Br.col(j)).sum());// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
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
    
    matrix <Type> e_eta;
    //components for reduced rank regression terms
    if((num_RR>0) && (random(2)<1)){
      //predictor coefficients RRR.  num_RR comes after num_lv_c
      //Since later dimensions are more likely to have less residual variance
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      
      //quadratic terms for fixed-effects only RRR
      //-num_lv to ensure that we pick num_RR from the middle
      if(quadratic>0){
        Eigen::DiagonalMatrix<Type,Eigen::Dynamic> D_RR(num_RR);
        D_RR.setZero();
        
        //quadratic coefficients for RRR
        if(lambda2.cols()==1){
          for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
            D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,0));
          }
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
          
        }else{
          for (int j=0; j<p;j++){
            for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
              D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,j));
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
              D(j).diagonal()(q) = fabs(lambda2(q-1-num_lv,0)); //common tolerances model
            }
          }
        }else{
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c+1); q<nlvr; q++){
              D(j).diagonal()(q) = fabs(lambda2(q-1-num_lv,j)); //full quadratic model
            }
          }
        }
        
      }else{
        if(lambda2.cols()==1){
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c); q<nlvr; q++){
              D(j).diagonal()(q) = fabs(lambda2(q-num_lv,0)); //common tolerances model
            }
          }
        }else{
          for (int j=0; j<p; j++){
            for (int q=(num_lv+num_lv_c); q<nlvr; q++){
              D(j).diagonal()(q) = fabs(lambda2(q-num_lv,j)); //full quadratic model
            }
          }
        }
      }
    }
    
    
    // Structured Row/Site effects
    if(((random(0)>0) && (nlvr==(num_lv+num_lv_c))) && (rstruc>0)){
      // Group specific random row effects:
      if(rstruc == 1){
        if(cstruc==0){
          for (int j=0; j<p;j++){
            cQ.col(j) += 0.5*(dr*Ar.matrix());
            eta.col(j) += dr*r0;
          }
          for (int i=0; i<nr; i++) {//i<n //!!!
            nll -= 0.5*(1 + log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2) - 2*log(sigma))*random(0); ///(n*p)
          }
        } else {
          // group specific random row effects, which are correlated between groups
          int j,d,r;
          
          matrix<Type> Sr(nr,nr);
          matrix <Type> SRI(nr,nr);
          if(cstruc==1){// AR1 covariance
            Sr = gllvm::corAR1(sigma, log_sigma(1), nr);
          } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
            Sr = gllvm::corCS(sigma, log_sigma(1), nr);
          } else {
            DiSc.setZero();
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
          vector <Type> ArmDiag(nr);
          for (d=0; d<(nr); d++){
            Arm(d,d)=Ar(d);
          }
          
          if((lg_Ar.size()>nr) && (Astruc>0)){ // unstructured Var.cov
            int k=0;
            for (d=0; d<(nr); d++){
              for (r=d+1; r<(nr); r++){
                Arm(r,d)=lg_Ar(nr+k);
                k++;
              }}
          }
          
          for (j=0; j<p;j++){
            cQ.col(j) += 0.5*(dr*(Arm*Arm.transpose()).diagonal().matrix());
            eta.col(j) += dr*r0;
          }
          ArmDiag = Arm.diagonal();
          SRI = atomic::matinv(Sr);
          nll -= (ArmDiag.log()).sum()- 0.5*((SRI*(Arm*Arm.transpose())).trace()-(r0.transpose()*(SRI*r0)).sum());// /(n*p)log(det(Ar_i))-sum(trace(Sr^(-1)Ar_i))*0.5 + ar_i*(Sr^(-1))*ar_i
          
          nll -= 0.5*(nr-log(Sr.determinant()));
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
          DiSc.setZero();
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
        vector<matrix<Type>> Arm(nr);
        for(int i=0; i<nr; i++){
            Arm(i).resize(times,times);
            Arm(i).setZero();
        }
          
        for(i=0; i<nr; i++){
          for (d=0; d<(times); d++){
            Arm(i)(d,d)=Ar(i*times+d);
          }
        }
        if((lg_Ar.size()>(nr*times)) && (Astruc>0)){ // unstructured Var.cov
          int k=0;
          for (d=0; d<(times); d++){
            for (r=d+1; r<(times); r++){
              for(int i=0; i<nr; i++){//i<nr
                Arm(i)(r,d)=lg_Ar(nr*times+k*nr+i);
                // Arm(d,r,i)=Arm(r,d,i);
              }
              k++;
            }}
        }
        
        for (j=0; j<p;j++){
          for (i=0; i<nr; i++) {
            for (d=0; d<(times); d++){
              cQ(i*times + d,j) += 0.5*(Arm(i).row(d)*Arm(i).row(d).transpose()).sum(); //Arm(d,d,i);
            }
          }
          eta.col(j).array() += r0.array();
        }
        r0.resize(times, nr);
        matrix <Type> SRI = atomic::matinv(Sr);
        for (i=0; i<nr; i++) {
          nll -= log(Arm(i).determinant()) + 0.5*( - (SRI*Arm(i)*Arm(i).transpose()).trace()-((r0.col(i).matrix()).transpose()*(SRI*(r0.col(i).matrix()))).sum());
          // log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
        }
        nll -= 0.5*nr*(times - log(Sr.determinant()));
      }
      // eta += dr*r0;
    }
    
    // Correlated LVs
    if(num_corlv>0) { //CorLV
      int i,j,d;
      int arank = 2;
      matrix<Type> AQ(num_corlv,num_corlv);
      AQ.setZero(); AQ.diagonal().fill(1.0);
      
      if(ucopy.rows() == nu){
        
        if(cstruc==0){
          vector<matrix<Type>> Alvm(nu);
          
          for(int d=0; d<nu; d++){
            Alvm(d).resize(num_corlv,num_corlv);
            Alvm(d).setZero();
          }
          
          eta += (dr*ucopy)*newlamCor;
          
          // Variational covariance for row effects
          for (int q=0; q<(num_corlv); q++){
            for (d=0; d<(nu); d++){
              Alvm(d)(q,q)=exp(Au(q*nu+d));
            }
          }
          if((Astruc>0) && (Au.size()>((num_corlv)*nu))){//unstructured cov
            int k=0;
            for (int c=0; c<(num_corlv); c++){
              for (int r=c+1; r<(num_corlv); r++){
                for(d=0; d<nu; d++){
                  Alvm(d)(r,c)=Au(nu*num_corlv+k*nu+d);
                }
                k++;
              }}
          }
          
          
          for (d=0; d<nu; d++) {
            nll -= log(Alvm(d).determinant()) + 0.5*( - (Alvm(d)*Alvm(d).transpose()).trace() - (ucopy.row(d).matrix()*ucopy.row(d).matrix().transpose()).sum());
            // for (d=0; d<nu; d++)
            for (j=0; j<p;j++){
              cQ.col(j) += 0.5*dr.col(d)*((newlamCor.col(j).transpose()*(Alvm(d)*Alvm(d).transpose()))*newlamCor.col(j));
            }
          }
          nll -= 0.5*(nu*num_corlv);
          
        } else {
          vector<matrix<Type> > Slv(num_corlv);
          for(int q=0; q<num_corlv; q++){
            Slv(q).resize(nu,nu);
            Slv(q).setZero();
          }
          
          matrix<Type> Slvinv(nu,nu);
          // matrix<Type> Slv(nu,nu);
          matrix<Type> uq(nu,1);
          eta += (dr*ucopy)*newlamCor;
          
          if(Astruc<3){
            vector<matrix<Type>> Alvm(num_corlv);
            
            for(int d=0; d<nu; d++){
              Alvm(d).resize(nu,nu);
              Alvm(d).setZero();
            }
            
            // matrix<Type> Alvm(nu,nu);
            for(int q=0; q<num_corlv; q++){
              // site specific LVs, which are correlated between groups
              // Slv.setZero();
              uq = ucopy.col(q);
              
              // group specific lvs
              if(cstruc==1){// AR1 covariance
                Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
              } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
                Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
              } else {
                DiSc.setZero();
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
                Alvm(q)(d,d)=exp(Au(q*nu+d));
              }
              
              if((Astruc>0) && (Au.size() > nu*num_corlv)){//reduced rank cov
                // if(Au.size()>(times*nu*num_corlv)){}
                int k=0;
                if(Astruc==1){
                  for (d=0; d<nu; d++){
                    for (int r=d+1; r<(nu); r++){
                      Alvm(q)(r,d)=Au(nu*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                } else if(Astruc==2) {
                  arank = NN.rows();
                  // arank = NN.cols();
                  // for (d=0; (d<nu); d++){
                  for (int r=0; r<(arank); r++){
                    Alvm(q)(NN(r,0)-1,NN(r,1)-1)=Au(nu*num_corlv+k*num_corlv+q);
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
              }
              
              for (j=0; j<p;j++){
                cQ.col(j) += 0.5*pow(newlamCor(q,j),2)*(dr*(Alvm(q)*Alvm(q).transpose()).diagonal().matrix());
              }
              Slvinv = atomic::matinv(Slv(q));
              nll -= log(Alvm(q).determinant()) + 0.5*(- (Slvinv*Alvm(q)*Alvm(q).transpose()).trace()-( uq.transpose()*(Slvinv*uq) ).sum());
              // nll -= log(Alvm.col(q).matrix().determinant()) + 0.5*(- (atomic::matinv(Slv(q))*(Alvm.col(q).matrix()*Alvm.col(q).matrix().transpose())).trace()-( uq.transpose()*(atomic::matinv(Slv(q))*uq) ).sum());
              
              nll -= 0.5*(nu-log(Slv(q).determinant()));
              
            }
          } else if(num_corlv>1){
            matrix<Type> Alvm(nu,nu);
            Alvm.setZero();
            // Kronecker Variational covariance
            for (d=0; d<(nu); d++){
              Alvm(d,d)=exp(Au(d));
            }
            //reduced rank cov
            // if(Au.size()>(times*nu*num_corlv)){}
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
              cQ.col(j) += 0.5*(dr*Alvm.diagonal())*((newlamCor.col(j).transpose()*AQ)*newlamCor.col(j));
            }
            nll -= 0.5*num_corlv*log(Alvm.determinant()) + 0.5*nu*log(AQ.determinant()) + 0.5*num_corlv*nu;
            //
            for(int q=0; q<num_corlv; q++){
              // site specific LVs, which are correlated between groups
              // Slv.setZero();
              uq = ucopy.col(q);
              
              // group specific lvs
              if(cstruc==1){// AR1 covariance
                Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
              } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
                Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
              } else {
                DiSc.setZero();
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
              nll -= 0.5*(- AQ(q,q)*(Slvinv*Alvm).trace()-( uq.transpose()*(Slvinv*uq) ).sum());
              // nll -= 0.5*(- AQ(q,q)*(atomic::matinv(Slv(q))*Alvm).trace()-( uq.transpose()*(atomic::matinv(Slv(q))*uq) ).sum());
              nll -= -0.5*log(Slv(q).determinant());
            }
          }
        }
      } else {
        
        eta += ucopy*newlamCor;
        vector<matrix<Type> > Slv(num_corlv);
        for(int q=0; q<num_corlv; q++){
          Slv(q).resize(times,times);
          Slv(q).setZero();
        }
        // matrix<Type> Slv(times,times);
        matrix<Type> Slvinv(times,times);
        
        // int acol = times;
        // if(Astruc>0){
        //   acol = arank;
        // }
        //
        if(Astruc<3){
          
          vector<matrix<Type>> Alvm(num_corlv);
          
          for(int d=0; d<num_corlv; d++){
            Alvm(d).resize(times*nu,times*nu);
            Alvm(d).setZero();
          }
          
          // array<Type> Alvm(times,times,nu);
          // matrix<Type> uq(times*nu,1);
          
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            
            // Variational covariance for row effects
            //diagonal
            for(i=0; i<nu; i++){
              for (d=0; d<(times); d++){
                Alvm(q)(i*times+d,i*times+d)=exp(Au(q*n+i*times+d));
              }
            }
            
            if((Astruc>0) && (Au.size() > nu*times*num_corlv)){//reduced rank cov
              // if(Au.size()>(times*nu*num_corlv)){}
              int k=0;
              if(Astruc==1){
                for(i=0; i<nu; i++){
                  for (d=0; ((d<arank) && (d<times)); d++){
                    // for (d=0; (d<times); d++){//(num_lv+num_lv_c)*n+k*n+q
                    // for (int r=d; r<(times); r++){
                    // if(r==d){
                    // Alvm(r,d,i)=exp(Au(k*num_corlv+q));
                    // } else {
                    // Alvm(r,d,i)=Au(k*num_corlv+q);
                    // }
                    for (int r=d+1; r<(times); r++){
                      Alvm(q)(i*times+r,i*times+d)=Au(nu*times*num_corlv+k*num_corlv+q);
                      // Alvm(r,d,i)=Au(nu*times*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                }
              } else if(Astruc==2) {
                arank = NN.rows();
                for(i=0; i<nu; i++){
                  for (int r=0; r<(arank); r++){
                    Alvm(q)(i*times+NN(r,0)-1,i*times+NN(r,1)-1)=Au(nu*times*num_corlv+k*num_corlv+q);
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
            }
            
            
            for (j=0; j<p;j++){
              for (i=0; i<(times*nu); i++) {
                cQ(i,j) += 0.5*pow(newlamCor(q,j),2)*(Alvm(q).row(i)*Alvm(q).row(i).transpose()).sum();
              }
            }
            nll -= log(Alvm(q).determinant());
            
            // Slv.setZero();
            // Alvm.setZero();
            // uq = ucopy.col(q);
            
            // Define covariance matrix
            if(cstruc==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else {
              DiSc.setZero();
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
            matrix <Type> Alvmblock;
            matrix <Type> ucopyblock;
            for (i=0; i<nu; i++) {
              Alvmblock = Alvm(q).block(i*times,i*times,times,times)*Alvm(q).block(i*times,i*times,times,times).transpose();
              ucopyblock = ucopy.block(i*times,q,times,1);
              nll -=  0.5*(- (Slvinv*Alvmblock).trace()-(ucopyblock.transpose()*Slvinv*ucopyblock).sum());
              // nll -= log(Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix().determinant()) - 0.5*((atomic::matinv(Slv(q))*(Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix()*Alvm.col(q).matrix().block(i*times,i*times,times,times).matrix().transpose())).trace()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(atomic::matinv(Slv(q))*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // nll -= 0.5*(log((Alvm.col(i).matrix()*Alvm.col(i).matrix().transpose()).determinant()) - (Slv(q).inverse()*(Alvm.col(i).matrix()*Alvm.col(i).matrix().transpose())).trace()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(Slv(q).inverse()*(ucopy.block(i*times,q,times,1).matrix()))).sum());
              // log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
            }
            
          }
          
        } else if(num_corlv>1){
          // Kron A=AQ*Alvm
          matrix<Type> Alvm(times*nu,times*nu);
          Alvm.setZero();
          // Variational covariance
          //diagonal
          for(i=0; i<nu; i++){
            for (d=0; d<(times); d++){
              Alvm(i*times+d,i*times+d)=exp(Au(i*times+d));
            }
          }
          
          //reduced rank cov
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
            cQ.col(j) += 0.5*(Alvm*Alvm.transpose()).diagonal().matrix()*((newlamCor.col(j).transpose()*(AQ*AQ.transpose()))*newlamCor.col(j));
            // // cQ.col(j) += 0.5*Alvm.diagonal().matrix()*((newlamCor.col(j).transpose()*AQ)*newlamCor.col(j));
          }
          nll -= num_corlv*log(Alvm.determinant()) + times*nu*log(AQ.determinant()) + 0.5*num_corlv*times*nu;
          // nll -= 0.5*num_corlv*log(Alvm.determinant()) + 0.5*times*nu*log(AQ.determinant()) + 0.5*num_corlv*times*nu;
          //
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            // Slv.setZero();
            
            // Define covariance matrix
            if(cstruc==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else {
              DiSc.setZero();
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
            matrix <Type> Alvmblock;
            matrix <Type> ucopyblock;
            for (i=0; i<nu; i++) {
              Alvmblock = Alvm.col(q).matrix().block(i*times,i*times,times,times)*Alvm.col(q).matrix().block(i*times,i*times,times,times).transpose();
              ucopyblock = ucopy.block(i*times,q,times,1);
              nll -=  0.5*(- (AQ.row(q)*AQ.row(q).transpose()).sum()*(Slvinv*Alvmblock).trace() - (ucopyblock.transpose()*(Slvinv*ucopyblock)).sum());
              // nll -=  0.5*(- AQ(q,q)*(atomic::matinv(Slv(q))*Alvm.block(i*times,i*times,times,times).matrix()).trace()-((ucopy.block(i*times,q,times,1).matrix()).transpose()*(atomic::matinv(Slv(q))*(ucopy.block(i*times,q,times,1).matrix()))).sum());
            }
            
          }
        }
        
      }
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
              eta(i,j) -=  x_lv.row(i)*b_lv2*D(j)*(x_lv.row(i)*b_lv2).transpose();
            }
          }
        }
      }
      lam = u*newlam;
      
      // Update cQ for non quadratic latent variable model and
      // also takes this route if there are quadratic constrained LVs with random row-effect
      if((quadratic < 1) || ((nlvr==1 && random(2)<0 && num_RR>0))){
        
        //Binomial, Gaussian, Ordinal
        for (int i=0; i<n; i++) {
          for (int j=0; j<p;j++){
            cQ(i,j) += 0.5*(newlam.col(j).transpose()*A(i)*A(i).transpose()*newlam.col(j)).sum();
          }
        }
        eta += lam;
      }
      // do not take this route not with quadratic model, constrained LVs and random row-effects.
      if(( (quadratic>0 && nlvr > 0) && (num_lv+num_lv_c)>0 ) || ( quadratic>0 && random(2) > 0 )){
        //quadratic model approximation
        //Poisson
        e_eta = matrix <Type> (n,p);
        if(family==0){
          matrix<Type> Binv(nlvr,nlvr);
          matrix<Type> Cinv(nlvr,nlvr);
          matrix<Type> Acov(nlvr,nlvr);
          vector<Type> AcholDiag(nlvr);
          Type detB;
          Type detA;
          matrix <Type> vBinvv(n,p);
          matrix <Type> Id(nlvr,nlvr);
          Id.setZero();
          Id.diagonal().fill(1.0);
          // matrix <Type> Id = Matrix<Type, Dynamic, Dynamic>::Identity(nlvr, nlvr);
          for (int i=0; i<n; i++) {
            Acov = A(i)*A(i).transpose();
            // Qliu =  A(i).solve(u.row(i).transpose());
            for (int j=0; j<p;j++){
              //does not follow calculation from van der Veen et al. 2021
              //but prevents Acov^-1 via woodbury matrix identity
              Cinv = (Id + 2*Acov*D(j)).inverse();
              Binv = Acov-2*Cinv*Acov*D(j)*Acov;
              //this calculation prevents having to explicitly invert A*A^t, or having to invert A(i).
              vBinvv(i,j) = (newlam.col(j).transpose()*Acov*newlam.col(j)).sum()+2*(newlam.col(j).transpose()*u.row(i).transpose()).sum()-2*(newlam.col(j).transpose()*Cinv*Acov*D(j)*Acov*newlam.col(j)).sum()-
              4*(newlam.col(j).transpose()*Cinv*Acov*D(j)*u.row(i).transpose()).sum()-2*(u.row(i)*D(j)*Cinv*u.row(i).transpose()).sum();
              
              if((random(2)<1) && (num_lv_c>0)){
              //last term is extra for concurrent ordination
               vBinvv(i,j) += -4*(newlam.col(j).transpose()*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum()-4*(u.row(i)*(Id-2*Cinv*D(j)*Acov)*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum()+4*(x_lv.row(i)*b_lv2*D(j)*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum();
              }
              
              AcholDiag = A(i).diagonal();
              detB = atomic::logdet(Binv);//log-determinant of inverse, use as -detB
              detA = (AcholDiag.log()).sum(); //log-determinant of cholesky
              e_eta(i,j) = exp(cQ(i,j) + eta(i,j) + 0.5*(vBinvv(i,j)+detB)-detA); //add all the other stuff to the quadratic approximation
              eta(i,j) += lam(i,j) - (u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*A(i)*A(i).transpose()).trace();
              
              if((random(2)<1) && (num_lv_c>0)){
                eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
        //NB, gamma, exponential
        if((family==1)||(family==4)||(family==8)){
          matrix<Type> Binv(nlvr,nlvr);
          matrix<Type> Cinv(nlvr,nlvr);
          matrix<Type> Acov(nlvr,nlvr);
          vector<Type> AcholDiag(nlvr);
          Type detB;
          Type detA;
          matrix <Type> vBinvv(n,p);
          matrix <Type> Id(nlvr,nlvr);
          Id.setZero();Id.diagonal().fill(1.0);
          for (int i=0; i<n; i++) {
            Acov = A(i)*A(i).transpose();
            for (int j=0; j<p;j++){
              //does not follow calculation from van der Veen et al. 2021
              //but prevents Acov^-1 via woodbury matrix identity
              Cinv = (Id - 2*Acov*D(j)).inverse();
              Binv = Acov+2*Cinv*Acov*D(j)*Acov;
              //this calculation prevents having to explicitly invert A*A^t, or having to invert A(i).
              vBinvv(i,j) = (newlam.col(j).transpose()*Acov*newlam.col(j)).sum()-2*(newlam.col(j).transpose()*u.row(i).transpose()).sum()+2*(newlam.col(j).transpose()*Cinv*Acov*D(j)*Acov*newlam.col(j)).sum()-
              4*(newlam.col(j).transpose()*Cinv*Acov*D(j)*u.row(i).transpose()).sum()+2*(u.row(i)*D(j)*Cinv*u.row(i).transpose()).sum();

              if((random(2)<1) && (num_lv_c>0)){
                //last term is extra for concurrent ordination
                vBinvv(i,j) += -4*(newlam.col(j).transpose()*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum()+4*(u.row(i)*(Id+2*Cinv*D(j)*Acov)*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum()+4*(x_lv.row(i)*b_lv2*D(j)*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).sum();
              }
              
              AcholDiag = A(i).diagonal();
              detB = log((Binv.llt().matrixL()).determinant());
              detA = (AcholDiag.log()).sum(); //log-determinant of cholesky
              e_eta(i,j) = exp(-eta(i,j) - cQ(i,j) + 0.5*(vBinvv(i,j))-detA+detB);
              eta(i,j) += lam(i,j) - (u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*A(i)*A(i).transpose()).trace();
              
              if((random(2)<1) && (num_lv_c>0)){
                eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
        
        // Binomial, Gaussian, Ordinal
        if((family==2)||(family==3)||(family==7)){
          for (int i=0; i<n; i++) {
            for (int j=0; j<p;j++){
              if((random(2)>0) || ((num_lv>0)&(num_lv_c==0))){
                cQ(i,j) += 0.5*(newlam.col(j)*newlam.col(j).transpose()*A(i)*A(i).transpose()).trace() + (D(j)*A(i)*A(i).transpose()*D(j)*A(i)*A(i).transpose()).trace() +2*(u.row(i)*D(j)*A(i)*A(i).transpose()*D(j)*u.row(i).transpose()).sum() - 2*(u.row(i)*D(j)*A(i)*A(i).transpose()*newlam.col(j)).sum();
              }else if(random(2)<1){
                //extra terms for concurrent ordination
                cQ(i,j) += 0.5*((newlam.col(j)-2*D(j)*(x_lv.row(i)*b_lv2).transpose())*(newlam.col(j)-2*D(j)*(x_lv.row(i)*b_lv2).transpose()).transpose()*A(i)*A(i).transpose()).trace() + (D(j)*A(i)*A(i).transpose()*D(j)*A(i)*A(i).transpose()).trace() +2*(u.row(i)*D(j)*A(i)*A(i).transpose()*D(j)*u.row(i).transpose()).sum() - 2*(u.row(i)*D(j)*A(i)*A(i).transpose()*(newlam.col(j)-2*D(j)*(x_lv.row(i)*b_lv2).transpose())).sum();
              }
              eta(i,j) += lam(i,j) - (u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*A(i)*A(i).transpose()).trace();
              if((num_lv_c>0) && (random(2)<1)){
                eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
      }
    }
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
    } else if((family == 1) && (method<1)){//NB VA
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
      
    } else if ((family == 1) && (method>1)) { // NB EVA
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
    } else if((family == 2) && (method<1)) {//binomial probit VA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
          nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j)))) - cQ(i,j);
        }
      }
    } else if ((family == 2) && (method>1)) { // Binomial EVA
      if (extra(0) == 0) { // logit
        Type mu_prime;
        CppAD::vector<Type> z(4);
        
        for (int i=0; i<n; i++) {
          for (int j=0; j<p; j++) {
            // nll -= gllvm::dbinom_logit_eva(y(i,j), eta(i,j), cQ(i,j));
            
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            
            z[0] = eta(i,j);
            z[1] = 0;
            z[2] = 1/(1+exp(-z[0]));
            z[3] = exp(z[0])/(exp(z[0])+1);
            
            mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
            mu_prime = mu(i,j) * (1-mu(i,j));
            nll -= y(i,j) * eta(i,j) + log(1-mu(i,j));
            nll += mu_prime*cQ(i,j);
          }
        }
      } else if (extra(0) == 1) { // probit
        Type etaP;
        for (int i=0; i<n; i++) {
          for (int j=0; j<p; j++) {
            etaP = pnorm_approx(Type(eta(i,j)));   //pnorm funktion approksimaatio
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
      zetanew.setZero();
      
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
      zetanew.setZero();
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
      Type mu_prime;
      Type mu_prime2;
      CppAD::vector<Type> z;
      if(extra(0)==0){
        z = CppAD::vector<Type> (4);
      }
      CppAD::vector<Type> a(2);
      CppAD::vector<Type> b(2);
      CppAD::vector<Type> aa;
      CppAD::vector<Type> bb;
      Type dig_a;
      Type dig_b;
      Type trig_a;
      Type trig_b;
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          // define mu, mu' and mu''
          mu(i,j) = 0.0;
          mu_prime = 0.0;
          mu_prime2 = 0.0;
          if (extra(0) == 0) { // logit
            
            z[0] = eta(i,j);
            z[1] = 0;
            z[2] = 1/(1+exp(-z[0]));
            z[3] = exp(z[0])/(exp(z[0])+1);
            
            mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
            mu_prime = mu(i,j) * (1-mu(i,j));
            mu_prime2 = mu_prime * (1-2*mu(i,j));
            
          } else if (extra(0) == 1) { // probit
            mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
            mu_prime = dnorm(eta(i,j), Type(0), Type(1));
            mu_prime2 = (-eta(i,j))*mu_prime;
          }
          a[0] = mu(i,j)*iphi(j);
          a[1] = 1;
          b[0] = (1-mu(i,j))*iphi(j);
          b[1] = 1;
          aa = a;
          bb = b;
          aa[1] = 2;
          bb[1] = 2;
          dig_a = Type(atomic::D_lgamma(a)[0]);
          dig_b = Type(atomic::D_lgamma(b)[0]);
          trig_a = Type(atomic::D_lgamma(aa)[0]);
          trig_b = Type(atomic::D_lgamma(bb)[0]);
          
          nll -= dbeta(squeeze(y(i,j)), Type(a[0]), Type(b[0]), 1);
          nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
          nll -= iphi(j) * mu_prime2 * (log(squeeze(y(i,j))) - log(1-squeeze(y(i,j)))) * cQ(i,j);
          
        }
      }
    }
    // nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
    
  return nll;
}