#include <TMB.hpp>
#include "distrib.h"
#define TMB_LIB_INIT R_init_gllvmLA

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
  
  // PARAMETER_VECTOR(scaledc);// scale parameters for dc, of length of dc.cols()
  PARAMETER_VECTOR(zeta); // ordinal family param
  
  DATA_VECTOR(extra); // extra values, power of 
  DATA_INTEGER(model);// which model, basic or 4th corner
  DATA_IVECTOR(random);//(0)1=random, (0)0=fixed row params, for Br: (1)1 = random slopes, (1)0 = fixed, for b_lv: (2)1 = random slopes, (2)0 = fixed slopes
  DATA_INTEGER(zetastruc); //zeta param structure for ordinal model
  DATA_INTEGER(rstruc); //Type for random rows. default = 0, when same as u:s. If 1, dr0 defines the structure. If 2, Points within groups has covariance struct defined by cstruc
  DATA_INTEGER(times); //number of time points
  DATA_INTEGER(cstruc); //correlation structure for row.params, 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm
  DATA_MATRIX(dc); //coordinates for sites, used for exponentially decaying cov. struc
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
  
  vector<Type> iphi = exp(lg_phi);
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
  vector<matrix<Type>> Sigmab_lv;
  
  if(random(2)>0){
    vector<matrix<Type>> Sigmab_lv(sbl3);
    sigmab_lv = exp(sigmab_lv);
    sigmab_lv *= sigmab_lv;
    
    if(sigmab_lv.size()==(num_lv_c+num_RR)){//randomB=="LV", Sigma_q = sigma_q I_klv
      Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Sigmab_lvtemp(sbl12);
      for (int q=0; q<(num_RR+num_lv_c); q++){
        Sigmab_lv(q) = Sigmab_lvtemp;
        Sigmab_lv(q).diagonal().array() = sigmab_lv(q);
      }
    }else if(sigmab_lv.size()==Klv){//randomB=="P", Sigma_klv = sigma_klv I_d
      Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Sigmab_lvtemp(sbl12);
      for (int klv=0; klv<Klv; klv++){
        Sigmab_lv(klv) = Sigmab_lvtemp;
        Sigmab_lv(klv).diagonal().array() = sigmab_lv(klv);
      }
    }else if(sigmab_lv.size()==Type(1)){
      Eigen::DiagonalMatrix<Type,Eigen::Dynamic> Sigmab_lvtemp(sbl12);
      for (int klv=0; klv<Klv; klv++){
        Sigmab_lv(klv) = Sigmab_lvtemp;
        Sigmab_lv(klv).diagonal().array() = sigmab_lv(0);
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
  
  using namespace density;

  if(random(2)>0){
      // REPORT(Sigmab_lv); //!!!!
      //MVNORM_t<Type> mvnorm(Sigmab_lv);
      if(sbl3==Klv){
        for (int klv=0; klv<Klv; klv++) {
          nll += MVNORM(Sigmab_lv(klv))(b_lv.row(klv));
        }
      }else{
        for (int q=0; q<(num_lv_c+num_RR); q++) {
          nll += MVNORM(Sigmab_lv(q))(b_lv.col(q));
        }
      }
      
    }
    
    // REPORT(ucopy);
    // REPORT(num_corlv);
    // REPORT(nu);
    // REPORT(dr);
    // REPORT(cstruc);
    
    matrix<Type> etaH(n,p); 
    etaH.setZero();
    
    //For fixed-effects RRR with and without quadratic term
    if(num_RR>0){
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      if(quadratic>0){
        Eigen::DiagonalMatrix<Type,Dynamic> D_RR(num_RR);
        D_RR.setZero();
        if(lambda2.cols()==1){
          for (int d=0; d<num_RR;d++){
            D_RR.diagonal()(d) = fabs(lambda2(d,0));
          }
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose();
            }
          }
          
        }else{
          for (int j=0; j<p;j++){
            for (int d=0; d<num_RR;d++){
              D_RR.diagonal()(d) = fabs(lambda2(d,j));
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
      vector<Type> Brcol;
      for (int j=0; j<p;j++){
        Brcol = Br.col(j);
        nll += VECSCALE(neg_log_MVN,sdsv)(Brcol);
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
    if(((random(0)>0) && (nlvr==(num_lv+num_lv_c))) && (rstruc>0)){
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
            DiSc.setZero();
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
          DiSc.setZero();
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
          Slv.setZero();
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
            Slv.setZero();
            
            if(cstruc==1){// AR1 covariance
              Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
            } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
              Slv = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
            } else {
              DiSc.setZero();
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
          Slv.setZero();
          // Define covariance matrix
          if(cstruc==1){// AR1 covariance
            Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
          } else if(cstruc==3) {// Compound Symm  if(cstruc==3)
            Slv = gllvm::corCS(Type(1), rho_lvc(q,0), times);
          } else {
            DiSc.setZero();
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
      if((num_RR>0) && (nlvr == 0) && random(2)>1){
        //use dnbinom_robust in this case - below code does not function well
        //for constrained ordination without any random-effects
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
          nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
          }
        }
      }else{
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            nll -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
          }
        } 
      }
      } else if(family==2) {//binomial family
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
      } else if((family==7) && (zetastruc == 1)){//ordinal, only here for models without random-effects
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
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
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
  
  return nll;
}