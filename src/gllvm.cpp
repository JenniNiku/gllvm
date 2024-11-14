#define TMB_LIB_INIT R_init_gllvm
#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#include<math.h>
#include "distrib.h"
#include "utils.h"

template<class Type>
struct dclist : vector<matrix<Type>> {
  dclist(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

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
  DATA_IMATRIX(csb_lv);
  DATA_MATRIX(xr); // design matrix for fixed row effects
  DATA_MATRIX(xb); // design matrix for random species effects, (n, spnr)
  DATA_SPARSE_MATRIX(dr0); // design matrix for rows, ( n, nr)
  DATA_SPARSE_MATRIX(dLV); // design matrix for latent variables, (n, nu)
  DATA_IMATRIX(cs); // matrix with indices for correlations of random species effects
  // DATA_SPARSE_MATRIX(colMat); //lower cholesky of column similarity matrix
  DATA_STRUCT(colMatBlocksI, dclist); //first entry is the number of species in colMat, rest are the blocks of colMat
  DATA_IMATRIX(nncolMat);
  DATA_VECTOR(Abranks);
  DATA_MATRIX(offset); //offset matrix
  DATA_IVECTOR(Ntrials);
  
  PARAMETER_MATRIX(r0f); // fixed site/row effects
  PARAMETER_MATRIX(r0r); // random site/row effects
  PARAMETER_MATRIX(b); // matrix of species specific intercepts and coefs
  // PARAMETER_MATRIX(bH); // matrix of species specific intercepts and coefs for beta hurdle model
  PARAMETER_MATRIX(B); // coefs of 4th corner model for TMBtrait, RE means for gllvm.TMB
  PARAMETER_MATRIX(Br); // random slopes for env species
  PARAMETER_MATRIX(b_lv); //slopes for RRR and constrained ord, VA means for random slopes
  //Left columns are for constrained ordination, Right for RRR
  PARAMETER_VECTOR(sigmaLV);//SD for LV
  PARAMETER_VECTOR(lambda); // lv loadings
  PARAMETER_MATRIX(lambda2);// quadratic lv loadings
  // PARAMETER_MATRIX(thetaH);// hurdle model lv loadings
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phiZINB);//extra param for ZINB
  PARAMETER_VECTOR(lg_phi); // dispersion params/extra zero probs for ZIP
  PARAMETER_VECTOR(sigmaB); // sds for random species effects
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
  DATA_INTEGER(randomB) //0 = P,single,iid and 1 = LV
  PARAMETER_VECTOR(Au); // variational covariances for u
  PARAMETER_VECTOR(lg_Ar); // variational covariances for r0r
  PARAMETER_VECTOR(Abb);  // variational covariances for Br
  // PARAMETER_VECTOR(scaledc);// scale parameters for dc, of length of dc.cols()
  PARAMETER_VECTOR(Ab_lv); //variational covariances for b_lv
  PARAMETER_VECTOR(zeta); // ordinal family param
  
  PARAMETER(ePower);
  DATA_VECTOR(extra); // extra values, power of 
  DATA_INTEGER(method);// 0=VA, 1=LA, 2=EVA
  DATA_INTEGER(Abstruc); //0 = diagonal, blockdiagonal, 1 = MNdiagonal, MNunstructured, 2 = diagonalCL2, 3 = diagonalCL1, CL1, 4 = CL2, 5 = unstructured
  DATA_INTEGER(model);// which model, basic or 4th corner
  DATA_IVECTOR(random);//(0)1=random, (0)0=fixed row params, for Br: (1)1 = random slopes, (1)0 = fixed, for b_lv: (2)1 = random slopes, (2)0 = fixed slopes, for Br: (3) 1 = random
  DATA_INTEGER(zetastruc); //zeta param structure for ordinal model
  DATA_IVECTOR(nr); // number of observations in each random row effect
  DATA_INTEGER(times); // number of time points, for LVs
  DATA_IVECTOR(cstruc); //correlation structure for row.params 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm, 4= Matern
  DATA_INTEGER(cstruclv); //correlation structure for LVs 0=indep sigma*I, 1=ar1, 2=exponentially decaying, 3=Compound Symm, 4= Matern
  DATA_STRUCT(dc, dclist); //coordinates for sites, used for exponentially decaying cov. struc
  DATA_MATRIX(dc_lv); //coordinates for sites, used for exponentially decaying cov. struc
  DATA_INTEGER(Astruc); //Structure of the variational covariance, 0=diagonal, 1=RR, (2=sparse cholesky not implemented yet)
  DATA_IMATRIX(NN); //nearest neighbours,
  
  int Klv = x_lv.cols();
  int n = y.rows();
  int p = y.cols();
  // int nt =n;
  int nu =n; //CorLV
  
  if(num_corlv>0){ //CorLV
    nu = dLV.cols();
  }
  
  vector<Type> iphi = exp(lg_phi);
  
  // Set first row param to zero, if row effects are fixed
  // if(random(0)==0 && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){  r0f(0,0) = 0;}
  int nlvr = num_lv+num_lv_c;//treating constr. ord random slopes as a LV, to use existing infrastructure for integration
  
  matrix<Type> ucopy = u;
  if(num_corlv>0){
    nlvr=0; num_lv=0; num_lv_c=0;
    quadratic=0;
    num_RR=0;
  }
  
  // Distance matrix calculated from the coordinates for LVs
  matrix<Type> DiSc_lv(dc_lv.cols(),dc_lv.cols()); DiSc_lv.fill(0.0);
  matrix<Type> dc_scaled_lv(dc_lv.rows(),dc_lv.cols()); dc_scaled_lv.fill(0.0);
  // matrix<Type> DistM(dc.rows(),dc.rows());
  // if(((num_corlv>0) || (((random(0)>0) & (nlvr==(num_lv+num_lv_c))) & (rstruc>0))) & ((cstruc(0)==2) || (cstruc(0)>3))){
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
  
  matrix<Type> eta(n,p);
  eta.setZero();
  matrix<Type> lam(n,p);
  lam.setZero();
  // Type nll = 0;
  parallel_accumulator<Type> nll(this); // initial value of log-likelihood
  
  matrix<Type> RRgamma(num_RR,p);
  RRgamma.setZero();
  
  matrix <Type> Delta(nlvr,nlvr);
  Delta.setZero();
  
  matrix <Type> newlam(nlvr,p);
  newlam.setZero();  
  
  //K*K*d or d*d*K
  int sbl12 = Klv;
  int sbl3 = num_lv_c + num_RR;

  if(randomB>0){
    sbl12 = num_lv_c + num_RR;;
    sbl3 = Klv;
  }
  
  vector<matrix<Type>> Sigmab_lv(sbl3);
  
  if(random(2)>0){
    for (int q=0; q<sbl3; q++){
      Sigmab_lv(q).resize(sbl12,sbl12);
      Sigmab_lv(q).setIdentity();
    }

    if(csb_lv.cols()<2){
    if(randomB<1){//Sigma_q = sigma_q I_klv
      //randomB="P","single","iid"
      if(sigmab_lv.size()>1){
      for (int q=0; q<(num_lv_c+num_RR); q++){
      Sigmab_lv(q).diagonal().array() = exp(sigmab_lv)*exp(sigmab_lv);
      }
      }else{
        for (int q=0; q<(num_lv_c+num_RR); q++){
          Sigmab_lv(q).diagonal().array() = exp(sigmab_lv(0))*exp(sigmab_lv(0));
        }
      }
      // matrix<Type> Sigmab_lvtemp(sbl12,sbl12);
      // Sigmab_lvtemp.setZero();
      // for (int q=0; q<sbl3; q++){
      //   Sigmab_lv(q) = Sigmab_lvtemp;
      //   Sigmab_lv(q).diagonal().array() = sigmab_lv(q);
      // }
    }else if(randomB>0){
      //randomB="LV"
        Sigmab_lv(0).diagonal().array() *= exp(sigmab_lv)*exp(sigmab_lv);
    }
    }else{
      matrix<Type> sds = Eigen::MatrixXd::Zero(x_lv.cols(),x_lv.cols());
      sds.diagonal() = exp(sigmab_lv.segment(0,x_lv.cols()));
      vector<Type>corsb_lv((x_lv.cols()*x_lv.cols()-x_lv.cols())/2);
      corsb_lv.fill(0.0);
      matrix<Type>Sigmab_lvL(x_lv.cols(),x_lv.cols());
      Sigmab_lvL.setIdentity();
      if(csb_lv.cols()>1){
        //need a vector with covariances and zeros in the right places
        for(int i=0; i<csb_lv.rows(); i++){
          corsb_lv((csb_lv(i,0) - 1) * (csb_lv(i,0) - 2) / 2 + csb_lv(i,1)-1) = sigmab_lv(x_lv.cols()+i);
        }
        Sigmab_lvL = sds*gllvmutils::constructL(corsb_lv);
    }
      for (int q=0; q<(num_lv_c+num_RR); q++){
      Sigmab_lv(q) = Sigmab_lvL;
      }
  }
  }
  
  if((nlvr>0)||(num_RR>0)){
    
    if(nlvr>0){
      newlam.row(0).fill(1.0);
      if((num_lv+num_lv_c)>0){
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
  if((method<1) || (method>1)){
    // add offset
    eta += offset;
    // add fixed row effects
    // if(r0f.size() == n && (random(0)==0) && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){
    //   eta += r0f.replicate(1,p);
    if(xr.rows()==n){
      eta += (xr*r0f).replicate(1,p);
    }
    
    matrix<Type> cQ(n,p);
    cQ.setZero();
    
    vector<matrix<Type>> A(n);
    
    if( (random(2)>0) && (num_RR>0) && (quadratic>0)){
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
      
      // Add VA terms to logL
      //Go this route if no random Bs
      matrix <Type> Atemp(nlvr,nlvr);
      for(int i=0; i<n; i++){
        Atemp = A(i).topLeftCorner(nlvr,nlvr);//to exlcude the 0 rows & columns for num_RR
        nll -= Atemp.diagonal().array().log().sum() - 0.5*((Atemp*Atemp.transpose()).trace()+(u.row(i)*u.row(i).transpose()).sum());
      }
      nll -= 0.5*n*nlvr;
      
      //scale LVs with standard deviations, as well as the VA covariance matrices
      u *= Delta;
      
      if((num_RR*random(2))>0 && (quadratic)>0){
        Delta.conservativeResize(nlvr+num_RR,nlvr+num_RR);
        for(int d=nlvr; d<(nlvr+num_RR); d++){
          Delta.col(d).setZero();
          Delta.row(d).setZero();
        }
      }
      
      for (int i=0; i<n; i++) {
        A(i) = Delta*A(i);
      }
      
    }
    
    //random slopes for constr. ord.
    vector<matrix<Type>> Ab_lvcov;  //covariance of LVs due to random slopes
    if((random(2)>0) && ((num_RR+num_lv_c)>0)){
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
            }
            k++;
          }}
      }

      //VA likelihood parts for random slope
      if(csb_lv.cols()<2){
      for(int q=0; q<sbl3; q++){
        if(randomB<1){
          nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(Sigmab_lv(q).diagonal().cwiseInverse().asDiagonal()*AB_lv(q)*AB_lv(q).transpose()).trace()-0.5*(b_lv.col(q).transpose()*Sigmab_lv(q).diagonal().cwiseInverse().asDiagonal()*b_lv.col(q)).sum());
          nll -= 0.5*(sbl12- Sigmab_lv(q).diagonal().array().log().sum());
        }
        if(randomB>0)nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(Sigmab_lv(0).diagonal().cwiseInverse().asDiagonal()*AB_lv(q)*AB_lv(q).transpose()).trace()-0.5*(b_lv.row(q)*Sigmab_lv(0).diagonal().cwiseInverse().asDiagonal()*b_lv.row(q).transpose()).sum());// log(det(A_bj))-sum(trace(S^(-1)A_bj))*0.5 + a_bj*(S^(-1))*a_bj
       }
      if(randomB>0)nll -= sbl3*0.5*(sbl12- Sigmab_lv(0).diagonal().array().log().sum());
        }else{
          matrix<Type>Iblv = Eigen::MatrixXd::Identity(x_lv.cols(),x_lv.cols());
          matrix<Type>Sigmab_lvI = (Sigmab_lv(0)*Sigmab_lv(0).transpose()).ldlt().solve(Iblv);
          for(int q=0; q<(num_lv_c+num_RR); q++){
          nll -= (AB_lv(q).diagonal().array().log().sum() - 0.5*(Sigmab_lvI*AB_lv(q)*AB_lv(q).transpose()).trace()-0.5*(b_lv.col(q).transpose()*Sigmab_lvI*b_lv.col(q)).sum());
          nll -= 0.5*Klv- Sigmab_lv(0).diagonal().array().log().sum();
        }
        
      }
      
      //resize ab_lvcov to correct size
      Ab_lvcov  = vector<matrix<Type>> (n);
      for(int i=0; i<n; i++){
        Ab_lvcov(i).resize(num_RR+num_lv_c,num_RR+num_lv_c);
        Ab_lvcov(i).setZero();
      }
      
      //fill ab_lvcov
      if(randomB>0){//variance per LV
        for(int klv=0; klv<Klv; klv++){
        matrix<Type>Ablv = AB_lv(klv)*AB_lv(klv).transpose();
        for(int i=0; i<n; i++){
            Ab_lvcov(i) += x_lv(i,klv)*x_lv(i,klv)*Ablv;
          }
        }
      }else{
      for(int q=0; q<(num_RR+num_lv_c); q++){
        matrix<Type>Ablv = AB_lv(q)*AB_lv(q).transpose();
        for(int i=0; i<n; i++){
            Ab_lvcov(i)(q,q) = (x_lv.row(i)*Ablv*x_lv.row(i).transpose()).sum();
          }
        }
      }
      
      if(num_lv_c>0){
        RRgamma.conservativeResize(num_RR+num_lv_c,Eigen::NoChange);
        if(num_RR>0)RRgamma.bottomRows(num_RR) = RRgamma.topRows(num_RR); 
        RRgamma.topRows(num_lv_c) = newlam.topRows(num_lv_c);
      }
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          cQ(i,j) += 0.5*(RRgamma.col(j).transpose()*Ab_lvcov(i)*RRgamma.col(j)).value();
        }
      }
      if(quadratic<1){
        eta += x_lv*b_lv*RRgamma;//for the quadratic model this component is added below
        
      }else if(quadratic>0){
        //now rebuild A and u with covariances for random slopes so that existing infrastructure below can be used
        //in essence, q(XBsigmab_lv + eDelta) ~ N(uDelta + \sum \limits^K X_ik b_lv_k , Delta A Delta + \sum \limits^K X_ik^2 AB_lv_k )
        //so build u and A accordingly (and note covariance due to Bs if num_lv_c and num_RR > 0)
        if(num_lv_c>0){
          u.leftCols(num_lv_c) += x_lv*b_lv.leftCols(num_lv_c);
        }
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
          newlam.bottomRows(num_RR) = RRgamma.bottomRows(num_RR);
          u.rightCols(num_RR) += x_lv*b_lv.rightCols(num_RR);
        }
        
        if((nlvr-num_RR-num_lv_c)>0){
          //rebuild Ab_lvcov to fit A below
          matrix<Type> tempRR(num_RR,num_RR);
          matrix<Type> tempCN(num_lv_c,num_lv_c);
          matrix<Type> tempRRCN(num_RR,num_lv_c);
          
          for(int i=0; i<n; i++){
            if(num_RR>0){
              tempRR = Ab_lvcov(i).bottomRightCorner(num_RR,num_RR);
            }
            if(num_lv_c>0){
              tempCN = Ab_lvcov(i).topLeftCorner(num_lv_c,num_lv_c);
            }
            if((num_RR+num_lv_c)>0){
              tempRRCN = Ab_lvcov(i).bottomLeftCorner(num_RR,num_lv_c);
            }
            //resize to fit A
            Ab_lvcov(i).conservativeResize(nlvr,nlvr);
            Ab_lvcov(i).setZero();
            
            //re-assign
            //place num_RR in back
            if(num_RR>0){
              Ab_lvcov(i).bottomRightCorner(num_RR,num_RR) = tempRR;
            }
            //num_lv_c is in front, but after a potential intercept
            if((num_lv_c)>0){
              Ab_lvcov(i).topLeftCorner(num_lv_c,num_lv_c) = tempCN;
            }
            //assign covariances of random slopes. There is no covariance if slb3==num_lv_c+num_RR, randomB == 0
            if((num_RR>0)&&(num_lv_c>0)&&(randomB>0)){
              Ab_lvcov(i).block(nlvr-num_RR,nlvr-num_lv_c-num_RR-num_lv,num_RR,num_lv_c) = tempRRCN;
              Ab_lvcov(i).block(nlvr-num_lv_c-num_RR-num_lv,nlvr-num_RR,num_lv_c,num_RR) = tempRRCN.transpose();
            }
            
          }
        }
      }
      
    }
    
    //components for reduced rank regression terms
    if((num_RR>0) && (random(2)<1)){
      //predictor coefficients RRR.  num_RR comes after num_lv_c
      //Since later dimensions are more likely to have less residual variance
      eta += x_lv*b_lv.rightCols(num_RR)*RRgamma;
      
      //quadratic terms for fixed-effects only RRR
      //-num_lv to ensure that we pick num_RR from the middle
      if(quadratic>0){
        matrix<Type> D_RR(num_RR,num_RR);
        D_RR.setZero();
        
        //quadratic coefficients for RRR
        if(lambda2.cols()==1){
          for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
            D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,0));
          }
          // for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            eta.row(i).array() -=  (x_lv.row(i)*b_lv.rightCols(num_RR)*D_RR*(x_lv.row(i)*b_lv.rightCols(num_RR)).transpose()).value();
          }
          // }
        }else{
          for (int j=0; j<p;j++){
            D_RR.setZero();
            for (int d=num_lv_c; d<(num_lv_c+num_RR);d++){
              D_RR.diagonal()(d-num_lv_c) = fabs(lambda2(d,j));
            }
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv.rightCols(num_RR)*D_RR*(x_lv.row(i)*b_lv.rightCols(num_RR)).transpose();
            }
            
          }
        }
        
      }
    }
    
    if((random(1)>0) | (random(3)>0)){
      if(random(1)>0){
        // random slopes in TMBtrait.R
        eta += xb*Br;  
      }else if(random(3)>0){
        //random slopes in gllvm.TMB
        if(B.rows()==xb.cols()){
          //RE means or community-level effects
          eta += (xb*B).replicate(1,p);
        }
        eta += xb*Br;
      }
      matrix <Type> Spr(xb.cols(),xb.cols());
      matrix <Type> SprI(xb.cols(),xb.cols());
      matrix <Type> SprIL(xb.cols(),xb.cols());
      
      Type logdetSpr = 0;
      
      int l = xb.cols();
      if(random(1)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB);
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL(l,l);
        SprL.fill(0.0);
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaij(i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      if(random(3)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB.segment(0,xb.cols()));
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL;
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaB(xb.cols()+i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      
      // vector<matrix<Type>> colCorMatIblocks(colMatBlocksI.size()-1);
      // Type logdetColCorMat =0;
      vector <Type> rhoSP(1);
      if(random(1)>0 && sigmaB.size()>xb.cols()){
        rhoSP.resize(sigmaB.size()-xb.cols());
      }else if(random(3)>0){
        rhoSP.resize(sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1));
      }
      
      if(colMatBlocksI(0)(0,0) == p){
        rhoSP.fill(1.0);
        if(sigmaB.size()>(xb.cols()+cs.rows()*(cs.cols()>1))){
          rhoSP = exp(-exp(sigmaB.segment(xb.cols()+cs.rows()*(cs.cols()>1),sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1))));
          if(nncolMat.rows()<p){
            //need to cap this on the lower end for numerical stability
            for(int re=0; re<rhoSP.size(); re++){
              rhoSP(re) = CppAD::CondExpLt(rhoSP(re), Type(1e-12), Type(1e-12), rhoSP(re));
            }
          }
        }
      }
      
      //Variational components
      int ncov = xb.cols();
      
      if(Abstruc == 0){
        //((Abb.size() ==(p*ncov)) || (Abb.size() == (p*ncov+p*ncov*(ncov-1)/2))) && Abranks(0) == 0
        // Ab.struct == "diagonal" or "blockdiagonal", i.e., independence over species
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        // vector<Type> SprIAtrace(p);
        vector<matrix <Type>> SArm(p);
        
        for(int j=0; j<p; j++){//limit the scope of SArmLst(j) to this iteration for memory efficiency
          SArm(j).resize(ncov,ncov);
          SArm(j).setZero();
          for (int d=0; d<ncov; d++){ // diagonals of varcov
            SArm(j)(d,d)=exp(Abb(sdcounter));
            sdcounter++;
          }
          
          if((Abb.size()>(p*ncov))){ // unstructured block Var.cov
            for (int d=0; d<(ncov); d++){
              for (int r=d+1; r<(ncov); r++){
                SArm(j)(r,d)=Abb(covscounter);
                covscounter++;
              }}
          }
          
          nll -= SArm(j).diagonal().array().log().sum();
          cQ.col(j) += 0.5*((xb*SArm(j)*SArm(j).transpose()).cwiseProduct(xb)).rowwise().sum();
        }
        // tr(S⁻¹A)+Br'S⁻¹Br+logdet(S)
        if((colMatBlocksI(0)(0,0)==p) && (rhoSP.size()==1)){
          int sp = 0;
          if(nncolMat.rows()<p){
            for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
              matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
              colCorMatI.setZero();
              gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
              nll -= -0.5*(Br.middleCols(sp, colMatBlocksI(cb+1).cols())*colCorMatI*Br.middleCols(sp, colMatBlocksI(cb+1).cols()).transpose()*SprI).trace();
              
              for(int j=0; j<(colMatBlocksI(cb+1).cols()); j++){
                nll -= -0.5*colCorMatI(j,j)*(SprI*SArm(sp+j)*SArm(sp+j).transpose()).trace();  
              }
              sp += colMatBlocksI(cb+1).cols();
            }
          }else{
            for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
              Eigen::SparseMatrix<Type> colCorMatUI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
              gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb+1).cols()));
              Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
              nll -= -0.5*(Br.middleCols(sp, colMatBlocksI(cb+1).cols())*colCorMatI*Br.middleCols(sp, colMatBlocksI(cb+1).cols()).transpose()*SprI).trace();
              vector<Type>colCorMatIdiag = colCorMatI.diagonal();
              for(int j=0; j<(colMatBlocksI(cb+1).cols()); j++){
                nll -= -0.5*colCorMatIdiag(j)*(SprI*SArm(sp+j)*SArm(sp+j).transpose()).trace();  
              }
              sp += colMatBlocksI(cb+1).cols();
            }
          }
          
          //determinant
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
          
        }else if((colMatBlocksI(0)(0,0)==p) && (rhoSP.size()>1)){
          // calculate trace term without building the covariance matrix
          // requires a little bit of work with a temporary matrix..
          int sp = 0;
          for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
            int blocksize = colMatBlocksI(cb+1).cols();
            vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUI(ncov);
            
            for (int d=0; d<ncov; d++){
              colCorMatUI(d).resize(blocksize,blocksize);
              gllvmutils::nngp(colCorMatUI(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              Br.row(d).middleCols(sp, blocksize) *= colCorMatUI(d);
              colCorMatUI(d) = colCorMatUI(d).transpose(); //now lower triangular, needed to expose .col below
            }
            nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            
            matrix<Type>tempMat(ncov,ncov);
            for(int j=0; j<blocksize; j++){
              for (int d=0; d<ncov; d++){
                for (int d2=d+1; d2<ncov; d2++){
                  tempMat(d2,d) = (colCorMatUI(d).col(j).transpose()*colCorMatUI(d2).col(j)).sum()*SprI(d,d2);
                  tempMat(d,d2) = tempMat(d2,d);
                }
                tempMat(d,d) = (colCorMatUI(d).col(j).transpose()*colCorMatUI(d).col(j)).sum()*SprI(d,d);
              }
              nll -= -0.5*(tempMat*SArm(sp+j)*SArm(sp+j).transpose()).trace();
            }
            sp += colMatBlocksI(cb+1).cols();
          }
          
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }else{
          for(int j=0; j<p; j++){
            nll -= -0.5*(SprI*SArm(j)*SArm(j).transpose()).trace();
          }
          nll -= 0.5*(p*ncov-p*logdetSpr); 
          nll -= -0.5*(Br*Br.transpose()*SprI).trace();
        }
      }else if(Abstruc == 1){
        //(Abb.size()==(ncov+p-1+p*(p-1)/2)) || (Abb.size() == (p+ncov-1 + p*(p-1)/2 + ncov*(ncov-1)/2)) || (Abb.size() == (p+ncov-1+ncov*(ncov-1)/2+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum())) || (Abb.size() == (p+ncov-1+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "MNdiagonal" OR "MNunstructured" with blocks due to Phylogeny
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        matrix<Type>SArmR(ncov, ncov);//row covariance matrix
        SArmR.setZero();
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p+ncov-1;
        
        //row covariance matrix
        for (int d=0; d<(ncov); d++){
          SArmR.diagonal()(d)=exp(Abb(sdcounter));
          sdcounter++;
        }
        
        if((Abb.size()>(ncov+p-1+p*(p-1)/2)) || (Abb.size()>(p+ncov-1+(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))){ // unstructured row covariance
          for (int d=0; d<(ncov); d++){
            for (int r=d+1; r<(ncov); r++){
              SArmR(r,d)=Abb(covscounter);
              covscounter++;
            }}
        }
        
        // determinant of kronecker product of two matrices based on their cholesky factors: first part
        nll -= p*SArmR.diagonal().array().log().sum();
        SArmR *= SArmR.transpose();
        
        //for later use
        matrix<Type> xbSArmxb = ((xb*SArmR).cwiseProduct(xb)).rowwise().sum();
        
        int sp = 0;//tracking how many species we have had so far
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          if(Abranks(cb)<blocksize){
            //use sparse matrices if we can to speed things up with many species
            // Eigen::SparseMatrix<Type>SArmC(blocksize,blocksize);
            // std::vector<T> tripletList;
            matrix<Type> SArmC(blocksize, CppAD::Integer(Abranks(cb)));
            SArmC.setZero();
            Eigen::Vector<Type, Eigen::Dynamic>SArmCD(blocksize);//remaining variances
            SArmCD.setZero();
            // set-up column (block) covariance matrix
            int jstart = 0;
            if(cb==0){
              //first entry fixed for identifiability
              SArmC(0,0) = 0.3;
              jstart = 1;
            }
            
            for (int j=jstart; j<blocksize; j++){
              if(j<Abranks(cb)){
                SArmC(j,j) = exp(Abb(sdcounter));
                SArmCD(j) = 0;
              }else{
                SArmCD(j) = exp(Abb(sdcounter));
              }
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<blocksize; r++){
                SArmC(r,j) =  Abb(covscounter);
                covscounter++;
              }
            }
            // determinant of kronecker product of two matrices based on their cholesky factors: second part
            nll -= SArmR.rows()*(SArmC.diagonal().array().log().sum() + SArmCD.array().tail(blocksize-SArmC.cols()).log().sum());
            // matrix<Type>SArmP = SArmC*SArmC.transpose();
            // SArmP.diagonal() += SArmCD.array().pow(2).matrix();
            
            //tr(S⁻¹A) and Br'S⁻¹Br
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                //not sparse
                matrix<Type>colCorMatI(colMatBlocksI(cb+1).cols(),colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(((SArmC.transpose()*colCorMatI*SArmC).trace()+(SArmCD.array().pow(2)*colCorMatI.diagonal().array()).sum())*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else{
                //NN sparse approximation
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Type traceUUtSArmP = ((SArmC.transpose()*colCorMatUI)*(colCorMatUI.transpose()*SArmC)).trace();
                //remaining diagonal entries
                // equivalent to diag(SArmCDs^2)UU', but via the cholesky instead
                // diag(SArmCDs)U, or the squared norm of the rows (),
                // so that I do not need to compute the whole product
                for (int j=0; j<blocksize; j++){
                  vector<Type> temp = colCorMatUI.col(j).cwiseProduct(SArmCD);
                  traceUUtSArmP += temp.pow(2).sum();
                }
                // Type traceold = (colCorMatUI*colCorMatUI.transpose()*SArmP).trace();
                // REPORT(traceold);
                // REPORT(traceUUtSArmP);
                nll -= -0.5*(traceUUtSArmP*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
                
                // nll -= -0.5*((colCorMatUI*colCorMatUI.transpose()*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              //tr(S^-1A)
              vector<Eigen::SparseMatrix<Type>> colCorMatUIBlocks(ncov);
              for (int d=0; d<ncov;d++){
                colCorMatUIBlocks(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIBlocks(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                //Br'S⁻¹Br
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIBlocks(d);
              }
              //this is written in a way that keeps dimensions of the
              //objects small thanks to the reduced rank
              //which is much more friendly in memory
              //based on the formulation A_m = LL'+D, where D is diagonal given by SArmCD^2
              // trace: tr((LL'+D)UU')*Ar(d,d2)*SprI(d,d2), split over the addition for efficiency
              
              matrix<Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
              vector<Type> temp(blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmC.transpose()*colCorMatUIBlocks(d);
                for (int d2=d+1; d2<ncov;d2++){
                  // nll -= -(SArmP*colCorMatUIBlocks(d)*colCorMatUIBlocks(d2).transpose()).trace()*SArmR(d,d2)*SprI(d,d2);
                  //reduced rank part: L
                  nll -= -((tempMat*(colCorMatUIBlocks(d2).transpose()*SArmC)).trace())*SArmR(d,d2)*SprI(d,d2);
                  // apparently faster than "just" computing the diagonal entries via an elementwise product
                  // I suppose Eigen can do some optimization with this
                  vector<Type> diags(blocksize);
                  diags.setZero();
                  for (int j=0; j<blocksize; j++){
                    diags += colCorMatUIBlocks(d).col(j).cwiseProduct(colCorMatUIBlocks(d2).col(j));
                  }
                  nll -= -(SArmCD.array().pow(2)*diags).sum()*SArmR(d,d2)*SprI(d,d2);
                }
                //reduced rank part
                nll -= -0.5*(tempMat).rowwise().squaredNorm().sum()*SArmR(d,d)*SprI(d,d);
                Type trace = 0;
                //remaining diagonal entries
                //need to do this column-wise so I do not need to compute the whole product
                //similar to above, trace via squared norm
                for (int j=0; j<blocksize; j++){
                  temp = colCorMatUIBlocks(d).col(j).cwiseProduct(SArmCD);
                  trace += temp.pow(2).sum();
                }
                nll -= -0.5*trace*SArmR(d,d)*SprI(d,d);
              }
              //nll -= -0.5*(SprIL.transpose()*Br.middleCols(sp, blocksize)).rowwise().squaredNorm().sum();//is slower it seems
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
            //retrieve diagonal outside loop for more efficient memory handling
            //calling it inside the loop constructs a temporary for each i
            matrix<Type>SArmPdiag = SArmC.rowwise().squaredNorm().array()+SArmCD.array().pow(2);
            cQ.middleCols(sp, blocksize).noalias() += 0.5*xbSArmxb*SArmPdiag.transpose();
            // for (int i=0; i<n;i++){
            //   cQ.row(i).middleCols(sp, blocksize) += 0.5*SArmPdiag*xbSArmxb(i);
            // } 
            sp += blocksize;
          }else{
            matrix<Type>SArmC;
            SArmC.resize(blocksize,blocksize);//column VA covariance
            SArmC.setZero();
            
            int jstart = 0;
            // set-up column (block) covariance matrix
            if(cb==0){
              SArmC(0,0) = 0.3;//first entry fixed for identifiability
              jstart = 1;
            }
            
            for (int j=jstart; j<SArmC.cols(); j++){
              SArmC(j,j) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                SArmC(r,j) = Abb(covscounter);
                covscounter++;
              }
            }
            
            // determinant of kronecker product of two matrices based on their cholesky factors: second part
            nll -= SArmR.rows()*SArmC.diagonal().array().log().sum();
            matrix<Type>SArmP = SArmC*SArmC.transpose();
            
            //tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                //not sparse
                matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*((colCorMatI*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else{
                //NN sparse approximation
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                nll -= -0.5*((colCorMatUI*colCorMatUI.transpose()*SArmP).trace()*(SprI*SArmR).trace()+(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              vector<Eigen::SparseMatrix<Type>> colCorMatUIBlocks(ncov);
              for (int d=0; d<ncov;d++){
                colCorMatUIBlocks(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIBlocks(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              }
              //tr(S⁻¹A)
              matrix<Type>tempMat(blocksize,blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmP*colCorMatUIBlocks(d);
                for (int d2=d+1; d2<ncov;d2++){
                  nll -= -(tempMat*colCorMatUIBlocks(d2).transpose()).diagonal().sum()*SArmR(d,d2)*SprI(d,d2);
                }
                nll -= -0.5*(tempMat*colCorMatUIBlocks(d).transpose()).diagonal().sum()*SArmR(d,d)*SprI(d,d);
              
              //Br'S⁻¹Br
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIBlocks(d);
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
            //retrieve diagonal outside loop for more efficient memory handling
            //calling it inside the loop constructs a temporary for each i
            matrix<Type> SArmPdiag = SArmP.diagonal();
            // for (int i=0; i<n;i++){
            //   cQ.row(i).middleCols(sp, blocksize) += 0.5*SArmPdiag*xbSArmxb(i);
            // } 
            cQ.middleCols(sp, blocksize).noalias() += 0.5*xbSArmxb*SArmPdiag.transpose();
            sp += blocksize;
          }
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 2){
        //(Abb.size() == ((ncov*p +ncov*p*(p-1)/2))) || (Abb.size()==(ncov*p+sum(nsp)*(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "diagonalCL2" with block structure due to phylogeny
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        
        int sp = 0;//keeping track of #species we've had due to blocks
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          matrix<Type>colCorMatI(blocksize,blocksize);
          vector<Eigen::SparseMatrix<Type>> colCorMatUIs(ncov);
          
          if((rhoSP.size()) == 1){
            gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
            if(nncolMat.rows()<p){
              nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }else if(nncolMat.rows()==p){
              colCorMatUIs(0).resize(blocksize, blocksize);
              gllvmutils::nngp(colCorMatUIs(0), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
              //sparse approximation to inverse
              nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatUIs(0)*colCorMatUIs(0).transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
            }
          }else{
            //NN sparse approximation
            //Br'S⁻¹Br'
            for (int d=0; d<(ncov); d++){
              colCorMatUIs(d).resize(blocksize, blocksize);
              gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
              Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
            }
            nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
          }
          
          if(Abranks(cb)<colMatBlocksI(0)(cb+1,0)){
            //can use truncated matrices
            for (int d=0; d<(ncov); d++){
              matrix<Type> SArm(blocksize,CppAD::Integer(Abranks(cb)));
              vector<Type> SArmD(blocksize);
              for (int j=0; j<blocksize; j++){ // diagonals
                if(j<Abranks(cb)){
                  SArm(j,j) = exp(Abb(sdcounter));
                  SArmD(j) = 0;
                }else{
                  SArmD(j) = exp(Abb(sdcounter));
                }
                
                sdcounter++;
              }
              
              for (int j=0; j<Abranks(cb); j++){
                for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                  SArm(r,j) = Abb(covscounter);
                  covscounter++;
                }
              }
              //determinant with first term reduced rank and second term remaining diagonal entries
              nll -= SArm.diagonal().array().log().sum() + SArmD.tail(blocksize-SArm.cols()).log().sum();
              
              matrix<Type> SArmpd = SArm*SArm.transpose();
              SArmpd.diagonal() += SArmD.pow(2).matrix();
              //remaining likelihood terms
              //tr(S⁻¹A)
              if((rhoSP.size() == 1)){
                if(nncolMat.rows()<p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatI*SArmpd).trace());
                }else{
                  //sparse approximation to inverse
                  nll -= -0.5*(SprI(d,d)*(colCorMatUIs(0)*colCorMatUIs(0).transpose()*SArmpd).trace());
                }
              }else{
                //sparse approximation to inverse
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArmpd).trace()*SprI(d,d);
              }
              
              matrix<Type>xbxb = (xb.col(d).cwiseProduct(xb.col(d))).rowwise().sum();
              cQ.middleCols(sp, blocksize) += 0.5*xbxb*SArmpd.diagonal().transpose();
            }
          }else{
            for (int d=0; d<(ncov); d++){
              matrix<Type> SArm(blocksize,blocksize);//limit the scope of SArm(d) to this iteration for memory efficiency due to potentially large d
              SArm.setZero();
              for (int j=0; j<blocksize; j++){ // diagonals
                SArm(j,j)=exp(Abb(sdcounter));
                sdcounter++;
              }
              
              for (int j=0; j<blocksize; j++){
                for (int r=j+1; r<(colMatBlocksI(0)(cb+1,0)); r++){
                  SArm(r,j)=Abb(covscounter);
                  covscounter++;
                }
              }
              
              nll -= SArm.diagonal().array().log().sum();
              
              SArm *= SArm.transpose();
              //tr(S⁻¹A)
              if((rhoSP.size() == 1)){
                if(nncolMat.rows()<p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatI*SArm).trace());
                }else if(nncolMat.rows()==p){
                  nll -= -0.5*(SprI(d,d)*(colCorMatUIs(0)*colCorMatUIs(0).transpose()*SArm).trace());
                }
              }else{
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArm).trace()*SprI(d,d);
              }
              
              matrix<Type>xbxb = (xb.col(d).cwiseProduct(xb.col(d))).rowwise().sum();
              cQ.middleCols(sp,blocksize) += 0.5*xbxb*SArm.diagonal().transpose();
            }
          }
          
          sp += blocksize;
          
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 3){
        //(Abb.size() == (p*ncov + p*(p-1)/2 + p-colCorMatIblocks.size())) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + p*ncov*(ncov-1)/2 + p*(p-1)/2)) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + p*ncov*(ncov-1)/2 + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum())) || (Abb.size() == (p*ncov + p-colCorMatIblocks.size() + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // Ab.struct == "CL1" OR "diagonalCL1". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        // always has p*p matrix with fixed diagonals, in combination with diagonal/blockdiagonal matrix for sp*ncov
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov+p-(colMatBlocksI.size()-1);
        int sp = 0;
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          //construct blockdiagonal list
          vector<matrix<Type>> SArmb(blocksize);
          
          for (int j=0; j<blocksize; j++){
            SArmb(j).resize(ncov,ncov);
            SArmb(j).setZero();
            for (int d=0; d<(ncov); d++){ // diagonals of varcov
              SArmb(j)(d,d) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            if((Abb.size()>(p*ncov + p-(colMatBlocksI.size()-1) + p*(p-1)/2)) || (Abb.size()>(p*ncov + p-(colMatBlocksI.size()-1) + (colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))){ // unstructured block Var.cov
              //only for "blockdiagonalsp"
              for (int d=0; d<(ncov); d++){
                for (int r=d+1; r<(ncov); r++){
                  SArmb(j)(r,d) = Abb(covscounter);
                  covscounter++;
                }}
            }
          }
          
          //construct p by p covariance matrix
          if(Abranks(cb)<blocksize){
            //construct as sparse matrix
            //use sparse matrices if we can to speed things up with many species
            matrix<Type>SArmC(blocksize,CppAD::Integer(Abranks(cb)));
            SArmC.setZero();
            vector<Type> SArmCD(blocksize);
            SArmCD.setZero();
            
            // set-up column (block) covariance matrix
            SArmC(0,0) = 1;//first entry fixed for identifiability
            for (int j=1; j<blocksize; j++){
              if(j<Abranks(cb)){
                SArmC(j,j) = exp(Abb(sdcounter));
              }else{
                SArmCD(j) = exp(Abb(sdcounter));
              }
              sdcounter++;
            }
            
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<blocksize; r++){
                SArmC(r,j) =  Abb(covscounter);
                covscounter++;
              }
            }
            
            //p*p matrix in the middle: correlation for identifiability
            matrix<Type> SArmP = SArmC*SArmC.transpose();
            SArmP.diagonal() += SArmCD.pow(2).matrix();
            SArmCD.array() /= SArmP.diagonal().array().cwiseSqrt();
            SArmC.array().colwise() /= SArmP.diagonal().array().cwiseSqrt();
            
            //determinant of this matrix is nsp*sum(log(diag(SArmC)))+2*sum(log(diags(SArmb)))
            nll -=  ncov*SArmC.diagonal().array().log().sum() + ncov*SArmCD.tail(blocksize-SArmC.cols()).log().sum();
            for (int j=0; j<blocksize; j++){
              nll -= SArmb(j).diagonal().array().log().sum();
            }
            
            SArmP = gllvmutils::cov2cor(SArmP);
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              matrix<Type> colCorMatI(blocksize, blocksize);
              if(nncolMat.rows()<p){
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
                matrix<Type>tempMat(ncov,ncov);
                for (int j=0; j<blocksize;j++){
                  tempMat = SprI*SArmb(j);
                  for (int j2=j+1; j2<blocksize;j2++){
                    // temp.diagonal().fill(SArmP(j,j2));
                    nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                  }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
                }
                
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                colCorMatI = colCorMatUI*colCorMatUI.transpose();
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }
              //write trace separately so we don't need to compute the whole product
              //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
              matrix<Type>tempMat(ncov,ncov);
              for (int j=0; j<blocksize;j++){
                tempMat = SprI*SArmb(j);
                for (int j2=j+1; j2<blocksize;j2++){
                  // temp.diagonal().fill(SArmP(j,j2));
                  nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
              }
            }else{
              vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
                colCorMatUIs(d) = colCorMatUIs(d).transpose();//to expose .col below
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              matrix<Type>tempMat(ncov,ncov);
              //this trace cannot be separated further due to element-wise product..
              for (int j=0; j<blocksize;j++){
                for (int j2=j+1; j2<blocksize;j2++){
                  for (int d=0; d<(ncov); d++){
                    for (int d2=0; d2<(ncov); d2++){
                      tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j2)).sum()*SprI(d,d2);
                    }
                  }
                  nll-= -(tempMat.cwiseProduct(SArmb(j)*SArmb(j2).transpose())*SArmP(j,j2)).trace();
                }
                for (int d=0; d<(ncov); d++){
                  for (int d2=d+1; d2<(ncov); d2++){
                    tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j)).sum()*SprI(d,d2);
                    tempMat(d2,d) = tempMat(d,d2);
                  }
                  tempMat(d,d) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d).col(j)).sum()*SprI(d,d);
                }
                nll-= -0.5*(tempMat.cwiseProduct(SArmb(j)*SArmb(j).transpose())).trace();
              }
            }
            
            //remaining likelihood terms
            matrix<Type>SArmB(ncov, ncov);
            for (int j=0; j<blocksize;j++){
              SArmB = SArmb(j)*SArmb(j).transpose();
              for(int i=0; i<n;i++){
                cQ(i, j+sp) += 0.5*((xb.row(i)*SArmB).cwiseProduct(xb.row(i))).sum();
              }
            }
          }else{
            //use sparse matrices if we can to speed things up with many species
            matrix<Type>SArmC(blocksize,blocksize);
            SArmC.setZero();
            // set-up column (block) covariance matrix
            
            SArmC(0,0) = 1;//first entry fixed to 1 for identifiability
            for (int j=1; j<SArmC.cols(); j++){
              SArmC(j,j) = exp(Abb(sdcounter));
              sdcounter++;
            }
            
            for (int j=0; j<SArmC.cols(); j++){
              for (int r=j+1; r<SArmC.cols(); r++){
                SArmC(r,j) = Abb(covscounter);
                covscounter++;
              }
            }
            matrix<Type> SArmP = SArmC*SArmC.transpose();
            SArmC *= SArmP.diagonal().cwiseInverse().cwiseSqrt().asDiagonal();//ID constraint
            
            //determinant of this matrix is 2*nsp*sum(log(diag(SArmC)))+2*sum(log(diags(SArmb)))
            nll -=  ncov*SArmC.diagonal().array().log().sum();
            for (int j=0; j<blocksize; j++){
              nll -= SArmb(j).diagonal().array().log().sum();
            }
            SArmP = gllvmutils::cov2cor(SArmP);
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              matrix<Type> colCorMatI(blocksize, blocksize);
              if(nncolMat.rows()<p){
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }else if(nncolMat.cols()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approxiamtion to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatUI*colCorMatUI.transpose()*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              }
              
              //write trace separately so we don't need to compute the whole product
              //Eigen::DiagonalMatrix<Type, Eigen::Dynamic> temp(ncov);//to hold repeated entries of SArmP, for kron(SArmP, diag(nsp))
              matrix<Type> tempMat(blocksize,blocksize);
              for (int j=0; j<blocksize;j++){
                tempMat = SprI*SArmb(j);
                for (int j2=j+1; j2<blocksize;j2++){
                    // temp.diagonal().fill(SArmP(j,j2));
                    nll -= -(tempMat*SArmb(j2).transpose()).trace()*SArmP(j,j2)*colCorMatI(j,j2);
                }
                nll -= -0.5*(tempMat*SArmb(j).transpose()).trace()*colCorMatI(j,j);
              }
            }else{
              vector<Eigen::SparseMatrix<Type, Eigen::ColMajor>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
                colCorMatUIs(d) = colCorMatUIs(d).transpose();//needed to expose .col below
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              matrix<Type>tempMat(ncov,ncov);
              //this trace cannot be separated further due to element-wise product..
              for (int j=0; j<blocksize;j++){
                for (int j2=j+1; j2<blocksize;j2++){
                  for (int d=0; d<(ncov); d++){
                    for (int d2=0; d2<(ncov); d2++){
                      tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j2)).sum()*SprI(d,d2);
                    }
                  }
                  nll-= -(tempMat.cwiseProduct(SArmb(j)*SArmb(j2).transpose())*SArmP(j,j2)).trace();
                }
                for (int d=0; d<(ncov); d++){
                  for (int d2=d+1; d2<(ncov); d2++){
                    tempMat(d,d2) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d2).col(j)).sum()*SprI(d,d2);
                    tempMat(d2,d) = tempMat(d,d2);
                  }
                  tempMat(d,d) = (colCorMatUIs(d).col(j).transpose()*colCorMatUIs(d).col(j)).sum()*SprI(d,d);
                }
                nll-= -0.5*(tempMat.cwiseProduct(SArmb(j)*SArmb(j).transpose())).trace();
              }
            }

            //remaining likelihood terms
            matrix<Type>SArmB(ncov, ncov);
            for (int j=0; j<blocksize;j++){
              SArmB = SArmb(j)*SArmb(j).transpose();
              for(int i=0; i<n;i++){
                cQ(i, j+sp) += 0.5*((xb.row(i)*SArmB).cwiseProduct(xb.row(i))).sum();
              }
            }
            
          }
          sp += blocksize;
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
        
      }else if(Abstruc == 4){
        // Ab.struct == "CL2". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        // always has k*k matrix with fixed diagonals, in combination with k p*p matrices
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = ncov-1+p*ncov;//(p-(colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()==1).count())*ncov;
        int sp = 0;
        
        matrix<Type> SArmR(ncov, ncov);
        SArmR.setZero();
        // row variance, fix first entry for identifiability (cov2cor)
        SArmR(0,0) = 1;
        for (int d=1; d<(ncov); d++){
          SArmR.diagonal()(d)=exp(Abb(sdcounter));
          sdcounter++;
        }
        
        for (int d=0; d<(ncov); d++){
          for (int r=d+1; r<(ncov); r++){
            SArmR(r,d)=Abb(covscounter);
            covscounter++;
          }}
        
        // to cholesky of correlation matrix
        SArmR.array().colwise() /= SArmR.array().rowwise().norm();
        //determinant first part
        nll -=  p*SArmR.diagonal().array().log().sum();
        //define correlation matrix for further down
        matrix<Type>SArmRc = SArmR*SArmR.transpose();
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          //construct blockdiagonal list
          vector<matrix<Type>> SArmPs(ncov);
          
          //construct p by p covariance matrix
          if(Abranks(cb)<blocksize){
            vector<vector<Type>> SArmCDs(ncov);
            
            for (int d=0; d<(ncov); d++){
              SArmPs(d).resize(blocksize,CppAD::Integer(Abranks(cb)));
              SArmPs(d).setZero();
              SArmCDs(d).resize(blocksize);
              SArmCDs(d).setZero();
              for (int j=0; j<blocksize; j++){
                if(j<Abranks(cb)){
                  SArmPs(d)(j,j) = exp(Abb(sdcounter));
                }else{
                  SArmCDs(d)(j) = exp(Abb(sdcounter));
                }
                sdcounter++;
              }
              for (int j=0; j<Abranks(cb); j++){
                for (int r=j+1; r<blocksize; r++){
                  SArmPs(d)(r,j) =  Abb(covscounter);
                  covscounter++;
                }
              }
            }
            
            //determinant of this matrix is 2*p*sum(log(diag(SArmR)))+2*sum(log(diags(SArmCs)))
            //determinant second part
            //expanding SArmPs so the remaining diagonal entries can be added to the reduced rank part
            for (int d=0; d<ncov; d++){
              nll -= SArmPs(d).diagonal().array().log().sum() + SArmCDs(d).tail(blocksize-SArmPs(d).cols()).log().sum();
            }
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                vector<Type> colCorMatIDiag = colCorMatI.diagonal();
                vector<Type> temp(blocksize);
                matrix <Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
                for (int d=0; d<ncov;d++){
                  tempMat = SArmPs(d).transpose()*colCorMatI;
                  temp = colCorMatIDiag*SArmCDs(d);
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(tempMat*SArmPs(d2)).trace()*SprI(d,d2)*SArmRc(d,d2);
                    nll -= -(temp*SArmCDs(d2)).sum()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(SArmPs(d).transpose()*colCorMatI*SArmPs(d)).trace()*SprI(d,d);
                  nll -= -0.5*(colCorMatIDiag*SArmCDs(d)*SArmCDs(d)).sum()*SprI(d,d)*SArmRc(d,d);
                }
              }else if(nncolMat.rows()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approximation to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                //this is the same as in matrix normal, except the m by m matrix is d-specific
                matrix<Type>tempMat(blocksize,blocksize);
                vector<Type> colCorMatIDiag = colCorMatI.diagonal();
                vector<Type>temp(blocksize);
                vector<Type> temp2(blocksize);
                for (int d=0; d<ncov;d++){
                  tempMat = SArmPs(d).transpose()*colCorMatUI;
                  temp2 = colCorMatIDiag*SArmCDs(d);
                  for (int d2=d+1; d2<ncov;d2++){
                    //reduced rank part: L
                    nll -= -(((tempMat)*(colCorMatUI.transpose()*SArmPs(d))).trace())*SArmRc(d,d2)*SprI(d,d2);
                    //remaining diagonal entries
                    nll -= -(temp2*SArmCDs(d2)).sum()*SArmRc(d,d2)*SprI(d,d2);
                  }
                  //reduced rank part
                  nll -= -0.5*(SArmPs(d).transpose()*colCorMatUI).rowwise().squaredNorm().sum()*SprI(d,d);
                  Type trace = 0;
                  //remaining diagonal entries
                  //need to do this column-wise so I do not need to compute the whole product
                  //similar to above, trace via squared norm
                  for (int j=0; j<blocksize; j++){
                    temp = colCorMatUI.col(j);
                    trace += (temp*SArmCDs(d)).pow(2).sum();
                  }
                  nll -= -0.5*trace*SprI(d,d);//no SArmRc because it has 1s on the diagonal
                  // nll -= -0.5*(SArmP*colCorMatUIBlocks(d).transpose()).trace()*SArmR(d,d)*SprI(d,d);
                }
              }
            }else{
              vector<Eigen::SparseMatrix<Type>> colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIs(d),colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
              }
              
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              //this is the same as in matrix normal, except the m by m matrix is d-specific
              matrix<Type> tempMat(CppAD::Integer(Abranks(cb)), blocksize);
              vector<Type> temp(blocksize);
              for (int d=0; d<ncov;d++){
                tempMat = SArmPs(d).transpose()*colCorMatUIs(d);
                for (int d2=d+1; d2<ncov;d2++){
                  //reduced rank part: L
                  nll -= -(tempMat*(colCorMatUIs(d2).transpose()*SArmPs(d2))).trace()*SArmRc(d,d2)*SprI(d,d2);
                  // apparently faster than "just" computing the diagonal entries via an elementwise product
                  // I suppose Eigen can do some optimization with this
                  vector<Type> diags(blocksize);
                  diags.setZero();
                  for (int j=0; j<blocksize; j++){
                    diags += colCorMatUIs(d).col(j).cwiseProduct(colCorMatUIs(d2).col(j));
                  }
                  nll -= -(SArmCDs(d)*SArmCDs(d2)*diags).sum()*SArmRc(d,d2)*SprI(d,d2);
                }
                //reduced rank part
                nll -= -0.5*(tempMat).rowwise().squaredNorm().sum()*SprI(d,d);
                Type trace = 0;
                //remaining diagonal entries
                //need to do this column-wise so I do not need to compute the whole product
                //similar to above, trace via squared norm
                for (int j=0; j<blocksize; j++){
                  temp = colCorMatUIs(d).col(j);
                  trace += (temp*SArmCDs(d)).pow(2).sum();
                }
                nll -= -0.5*trace*SprI(d,d);//no SArmRc because it has 1s on the diagonal..
                // nll -= -0.5*(SArmP*colCorMatUIBlocks(d).transpose()).trace()*SArmR(d,d)*SprI(d,d);
              }
            }
          
          //species-specific submatrices in A
            matrix<Type>Aspp(ncov,ncov);
            for (int j=0; j<blocksize;j++){
              //need to build the dxd species-specific matrix in A
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  Aspp(d2,d) = ((SArmPs(d).row(j).cwiseProduct(SArmPs(d2).row(j))).sum()+(SArmCDs(d)(j)*SArmCDs(d2)(j)))*SArmRc(d2,d);
                  Aspp(d,d2) = Aspp(d2,d);
                }
                Aspp(d,d) = ((SArmPs(d).row(j).cwiseProduct(SArmPs(d).row(j))).sum()+(SArmCDs(d)(j)*SArmCDs(d)(j)));
              }
              cQ.col(sp+j).noalias() += 0.5*(xb*Aspp).cwiseProduct(xb).rowwise().sum();
            }
          //instead compute cholesky factors of the species-specific submatrices in A
            // matrix<Type>Aspp(ncov, ncov*blocksize);
            // Aspp.setZero();
            // 
            // Type sum_diag = 0;
            // Type sum_off_diag = 0;
            // 
            // for (int j=0; j<blocksize;j++){
            // for (int d = 0; d < ncov; d++) {
            //   Type SArmCDsdj = SArmCDs(d)(j);
            //   sum_diag = SArmPs(d).row(j).array().pow(2).sum() + pow(SArmCDsdj, 2);
            //   if (d > 0) {
            //     sum_diag -= (Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d).cwiseProduct(Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d))).sum();
            //   }
            //   Aspp.block(0,j*ncov,ncov,ncov)(d, d) = sqrt(sum_diag);
            //   
            //   vector<Type> SArmPsdj = SArmPs(d).row(j);
            //   matrix<Type> Asppdd = Aspp.block(0,j*ncov,ncov,ncov).row(d).head(d);
            //   for (int d2 =d+1; d2 < ncov; d2++) {
            //       sum_off_diag = (SArmPs(d2).row(j).array() * SArmPsdj).sum() + SArmCDs(d2)(j) * SArmCDsdj;
            //       sum_off_diag *= SArmRc(d2, d);
            //       if (d > 0) {
            //         sum_off_diag -= (Aspp.block(0,j*ncov,ncov,ncov).row(d2).head(d).cwiseProduct(Asppdd)).sum();
            //       }
            //       Aspp.block(0,j*ncov,ncov,ncov)(d2, d) = sum_off_diag / Aspp.block(0,j*ncov,ncov,ncov)(d, d);
            //   }
            // }
            // }
          }else{
            for (int d=0; d<(ncov); d++){
              SArmPs(d).resize(blocksize,blocksize);
              SArmPs(d).setZero();
              
              for (int j=0; j<blocksize; j++){
                SArmPs(d)(j,j) = exp(Abb(sdcounter));
                sdcounter++;
              }
              for (int j=0; j<blocksize; j++){
                for (int r=j+1; r<blocksize; r++){
                  SArmPs(d)(r,j) =  Abb(covscounter);
                  covscounter++;
                }
              }
            }
            
            //determinant of this matrix is 2*p*sum(log(diag(SArmR)))+2*sum(log(diags(SArmCs)))
            //second part
            for (int d=0; d<ncov; d++){
              nll -= SArmPs(d).diagonal().array().log().sum();
            }
            
            //Br'S⁻¹Br and tr(S⁻¹A)
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                for (int d=0; d<ncov;d++){
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(colCorMatI*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(colCorMatI*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
                }
              }else if(nncolMat.rows()==p){
                Eigen::SparseMatrix<Type> colCorMatUI(colMatBlocksI(cb+1).cols(), colMatBlocksI(cb+1).cols());
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb+1).cols()));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                //sparse approximation to inverse
                nll -= -0.5*(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
                
                //write trace separately so we don't need to compute the whole product
                for (int d=0; d<ncov;d++){
                  for (int d2=d+1; d2<ncov;d2++){
                    nll -= -(colCorMatI*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                  }
                  nll -= -0.5*(colCorMatI*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
                }
              }
            }else{
              vector<Eigen::SparseMatrix<Type>>colCorMatUIs(ncov);
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                colCorMatUIs(d).resize(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUIs(d), colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUIs(d);
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //write trace separately so we don't need to compute the whole product
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  nll -= -(colCorMatUIs(d)*colCorMatUIs(d2).transpose()*SArmPs(d2)*SArmPs(d).transpose()).trace()*SprI(d,d2)*SArmRc(d,d2);
                }
                nll -= -0.5*(colCorMatUIs(d)*colCorMatUIs(d).transpose()*SArmPs(d)*SArmPs(d).transpose()).trace()*SprI(d,d);
              }
            }
            
            //remaining likelihood terms
            matrix<Type>tempMat(ncov,ncov);
            for (int j=0; j<blocksize;j++){
              //need to build the dxd species-specific matrix in A
              for (int d=0; d<ncov;d++){
                for (int d2=d+1; d2<ncov;d2++){
                  tempMat(d2,d) = SArmPs(d).col(j).cwiseProduct(SArmPs(d2).row(j).transpose()).sum()*SArmRc(d2,d);
                  tempMat(d,d2) = tempMat(d2,d);
                }
                tempMat(d,d) = SArmPs(d).col(j).cwiseProduct(SArmPs(d).row(j).transpose()).sum();
              }
              for (int i=0; i<n;i++){
                cQ(i,j+sp) += 0.5*((xb.row(i)*tempMat).cwiseProduct(xb.row(i))).sum();
              }
            }
          }
          sp += blocksize;
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
      }else if(Abstruc == 5){
        // (Abb.size() == ((ncov*p)*(ncov*p)-(ncov*p)*((ncov*p)-1)/2)) || (Abb.size() == (ncov*p+(sum(nsp)*colMatBlocksI(0).col(0).segment(1,colMatBlocksI.size()-1).array()*Abranks-Abranks*(Abranks-1)/2-Abranks).sum()))
        // matrix<Type> SArm(p*ncov,p*ncov);
        // Ab.struct == "unstructured". Only here as reference and should not be used unless "colMat" is present
        // ASSUMED THAT THERE IS A PHYLOGENY WHEN GOING HERE
        
        Type logdetColCorMat = 0;
        int sdcounter = 0;
        int covscounter = p*ncov;
        int sp = 0;
        
        typedef Eigen::Triplet<Type> T;
        
        for(int cb=0; cb<(colMatBlocksI.size()-1); cb++){
          int blocksize = colMatBlocksI(cb+1).cols();
          
          if(Abranks(cb)<blocksize){
            std::vector<T> tripletList;
            //can use sparse matrices
            Eigen::SparseMatrix<Type>SArm(blocksize*ncov,blocksize*ncov);
            //block structure :)
            
            for (int d=0; d<(blocksize*ncov); d++){ // diagonals of varcov
              tripletList.push_back(T(d,d,exp(Abb(sdcounter))));
              sdcounter++;
            }
            
            // unstructured block Var.cov
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(ncov*blocksize); r++){
                tripletList.push_back(T(r,j,Abb(covscounter)));
                covscounter++;
              }
            }
            
            SArm.setFromTriplets(tripletList.begin(),tripletList.end());
            
            nll -= SArm.diagonal().array().log().sum();
            
            matrix<Type> SArmb = SArm*SArm.transpose();
            
            //remaining likelihood terms
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type>SprMNI(blocksize*ncov, blocksize*ncov); // inverse of covariane
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                SprMNI = tmbutils::kronecker(colCorMatI,SprI);
                nll -= -0.5*((SprMNI*SArmb).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                Eigen::SparseMatrix<Type> SprMNI = tmbutils::kronecker(colCorMatI,tmbutils::asSparseMatrix(SprI));
                nll -= -0.5*((SprMNI*SArmb).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              matrix<Type>colCorMatUIs(ncov*blocksize,ncov*blocksize);
              colCorMatUIs.setZero();
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                Eigen::SparseMatrix<Type>colCorMatUI(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUI;
                colCorMatUIs.block(d*blocksize,d*blocksize,blocksize,blocksize) = colCorMatUI;
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //get colCorMat in same order as Ab
              vector<int> permVec(blocksize*ncov);
              Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>perm;
              int k = 0;
              for (int j=0; j<blocksize; j++){
                for (int d=0; d<ncov; d++){
                  permVec(k) = j+d*blocksize;
                  k++;
                }
              }
              perm.indices() = permVec;
              matrix<Type> Ip(blocksize,blocksize);
              Ip.setIdentity();
              Eigen::SparseMatrix<Type> kronSprI= tmbutils::kronecker(SprI,Ip).sparseView();
              Eigen::SparseMatrix<Type> SprMNI = tmbutils::asSparseMatrix(colCorMatUIs)*kronSprI*tmbutils::asSparseMatrix(colCorMatUIs).transpose();
              SprMNI = SprMNI.twistedBy(perm.transpose());
              nll -= -0.5*(SprMNI*SArmb).trace();
            }
            
            for (int j=0; j<blocksize;j++){
              //consume smore memory for some reason, but is faster in n
              cQ.col(j+sp) += 0.5*((xb*SArmb.block(j*ncov,j*ncov,ncov,ncov)).cwiseProduct(xb)).rowwise().sum();
            }
          }else{
            matrix<Type>SArm;
            //block structure :)
            SArm = matrix<Type>(blocksize*ncov,blocksize*ncov);
            SArm.setZero();
            for (int d=0; d<(blocksize*ncov); d++){ // diagonals of varcov
              SArm(d,d)=exp(Abb(sdcounter));
              sdcounter++;
            }
            
            // unstructured block Var.cov
            for (int j=0; j<Abranks(cb); j++){
              for (int r=j+1; r<(ncov*colMatBlocksI(0)(cb+1)); r++){
                SArm(r,j)=Abb(covscounter);
                covscounter++;
              }
            }
            
            nll -= SArm.diagonal().array().log().sum();
            
            SArm = SArm*SArm.transpose();
            
            //remaining likelihood terms
            if((rhoSP.size()==1)){
              if(nncolMat.rows()<p){
                matrix<Type> colCorMatI(blocksize, blocksize);
                gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0));
                matrix <Type> SprMNI = tmbutils::kronecker(colCorMatI,SprI);
                nll -= -0.5*((SprMNI*SArm).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }else if(nncolMat.rows()==p){
                //sparse approximation to inverse
                Eigen::SparseMatrix<Type> colCorMatUI(blocksize, blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, blocksize));
                Eigen::SparseMatrix<Type> colCorMatI = colCorMatUI*colCorMatUI.transpose();
                Eigen::SparseMatrix<Type>SprMNI = tmbutils::kronecker(colCorMatI,tmbutils::asSparseMatrix(SprI));
                nll -= -0.5*((SprMNI*SArm).trace()+(Br.middleCols(sp, blocksize)*colCorMatI*Br.middleCols(sp, blocksize).transpose()*SprI).trace());
              }
            }else{
              matrix<Type>colCorMatUIs(ncov*blocksize,ncov*blocksize);
              colCorMatUIs.setZero();
              //NN sparse approximation
              for (int d=0; d<(ncov); d++){
                Eigen::SparseMatrix<Type>colCorMatUI(blocksize,blocksize);
                gllvmutils::nngp(colCorMatUI, colMatBlocksI(cb+1), logdetColCorMat, rhoSP(d), nncolMat.middleCols(sp, blocksize));
                Br.row(d).middleCols(sp,blocksize) *= colCorMatUI;
                colCorMatUIs.block(d*blocksize,d*blocksize,blocksize,blocksize) = colCorMatUI;
              }
              nll -= -0.5*(Br.middleCols(sp, blocksize)*Br.middleCols(sp, blocksize).transpose()*SprI).trace();
              
              //get colCorMat in same order as Ab
              vector<int> permVec(blocksize*ncov);
              Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>perm;
              int k = 0;
              for (int j=0; j<blocksize; j++){
                for (int d=0; d<ncov; d++){
                  permVec(k) = j+d*blocksize;
                  k++;
                }
              }
              perm.indices() = permVec;
              matrix<Type> Ip(blocksize,blocksize);
              Ip.setIdentity();
              Eigen::SparseMatrix<Type> kronSprI= tmbutils::kronecker(SprI,Ip).sparseView();
              Eigen::SparseMatrix<Type> SprMNI = tmbutils::asSparseMatrix(colCorMatUIs)*kronSprI*tmbutils::asSparseMatrix(colCorMatUIs).transpose();
              SprMNI = SprMNI.twistedBy(perm.transpose());
              nll -= -0.5*(SprMNI*SArm).trace();
            }
            
            //remaining likelihood terms
            for (int j=0; j<blocksize;j++){
              //consume smore memory for some reason, but is faster in n
              cQ.col(j+sp) += 0.5*((xb*SArm.block(j*ncov,j*ncov,ncov,ncov)).cwiseProduct(xb)).rowwise().sum();
            }
          }
          sp += blocksize;
          
        }
        
        //det(S⁻¹)
        if((rhoSP.size() == 1)){
          nll -= 0.5*(p*ncov-p*logdetSpr-ncov*logdetColCorMat);
        }else{
          nll -= 0.5*(p*ncov-p*logdetSpr-logdetColCorMat);
        }
      }
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
    
    // Row/Site effects
    if((random(0)>0)){
      vector<Type> sigma = exp(log_sigma);
      eta += (dr0*r0r).replicate(1,p);//*matrix<Type>(Eigen::MatrixXd::Ones(1,p));
      int dccounter = 0; // tracking used dc entries
      int sigmacounter = 0; // tracking used sigma entries
      
      // One: build Arm, variational covariance matrix for all random effects as a list
      int sdcounter = 0;
      int covscounter = nr.sum();
      for(int re=0; re<nr.size();re++){
        
        //unstructured row cov
        if((lg_Ar.size()>nr.sum() && cstruc(re)>0)){ 
          // do not go here with iid RE
          // unstructured Var.cov
          matrix<Type> Arm(nr(re),nr(re));
          matrix<Type> Sr(nr(re), nr(re));
          Arm.setZero();Sr.setZero();
          
          for (int d=0; d<(nr(re)); d++){ // diagonals of varcov
            Arm(d,d)=exp(lg_Ar(sdcounter));
            sdcounter++;
          }
          
          for (int d=0; d<(nr(re)); d++){
            for (int r=d+1; r<(nr(re)); r++){
              Arm(r,d)=lg_Ar(covscounter);
              covscounter++;
            }}
          
          // add terms to cQ
          matrix<Type> ArmMat = Arm*Arm.transpose();
          cQ += (0.5*(dr0.middleCols(nr.head(re).sum(), nr(re))*ArmMat*dr0.middleCols(nr.head(re).sum(), nr(re)).transpose()).diagonal()).replicate(1,p);
          
          // We build the actual covariance matrix
          // This can straightforwardly be extended to estimate correlation between effects
          matrix <Type> invSr(nr(re),nr(re));invSr.setZero();
          Type logdetSr;
          if(cstruc(re) == 1){ // corAR1
            Sr = gllvm::corAR1(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
            sigmacounter+= 2;
          }else if(cstruc(re) == 3){ // corCS
            Sr = gllvm::corCS(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
            sigmacounter += 2;
          }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
            // Distance matrix calculated from the coordinates for rows
            matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
            matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
            DiSc.setZero();
            DiSc.diagonal().array() += 1/sigma(sigmacounter);
            sigmacounter++;
            dc_scaled = dc(dccounter)*DiSc;
            if(cstruc(re)==2){ // corExp
              Sr = gllvm::corExp(sigma(sigmacounter), Type(0), nr(re), dc_scaled);
              sigmacounter++;
            } else if(cstruc(re)==4) { // corMatern
              Sr = gllvm::corMatern(sigma(sigmacounter), Type(1), sigma(sigmacounter+1), nr(re), dc_scaled);
              sigmacounter += 2;
            }
            dccounter++;
          }
          //TMB's matinvpd function: inverse of matrix with logdet for free
          CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
          logdetSr = res[0];
          invSr = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
          
          if(re==0){
            nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(0,nr(re)).transpose()*(invSr*r0r.col(0).segment(0,nr(re)))).sum());
          }else{
            nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*ArmMat).trace()+(r0r.col(0).segment(nr.head(re).sum(),nr(re)).transpose()*(invSr*r0r.col(0).segment(nr.head(re).sum(),nr(re)))).sum());
          }
          // determinants of each block of the covariance matrix
          nll -= 0.5*(nr(re)-logdetSr);
        }else{
          Eigen::DiagonalMatrix<Type, Eigen::Dynamic> Arm(nr(re));
          matrix<Type> Sr(nr(re), nr(re));Sr.setZero();
          
          for (int d=0; d<(nr(re)); d++){ // diagonals of varcov
            Arm.diagonal()(d)=exp(lg_Ar(sdcounter));
            sdcounter++;
          }
          // add terms to cQ
          cQ += (0.5*(dr0.middleCols(nr.head(re).sum(), nr(re))*Arm*Arm*dr0.middleCols(nr.head(re).sum(), nr(re)).transpose()).eval().diagonal()).replicate(1,p);
          
          // We build the actual covariance matrix
          // This can straightforwardly be extended to estimate correlation between effects
          matrix <Type> invSr(nr(re), nr(re));invSr.setZero();
          Type logdetSr = 0;
          // diagonal row effect
          if(cstruc(re) == 0){
            // inverse and log determinant are straighforwardly available here
            logdetSr = 2*nr(re)*log(sigma(sigmacounter));
            sigmacounter++;
          }else if(cstruc(re) == 1){ // corAR1
            Sr = gllvm::corAR1(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
            sigmacounter+= 2;
          }else if(cstruc(re) == 3){ // corCS
            Sr = gllvm::corCS(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
            sigmacounter += 2;
          }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
            // Distance matrix calculated from the coordinates for rows
            matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
            matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
            DiSc.setZero();
            DiSc.diagonal().array() += 1/sigma(sigmacounter);
            sigmacounter++;
            dc_scaled = dc(dccounter)*DiSc;
            if(cstruc(re)==2){ // corExp
              Sr = gllvm::corExp(sigma(sigmacounter), Type(0), nr(re), dc_scaled);
              sigmacounter++;
            } else if(cstruc(re)==4) { // corMatern
              Sr = gllvm::corMatern(sigma(sigmacounter), Type(1), sigma(sigmacounter+1), nr(re), dc_scaled);
              sigmacounter += 2;
            }
            dccounter++;
          }
          if(cstruc(re)>0){
            //TMB's matinvpd function: inverse of matrix with logdet for free
            CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(Sr));
            logdetSr = res[0];
            invSr = atomic::vec2mat(res,Sr.rows(),Sr.cols(),1);
          }
          
          if(re==0){
            if(cstruc(re)==0){
              nll -= Arm.diagonal().array().log().sum() - 0.5*pow(sigma(sigmacounter-1), -2)*(Arm.diagonal().array().pow(2).sum()+(r0r.col(0).segment(0,nr(re)).transpose()*r0r.col(0).segment(0,nr(re))).sum());              
            }else{
              nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*Arm*Arm).trace()+(r0r.col(0).segment(0,nr(re)).transpose()*(invSr*r0r.col(0).segment(0,nr(re)))).sum());
            }
          }else{
            if(cstruc(re)==0){
              nll -= Arm.diagonal().array().log().sum() - 0.5*pow(sigma(sigmacounter-1), -2)*(Arm.diagonal().array().pow(2).sum()+(r0r.col(0).segment(nr.head(re).sum(),nr(re)).transpose()*r0r.col(0).segment(nr.head(re).sum(),nr(re))).sum());
            }else{
              nll -= Arm.diagonal().array().log().sum() - 0.5*((invSr*Arm*Arm).trace()+(r0r.col(0).segment(nr.head(re).sum(),nr(re)).transpose()*(invSr*r0r.col(0).segment(nr.head(re).sum(),nr(re)))).sum());
            }
          }
          // determinants of each block of the covariance matrix
          nll -= 0.5*(nr(re)-logdetSr);
        }
        
      }
    }
    
    
    // Correlated LVs
    if(num_corlv>0) { //CorLV
      int i,j,d;
      int arank = 2;
      matrix<Type> AQ(num_corlv,num_corlv);
      AQ.setZero(); AQ.diagonal().fill(1.0);
      
      if(ucopy.rows() == nu){
        eta += (dLV*ucopy)*newlamCor;
        
        if(cstruclv==0){
          vector<matrix<Type>> Alvm(nu);
          
          for(int d=0; d<nu; d++){
            Alvm(d).resize(num_corlv,num_corlv);
            Alvm(d).setZero();
          }
          
          // Variational covariance for row effects
          for (int q=0; q<(num_corlv); q++){
            for (int d=0; d<(nu); d++){
              Alvm(d)(q,q)=exp(Au(q*nu+d));
            }
          }
          if((Astruc>0) & (Au.size()>((num_corlv)*nu))){//unstructured cov
            int k=0;
            for (int c=0; c<(num_corlv); c++){
              for (int r=c+1; r<(num_corlv); r++){
                for(int d=0; d<nu; d++){
                  Alvm(d)(r,c)=Au(nu*num_corlv+k*nu+d);
                }
                k++;
              }}
          }
          
          for (d=0; d<nu; d++) {
            nll -= atomic::logdet(Alvm(d)) + 0.5*( - (Alvm(d)*Alvm(d).transpose()).trace() - (ucopy.row(d).matrix()*ucopy.row(d).matrix().transpose()).sum());  // nll -= atomic::logdet(Alvm.col(d).matrix()) + 0.5*( - (Alvm.col(d).matrix()*Alvm.col(d).matrix().transpose()).diagonal().sum() - (ucopy.row(d).matrix()*ucopy.row(d).matrix().transpose()).sum());
            for (j=0; j<p;j++){
              cQ.col(j) += 0.5*dLV.col(d)*((newlamCor.col(j).transpose()*(Alvm(d)*Alvm(d).transpose()))*newlamCor.col(j));
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
          
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated between sites/groups
            if(cstruclv==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
            } else if(cstruclv==3) {// Compound Symm  if(cstruclv==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
            } else {
              DiSc_lv.fill(0.0);
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv==2){// exp decaying
                Slv(q) = gllvm::corExp(Type(1), Type(0), nu, dc_scaled_lv);
              } else if(cstruclv==4) {// Matern
                Slv(q) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,dc_lv.cols())), nu, dc_scaled_lv);
              }
            }
            nll -= 0.5*(nu - atomic::logdet(Slv(q)));
          }
          
          if(Astruc<3){
            
            for(int q=0; q<num_corlv; q++){
              
              // u^T*Sinv*u
              Slvinv = atomic::matinv(Slv(q));
              nll -= - 0.5*( ucopy.col(q).transpose()*(Slvinv*ucopy.col(q)) ).sum();
              
              if(Astruc==0 ){//diagonal A cov
                vector<Type> Atemp(nu);
                Atemp.setZero();
                
                for (d=0; d<(nu); d++){
                  Atemp(d)=exp(Au(q*nu+d));
                  // - tr(Sinv*A)
                  nll -= - 0.5*Slvinv(d,d)*pow(Atemp(d),2);
                }
                // 0.5*lambda_qj*A_qii*lambda_qj
                for (j=0; j<p;j++){
                  cQ.col(j) = cQ.col(j) + 0.5*pow(newlamCor(q,j),2)*(dLV*(Atemp.array()*Atemp.array()).matrix());
                }
                // 0.5*logdet(A)
                nll -= log(Atemp.prod());
                
              } else if((Astruc>0)){
                matrix<Type> Atemp(nu, nu);
                Atemp.setZero();
                
                for (d=0; d<(nu); d++){
                  Atemp(d,d)=exp(Au(q*nu+d));
                }
                int k=0;
                if((Astruc==1) & (Au.size() > nu*num_corlv) ){ // unstructured variational covariance
                  for (d=0; d<nu; d++){
                    for (int r=d+1; r<(nu); r++){
                      Atemp(r,d)=Au(nu*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                } else if((Astruc==2) & (Au.size() > nu*num_corlv)) { // bdNN variational covariance
                  arank = NN.rows();
                  for (int r=0; r<(arank); r++){
                    Atemp(NN(r,0)-1,NN(r,1)-1)=Au(nu*num_corlv+k*num_corlv+q);
                    k++;
                  }
                }
                // REPORT(k);
                
                // 0.5*lambda_qj*A_qii*lambda_qj
                for (j=0; j<p;j++){
                  // cQ.col(j) += 0.5*pow(newlamCor(q,j),2)*(dLV*(Alvm(q)*Alvm(q).transpose()).diagonal().matrix());
                  cQ.col(j) += 0.5*pow(newlamCor(q,j),2)*(dLV*(Atemp*Atemp.transpose()).diagonal().matrix());
                }
                
                // 0.5*logdet(A) -0.5*tr(Sinv*A)
                nll -= log(Atemp.diagonal().prod()) + 0.5*(- (Slvinv*(Atemp*Atemp.transpose())).diagonal().sum());
              }
              
            }
            
          } else if((num_corlv>1) & (Astruc<6)){
            // UNN/Kronecker variational covariance
            matrix<Type> Alvm(nu,nu);
            Alvm.setZero();
            
            for (d=0; d<(nu); d++){
              Alvm(d,d)=exp(Au(d));
            }
            
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
            
            // logdet(A) for triang.mat = prod of diag. elements
            nll -= num_corlv*log(Alvm.diagonal().prod()) + nu*log(AQ.diagonal().prod()); // Moved right after Slv initialization: + 0.5*num_corlv*nu;
            // nll -= num_corlv*log(Alvm.determinant()) + nu*log(AQ.determinant()) + 0.5*num_corlv*nu;
            
            Alvm *= Alvm.transpose();
            AQ *= AQ.transpose();
            
            // tr(Sinv*A) + u^T*Sinv*u
            for(int q=0; q<num_corlv; q++){
              Slvinv = atomic::matinv(Slv(q));
              nll -= 0.5*(- AQ(q,q)*(Slvinv*Alvm).trace()-( ucopy.col(q).transpose()*(Slvinv*ucopy.col(q)) ).sum());
            }
            
            // 0.5*lambda_qj*A_qii*lambda_qj
            for (j=0; j<p;j++){
              cQ.col(j) += 0.5*(dLV*Alvm.diagonal())*((newlamCor.col(j).transpose()*AQ)*newlamCor.col(j));
            }
            
            REPORT(Alvm);
          }
          
        }
      } else {
        
        eta += ucopy*newlamCor;
        vector<matrix<Type> > Slv(num_corlv);
        for(int q=0; q<num_corlv; q++){
          Slv(q).resize(times,times);
          Slv(q).setZero();
        }
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
              int k=0;
              if(Astruc==1){
                for(i=0; i<nu; i++){
                  for (d=0; ((d<arank) && (d<times)); d++){
                    for (int r=d+1; r<(times); r++){
                      Alvm(q)(i*times+r,i*times+d)=Au(nu*times*num_corlv+k*num_corlv+q);
                      k++;
                    }
                  }
                }
              } else if(Astruc==2) { //bdNN var cov
                arank = NN.rows();
                for(i=0; i<nu; i++){
                  for (int r=0; r<(arank); r++){
                    Alvm(q)(i*times+NN(r,0)-1,i*times+NN(r,1)-1)=Au(nu*times*num_corlv+k*num_corlv+q);
                    k++;
                  }
                }
              }
            }
            
            
            for (j=0; j<p;j++){
              for (i=0; i<(times*nu); i++) {
                cQ(i,j) += 0.5*pow(newlamCor(q,j),2)*(Alvm(q).row(i)*Alvm(q).row(i).transpose()).sum();
              }
            }
            // 0.5*logdet(A)
            nll -= log(Alvm(q).diagonal().prod()); //log(Alvm(q).determinant());
            
            // Define covariance matrix
            if(cstruclv==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruclv==3) {// Compound Symm  if(cstruclv==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else {
              DiSc_lv.setZero();
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv==2){// exp decaying
                Slv(q) = gllvm::corExp(Type(1), Type(0), times, dc_scaled_lv);
              } else if(cstruclv==4) {// matern
                Slv(q) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,dc_lv.cols())), times, dc_scaled_lv);
              }
            }
            
            nll -= 0.5*nu*(times - atomic::logdet(Slv(q)));
            
            Slvinv = atomic::matinv(Slv(q));
            matrix <Type> Alvmblock(times,times);
            matrix <Type> ucopyblock(times,1);
            for (i=0; i<nu; i++) {
              Alvmblock.setZero();
              ucopyblock.setZero();
              Alvmblock = Alvm(q).block(i*times,i*times,times,times)*Alvm(q).block(i*times,i*times,times,times).transpose();
              ucopyblock = ucopy.block(i*times,q,times,1);
              nll -=  0.5*(- (Slvinv*Alvmblock).trace()-(ucopyblock.transpose()*Slvinv*ucopyblock).sum());
            }
            
          }
          
          REPORT(Alvm);
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
          nll -= num_corlv*log(Alvm.diagonal().prod()) + times*nu*log(AQ.diagonal().prod()) + 0.5*num_corlv*times*nu;
          // nll -= num_corlv*log(Alvm.determinant()) + times*nu*log(AQ.determinant()) + 0.5*num_corlv*times*nu;
          // Alvm *= Alvm.transpose();
          // AQ *= AQ.transpose();
          
          for (j=0; j<p;j++){
            cQ.col(j) += 0.5*(Alvm*Alvm.transpose()).diagonal().matrix()*((newlamCor.col(j).transpose()*(AQ*AQ.transpose()))*newlamCor.col(j));
          }
          
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated within groups
            // Slv.setZero();
            
            // Define covariance matrix
            if(cstruclv==1){// AR1 covariance
              Slv(q) = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
            } else if(cstruclv==3) {// Compound Symm  if(cstruclv==3)
              Slv(q) = gllvm::corCS(Type(1), rho_lvc(q,0), times);
            } else {
              DiSc_lv.setZero();
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv==2){// exp decaying
                Slv(q) = gllvm::corExp(Type(1), Type(0), times, dc_scaled_lv);
              } else if(cstruclv==4) {// matern
                Slv(q) = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,dc_lv.cols())), times, dc_scaled_lv);
              }
            }
            
            nll -= - 0.5*nu*atomic::logdet(Slv(q));
            Slvinv = atomic::matinv(Slv(q));
            matrix <Type> Alvmblock;
            matrix <Type> ucopyblock;
            for (i=0; i<nu; i++) {
              Alvmblock = Alvm.col(q).matrix().block(i*times,i*times,times,times)*Alvm.col(q).matrix().block(i*times,i*times,times,times).transpose();
              ucopyblock = ucopy.block(i*times,q,times,1);
              nll -=  0.5*(- (AQ.row(q)*AQ.row(q).transpose()).sum()*(Slvinv*Alvmblock).trace() - (ucopyblock.transpose()*(Slvinv*ucopyblock)).sum());
            }
            
          }
          REPORT(Alvm);
        }
        
        
      }
      REPORT(AQ);
    }
    
    vector<Eigen::DiagonalMatrix<Type, Eigen::Dynamic>> D(p);
    
    if(nlvr>0){
      matrix<Type> b_lv2(x_lv.cols(),nlvr);
      b_lv2.setZero();
      
      if((num_lv_c>0) && (random(2)<1)){
        //concurrent ordination terms
        //predictor coefficients for constrained ordination
        b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        eta += x_lv*b_lv2*newlam;
        
      }else if((nlvr>0) && (random(2)>0) && (quadratic > 0)){
        if(num_lv_c>0)b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        if(num_RR>0) b_lv2.rightCols(num_RR) = b_lv.rightCols(num_RR);
      }
      lam = u*newlam;
      
      // Update cQ for linear term
      //Binomial, Gaussian, Ordinal
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          cQ(i,j) += 0.5*(newlam.col(j).transpose()*A(i)*A(i).transpose()*newlam.col(j)).value();
        }
      }
      eta += lam;
      
      if(((quadratic>0) && (nlvr>0)) || ((quadratic>0) && (num_RR>0))){
        //quadratic coefficients for ordination
        //if random rows, add quadratic coefficients for num_RR to D otherwise
        //they go into D_RR below
        //The ordering here is num_lv_c-num_lv-num_RR so that the code works for
        //fixed-effects B and random effects B
        //The order we need to pick them from lambda2 is
        //num_lv_c-num_RR-num_lv however, to ensure everything on the R-side works
        if(((num_lv+num_lv_c+num_RR*random(2))>0)){
          
          for (int j=0; j<p; j++){
            D(j).resize(nlvr);
            D(j).setZero();
          }
          
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
          if(num_lv>0){
            if(lambda2.cols()==1){
              //make sure that num_lv is taken from the middle even with num_RR
              for (int j=0; j<p; j++){
                for (int q=(num_lv_c+num_RR*random(2)); q<(num_lv_c+num_RR*random(2)+num_lv); q++){
                  D(j).diagonal()(q-num_RR*random(2)) = fabs(lambda2(q,0)); //common tolerances model
                }
              }
            }else{
              for (int j=0; j<p; j++){
                for (int q=(num_lv_c+num_RR*random(2)); q<(num_lv_c+num_RR*random(2)+num_lv); q++){
                  D(j).diagonal()(q-num_RR*random(2)) = fabs(lambda2(q,j)); //full quadratic model
                }
              }
            }
          }
        }
        if((num_lv_c>0) && (random(2)<1)){
          //quadratic reduced rank term for concurrent ordination
          for (int j=0; j<p;j++){
            for (int i=0; i<n; i++) {
              eta(i,j) -=  x_lv.row(i)*b_lv2*D(j)*(x_lv.row(i)*b_lv2).transpose();
            }
          }
        }
        
        // do not take this route not with quadratic model, (fixed-effect) constrained LVs and random row-effects.
        if(((nlvr > 0) && (num_lv+num_lv_c)>0) || ((quadratic>0) && (random(2) > 0))){
          //quadratic model approximation
          
          matrix <Type> Acov(nlvr,nlvr);
          
          //Poisson, NB, gamma, exponential,ZIP
          if((family==0)||(family==1)||(family==4)||(family==6)||(family==8)||(family==11)){
            int sign = 1;
            //sign controls whether it's Poisson or other
            if((family>0) && (family != 6)){
              sign = 1;
            }else if((family==0)||(family==6)){
              sign = -1;
            }
            
            matrix<Type> Binv(nlvr,nlvr);
            matrix<Type> Cinv(nlvr,nlvr);
            Type logdetC;
            Type vBinvv;
            matrix <Type> BiQ(nlvr,nlvr);
            matrix <Type> Id(nlvr,nlvr);
            Id.setZero();Id.diagonal().fill(1.0);
            //this implementation does not follow calculation from van der Veen et al. 2021
            //but prevents Acov^-1 via woodbury matrix identity
            //see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
            for (int i=0; i<n; i++) {
              Acov = A(i)*A(i).transpose();
              if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
              for (int j=0; j<p;j++){
                Cinv.setZero();
                Binv.setZero();
                BiQ.setZero();
                Cinv = (Id - 2*sign*Acov*D(j)).inverse();
                Binv = Acov+2*sign*Cinv*Acov*D(j)*Acov;
                BiQ = Id+2*sign*D(j)*Cinv*Acov;//Q*Binv
                
                //the calculation generally prevents having to explicitly invert A*A^t, or having to invert A(i).
                vBinvv = (2*newlam.col(j).transpose()*Cinv*Acov*D(j)*(sign*Acov*newlam.col(j)-2*u.row(i).transpose())+2*sign*u.row(i)*D(j)*Cinv*u.row(i).transpose()).value();
                
                //extra cQ contribution for  XB,e cross term in concurrent model
                if((random(2)<1) && (num_lv_c>0)){
                  vBinvv += (-4*x_lv.row(i)*b_lv2*D(j)*Binv*newlam.col(j)+4*sign*u.row(i)*BiQ*D(j)*(x_lv.row(i)*b_lv2).transpose()+4*x_lv.row(i)*b_lv2*D(j)*Binv*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                  cQ(i,j) -= sign*2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
                }
                
                //-logdetA + logdetB = logdetQ + logdetB = logdetC = det(QB^-1)
                // partialPivLU because BiQ is asymmetric, and it ensures that it is invertible.
                logdetC = BiQ.partialPivLu().matrixLU().diagonal().array().log().sum();
                cQ(i,j) += 0.5*(vBinvv+logdetC) - sign*(D(j)*Acov).trace() - sign*(u.row(i)*D(j)*u.row(i).transpose()).sum();
              }
            }
          }
          // Binomial, Gaussian, Ordinal
          if((family==2)||(family==3)||(family==7)){
            for (int i=0; i<n; i++) {
              Acov = A(i)*A(i).transpose();
              if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
              for (int j=0; j<p;j++){
                cQ(i,j) += (D(j)*Acov*D(j)*Acov).trace() +2*(u.row(i)*D(j)*Acov*D(j)*u.row(i).transpose()).value() - 2*(u.row(i)*D(j)*Acov*newlam.col(j)).value();
                if((num_lv_c>0) && (random(2)<1)){
                  //extra terms for concurrent ordination
                  cQ(i,j) += (2*x_lv.row(i)*b_lv2*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose() -2*x_lv.row(i)*b_lv2*D(j)*Acov*newlam.col(j)+4*u.row(i)*D(j)*Acov*D(j)*(x_lv.row(i)*b_lv2).transpose()).value();
                }
              }
            }
          }
          
          for (int i=0; i<n; i++) {
            Acov = A(i)*A(i).transpose();
            if(random(2)>0 && (num_lv_c+num_RR)>0)Acov += Ab_lvcov(i);
            for (int j=0; j<p;j++){
              eta(i,j) += - (u.row(i)*D(j)*u.row(i).transpose()).sum() - (D(j)*A(i)*A(i).transpose()).trace();
              if((num_lv_c>0) && (random(2)<1)){
                eta(i,j) -= 2*u.row(i)*D(j)*(x_lv.row(i)*b_lv2).transpose();
              }
            }
          }
        }
      }
    }
    
    if(family==0){//poisson
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j)))nll -= dpois(y(i,j), exp(eta(i,j)+cQ(i,j)), true)-y(i,j)*cQ(i,j);
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0r(i)/sigma,2))*random(0);
      }
    } else if((family == 1) && (method<1)){//NB VA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          // nll -= Type(gllvm::dnegbinva(y(i,j), eta(i,j), iphi(j), cQ(i,j)));
          // if(!gllvmutils::isNA(y(i,j)))nll -= y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*cQ(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
          if(!gllvmutils::isNA(y(i,j))){
            nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
            Type log_term_1 = log1p(exp(eta(i,j)-lg_phi(j)));
            Type log_term_2 = log1p(exp(eta(i,j)-cQ(i,j)-lg_phi(j)));
            nll -= (y(i,j)+iphi(j))*(log_term_1-log_term_2)-(y(i,j)+iphi(j))*cQ(i,j);
          }
        }
      }
    } else if ((family == 1) && (method>1)) { // NB EVA
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j))){
            nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
            nll += (((iphi(j)+y(i,j)) / (iphi(j)+exp(eta(i,j)))) * exp(eta(i,j)) - ((iphi(j)+y(i,j))*pow(iphi(j)+exp(eta(i,j)),-2))*pow(exp(eta(i,j)),2)) * cQ(i,j);
          }
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
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
          if(!gllvmutils::isNA(y(i,j))){
            nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(j)-y(i,j));
            nll += cQ(i,j)*Ntrials(j);
            if(Ntrials(j)>1 && (Ntrials(j)>y(i,j))){
              nll -= lgamma(Ntrials(j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(j)-y(i,j)+1.);//norm.const.
            }
          }
        }
      }
      REPORT(mu);
      REPORT(eta);
      REPORT(cQ);
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
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
            mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
            
            mu_prime = mu(i,j) * (1-mu(i,j));
            if(!gllvmutils::isNA(y(i,j))){
              nll -= y(i,j) * eta(i,j) + log(1-mu(i,j));
              nll += mu_prime*cQ(i,j);
            }
          }
        }
      } else if (extra(0) == 1) { // probit
        Type etaP;
        for (int i=0; i<n; i++) {
          for (int j=0; j<p; j++) {
            if(!gllvmutils::isNA(y(i,j))){
              etaP = pnorm_approx(Type(eta(i,j)));   //pnorm funktion approksimaatio
              nll -= y(i,j)*log(etaP) + (1-y(i,j))*log(1-etaP); //
              Type etaD =  dnorm(Type(eta(i,j)), Type(0), Type(1), true);   // log normal density evaluated at eta(i,j)
              nll -= ((y(i,j)*(etaP*exp(etaD)*(-eta(i,j))-pow(exp(etaD),2))*pow(1-etaP,2) + (1-y(i,j))*((1-etaP)*exp(etaD)*eta(i,j)-pow(exp(etaD),2))*pow(etaP,2) )/(etaP*etaP*(etaP*etaP-2*etaP+1)))*cQ(i,j); //Tää toimii ok tähän etaD = (log=true)
            }
          }
        }
      }
    } else if(family==3) {//gaussian
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j)))nll -= (y(i,j)*eta(i,j) - 0.5*eta(i,j)*eta(i,j) - cQ(i,j))/(iphi(j)*iphi(j)) - 0.5*(y(i,j)*y(i,j)/(iphi(j)*iphi(j)) + log(2*iphi(j)*iphi(j))) - log(M_PI)/2;
        }
      }
    } else if(family==4) {//gamma
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j)))nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) )*iphi(j) + log(y(i,j)*iphi(j))*iphi(j) - log(y(i,j)) -lgamma(iphi(j));
        }
      }
    } else if(family==5){ // Tweedie EVA
      //Type ePower = extra(0);
      ePower = invlogit(ePower) + Type(1);
      for (int i=0; i<n; i++) {
        for (int j=0; j<p; j++) {
          if(!gllvmutils::isNA(y(i,j))){
            // Tweedie log-likelihood:
            nll -= dtweedie(y(i,j), exp(eta(i,j)), iphi(j), ePower, true);
            if (y(i,j) == 0) {
              // Hessian-trace part:
              nll += (1/iphi(j)) * (2-ePower)*exp(2*eta(i,j))*exp(-ePower*eta(i,j)) * cQ(i,j);
            } else if (y(i,j) > 0) {
              nll -= (1/iphi(j)) * (y(i,j)*(1-ePower)*exp((1-ePower)*eta(i,j)) - (2-ePower)*exp((2-ePower)*eta(i,j))) * cQ(i,j);
            }
          }
        }
      }
    } else if(family==6){ 
      iphi = iphi/(1+iphi);
      Type pVA;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            if(y(i,j)>0){
              nll -= log(1-iphi(j))+y(i,j)*eta(i,j)-exp(eta(i,j)+cQ(i,j))-lfactorial(y(i,j));
            }else{
              pVA = exp(log(-iphi(j)+1)-exp(eta(i,j)+cQ(i,j))-log((1-iphi(j))*exp(-exp(eta(i,j)+cQ(i,j)))+iphi(j)));
              pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
              pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
              nll -= log(iphi(j))-log(1-pVA);
            }
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
      if (method<1) { // VA
        for (int i=0; i<n; i++) {
          for(int j=0; j<p; j++){
            if(!gllvmutils::isNA(y(i,j))){
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
            }
            //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));//
          }
        }
      } else if (method>1) { // EVA
        if (extra(0)==0) { // logit
          for (int i=0; i<n; i++) {
            for(int j=0; j<p; j++){
              if(!gllvmutils::isNA(y(i,j))){
                int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
                //minimum category
                if(y(i,j)==1){
                  //nll -= -log1plus(exp(eta(i,j))-zetanew(j,0));
                  //nll -= -mfexp(dlogis(zetanew(j,0),eta(i,j),Type(1),1))*cQ(i,j);
                  nll -= -logspace_add(Type(0), eta(i,j)-zetanew(j,0));
                  nll -= -dlogis(zetanew(j,0), eta(i,j), Type(1), 0)*cQ(i,j);
                }else if(y(i,j)==ymaxj){
                  //maximum category
                  int idx = ymaxj-2;
                  //nll -= -log1plus(exp(zetanew(j,idx)-eta(i,j)));
                  //nll -= -mfexp(dlogis(zetanew(j,idx),eta(i,j),Type(1),1))*cQ(i,j);
                  nll -= -logspace_add(Type(0), zetanew(j,idx)-eta(i,j));
                  nll -= -dlogis(zetanew(j,idx), eta(i,j), Type(1), 0)*cQ(i,j);
                }else if(ymaxj>2){
                  for (int l=2; l<ymaxj; l++) {
                    if((y(i,j)==l) && (l != ymaxj)){
                      nll -= logspace_sub(-logspace_add(Type(0), eta(i,j)-zetanew(j,l-1)), -logspace_add(Type(0), eta(i,j)-zetanew(j,l-2)));
                      nll -= -dlogis(zetanew(j,l-2),eta(i,j),Type(1),0)*cQ(i,j);
                      nll -= -dlogis(zetanew(j,l-1),eta(i,j),Type(1),0)*cQ(i,j);
                    }
                  }
                }
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
      if (method<1) {
        for (int i=0; i<n; i++) {
          for(int j=0; j<p; j++){
            if(!gllvmutils::isNA(y(i,j))){
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
          }
          // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0r(i)/sigma,2))*random(0);
        }
      } else if (method>1) {
        if (extra(0)==0) {
          for (int i=0; i<n; i++) {
            for (int j=0; j<p; j++) {
              if (y(i,j)==1) { // min category
                nll -= -logspace_add(Type(0),eta(i,j)-zetanew(0));
                nll -= -dlogis(zetanew(0,0),eta(i,j),Type(1),0)*cQ(i,j);
              } else if(y(i,j)==ymax) { // max category
                int idx = ymax-2;
                nll -= -logspace_add(Type(0),zetanew(idx)-eta(i,j));
                nll -= -dlogis(zetanew(idx),eta(i,j),Type(1),0)*cQ(i,j);
              } else if(ymax>2) {
                for (int l=2; l<ymax; l++) {
                  if ((y(i,j)==l) && (l != ymax)) {
                    //nll(i,j) -= logspace_sub(-log1plus(exp(eta(i,j)-zetanew(0,l-1))),-log1plus(exp(eta(i,j)-zetanew(0,l-2))));
                    nll -= logspace_sub(-logspace_add(Type(0),eta(i,j)-zetanew(l-1)),-logspace_add(Type(0),eta(i,j)-zetanew(l-2)));
                    nll -= (-dlogis(zetanew(l-2),eta(i,j),Type(1),0) - dlogis(zetanew(l-1),eta(i,j),Type(1),0))*cQ(i,j);
                  }
                }
              }
            }
          }
        }
      }
      
    } else if(family==8) {// exp dist
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j))) nll -= ( -eta(i,j) - exp(-eta(i,j)+cQ(i,j))*y(i,j) );
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
          if(!gllvmutils::isNA(y(i,j))){
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
    } else if((family==10) && (method<1))  { // hurdle Beta VA-EVA hybrid
      int truep = (p/2);
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
        for (int j=0; j<truep; j++) {
          if(!gllvmutils::isNA(y(i,j))){
            // define mu, mu' and mu''
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            mu_prime2 = 0.0;
            if (extra(0) == 0) { // logit
              // mu(i,truep+j) = Type(CppAD::CondExpGe(eta(i,truep+j), type(0), 1/(1+exp(-eta(i,truep+j)) ), exp(eta(i,truep+j))/(exp(eta(i,truep+j))+1) ));
              z[0] = eta(i,truep+j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
              
              mu(i,truep+j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              
              z[0] = eta(i,j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
              
              mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              mu_prime = mu(i,j) * (1-mu(i,j));
              mu_prime2 = mu_prime * (1-2*mu(i,j));
              
            } else if (extra(0) == 1) { // probit
              mu(i,truep+j) = pnorm(eta(i,truep+j), Type(0), Type(1));
              mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
              mu_prime = dnorm(eta(i,j), Type(0), Type(1));
              mu_prime2 = (-eta(i,j))*mu_prime;
            }
            
            if(y(i,j)==0){
              nll -= log( 1.0 - mu(i,truep+j) ) - cQ(i,truep+j);
            } else{
              nll -= log( mu(i,truep+j) ) - cQ(i,truep+j);
              
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
      }
      
    } else if ((family==10) && (method>1)) { // hurdle beta EVA
      int truep = (p/2);
      Type mu_prime;
      Type mu_prime2;
      Type mu0_prime;
      Type mu0_prime2;
      
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
        for (int j=0; j<truep; j++) {
          if(!gllvmutils::isNA(y(i,j))){
            // define mu, mu' and mu''
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            mu_prime2 = 0.0;
            mu0_prime = 0.0;
            mu0_prime2 = 0.0;
            if (extra(0) == 0) { // logit
              // mu(i,truep+j) = Type(CppAD::CondExpGe(eta(i,truep+j), type(0), 1/(1+exp(-eta(i,truep+j)) ), exp(eta(i,truep+j))/(exp(eta(i,truep+j))+1) ));
              z[0] = eta(i,truep+j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
              
              mu(i,truep+j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              
              z[0] = eta(i,j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
              
              mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              mu_prime = mu(i,j) * (1-mu(i,j));
              mu_prime2 = mu_prime * (1-2*mu(i,j));
              
              mu0_prime = mu(i,truep+j) * (1-mu(i,truep+j));
              mu0_prime2 = mu0_prime * (1-2*mu(i,truep+j));
              
            } else if (extra(0) == 1) { // probit
              mu(i,truep+j) = pnorm(eta(i,truep+j), Type(0), Type(1));
              mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
              mu_prime = dnorm(eta(i,j), Type(0), Type(1));
              mu_prime2 = (-eta(i,j))*mu_prime;
              
              mu0_prime = dnorm(eta(i,truep+j), Type(0), Type(1));
              mu0_prime2 = (-eta(i,truep+j))*mu0_prime;
            }
            
            if(y(i,j)==0){
              nll -= log( 1.0 - mu(i,truep+j) );
              //nll -= -dlogis(Type(0), eta(i,truep+j), Type(1), 0)*cQ(i,truep+j);
              nll -= -(mu0_prime2 * (1-mu(i,truep+j)) + pow(mu0_prime,2))/pow(1-mu(i,truep+j),2) * cQ(i,truep+j);            
            } else{
              nll -= log( mu(i,truep+j) );
              //nll -= -dlogis(eta(i,truep+j), Type(0.0), Type(1), 0)*cQ(i,truep+j);
              nll -= (mu(i,truep+j)*mu0_prime2-pow(mu0_prime,2))/pow(mu(i,truep+j),2) * cQ(i,truep+j);
              
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
      }
    } else if(family==11){ // ZINB
      iphi = iphi/(1+iphi);
      vector<Type> iphiZINB = exp(lg_phiZINB);
      Type pVA;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            if(y(i,j)>0){
              nll -= log(1-iphi(j))+y(i,j)*(eta(i,j)-cQ(i,j)) - (y(i,j)+iphiZINB(j))*log(iphiZINB(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(y(i,j)+iphiZINB(j)) - iphiZINB(j)*cQ(i,j) + iphiZINB(j)*log(iphiZINB(j)) - lgamma(iphiZINB(j)) -lfactorial(y(i,j));
            }else{
              pVA = exp(log(1-iphi(j))- iphiZINB(j)*log(iphiZINB(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(iphiZINB(j)) - iphiZINB(j)*cQ(i,j) + iphiZINB(j)*log(iphiZINB(j)) - lgamma(iphiZINB(j))-log((1-iphi(j))*exp(- iphiZINB(j)*log(iphiZINB(j)+exp(eta(i,j)-cQ(i,j))) + lgamma(iphiZINB(j)) - iphiZINB(j)*cQ(i,j) + iphiZINB(j)*log(iphiZINB(j)) - lgamma(iphiZINB(j)))+iphi(j)));
              pVA = Type(CppAD::CondExpEq(pVA, Type(1), pVA-Type(1e-12), pVA));//check if pVA is on the boundary
              pVA = Type(CppAD::CondExpEq(pVA, Type(0), pVA+Type(1e-12), pVA));//check if pVA is on the boundary
              nll -= log(iphi(j))-log(1-pVA);
            }
          }
        }
      }
    } else if((family==12) && (method<1)) { // ordered Beta VA-EVA hybrid
      
      matrix <Type> zetanew(p,2);
      zetanew.setZero();
      for(int j=0; j<p; j++){
        zetanew(j,0)= zeta(j);
        if(zeta.size()>p) zetanew(j,1)= exp(zeta(p+j));
      }
      
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
          if(!gllvmutils::isNA(y(i,j))){
            // define mu, mu' and mu''
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            mu_prime2 = 0.0;
            // probit link
            if((y(i,j)==0)){
              // mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
              // nll -= log(pow(1.0 - pnorm(zetanew(j,1) - eta(i,j), Type(0), Type(1)), y(i,j)) * pow(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)),(1-y(i,j)))) - cQ(i,j);
              nll -= (1-y(i,j))*log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
            } else if((y(i,j)==1)){
              nll -= y(i,j)*log(1.0 - pnorm(zetanew(j,1) - eta(i,j), Type(0), Type(1)) ) - cQ(i,j); //
            } else{
              // if (extra(0) == 1) { // probit
              if(zeta.size()>p) {
                nll -= log(pnorm(zetanew(j,1) - eta(i,j), Type(0), Type(1)) - pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
              } else {
                nll -= log(1 - pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
              }
              mu(i,j) = pnorm(eta(i,j), Type(0), Type(1));
              mu_prime = dnorm(eta(i,j), Type(0), Type(1));
              mu_prime2 = (-eta(i,j))*mu_prime;
              // }
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
              
              nll -= dbeta(y(i,j), Type(a[0]), Type(b[0]), 1);
              nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
              nll -= iphi(j) * mu_prime2 * (log(y(i,j)) - log(1-y(i,j))) * cQ(i,j) ;
            }
            
          }
        }
      }
      
    } else if ((family==12) && (method>1)) {  // Ordered beta EVA
      matrix <Type> zetanew(p,2);
      zetanew.setZero();
      for(int j=0; j<p; j++){
        zetanew(j,0)= zeta(j);
        if(zeta.size()>p) zetanew(j,1)= exp(zeta(p+j));
      }
      
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
          if(!gllvmutils::isNA(y(i,j))){
            // define mu, mu' and mu''
            mu(i,j) = 0.0;
            mu_prime = 0.0;
            mu_prime2 = 0.0;
            if((y(i,j)==0)){
                nll -= -logspace_add(Type(0),eta(i,j)-zetanew(j,0));
                nll -= -dlogis(zetanew(j,0), eta(i,j), Type(1), 0)*cQ(i,j);
            } else if((y(i,j)==1)){
              nll -= -logspace_add(Type(0),zetanew(j,1)-eta(i,j));
              nll -= -dlogis(zetanew(j,1), eta(i,j), Type(1), 0)*cQ(i,j);
            } else{
              if(zeta.size()>p) {
                //nll -= log(pnorm(zetanew(j,1) - eta(i,j), Type(0), Type(1)) - pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
                nll -= logspace_sub(-logspace_add(Type(0),zetanew(j,0)-eta(i,j)), -logspace_add(Type(0),zetanew(j,1)-eta(i,j)));
                nll -= -dlogis(zetanew(j,0), eta(i,j), Type(1), 0)*cQ(i,j); 
                nll -= -dlogis(zetanew(j,1), eta(i,j), Type(1), 0)*cQ(i,j);
              } else {
                //nll -= log(1 - pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1))) - cQ(i,j); //
                nll -= eta(i,j) - zetanew(j,0) - logspace_add(Type(0), eta(i,j)-zetanew(j,0));
                nll -= -dlogis(zetanew(j,0), eta(i,j), Type(1), 0)*cQ(i,j);
              }
              CppAD::vector<Type> z(4);
              z[0] = eta(i,j);
              z[1] = 0;
              z[2] = 1/(1+exp(-z[0]));
              z[3] = exp(z[0])/(exp(z[0])+1);
          
              mu(i,j) = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
              mu_prime = mu(i,j) * (1-mu(i,j));
              mu_prime2 = mu_prime * (1-2*mu(i,j));
              
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
              
              nll -= dbeta(y(i,j), Type(a[0]), Type(b[0]), 1);
              nll -= ((-trig_a) * pow(iphi(j)*mu_prime, 2) - dig_a * iphi(j) * mu_prime2 - trig_b * pow(iphi(j)*mu_prime, 2) + dig_b * iphi(j) * mu_prime2) * cQ(i,j);
              nll -= iphi(j) * mu_prime2 * logit(y(i,j)) * cQ(i,j);
            }
            
          }
        }
      }
    }
    
    
  }else{
    using namespace density;
    if(random(2)>0){
      // REPORT(Sigmab_lv); //!!!!
      if(randomB>0){
        MVNORM_t<Type> mvnorm(Sigmab_lv(0));
        for (int klv=0; klv<Klv; klv++) {
          nll += mvnorm(b_lv.row(klv));
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
    // REPORT(dr0);
    // REPORT(cstruc);
    
    // matrix<Type> etaH(n,p); 
    // etaH.setZero();
    
    //For fixed-effects RRR with and without quadratic term
    if(num_RR>0){
      matrix<Type> b_lv3 = b_lv.rightCols(num_RR);
      eta += x_lv*b_lv3*RRgamma;
      if(quadratic>0){
        matrix<Type> D_RR(num_RR,num_RR);
        D_RR.setZero();
        if(lambda2.cols()==1){
          for (int d=0; d<num_RR;d++){
            D_RR.diagonal()(d) = fabs(lambda2(d,0));
          }
          for (int i=0; i<n; i++) {
            eta.row(i).array() -=  (x_lv.row(i)*b_lv3*D_RR*(x_lv.row(i)*b_lv3).transpose()).value();
          }
          
        }else{
          for (int j=0; j<p;j++){
            D_RR.setZero();
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
    // if(r0f.size() == n && (random(0)==0) && xr.rows() != n && r0f.rows() == n && r0f.cols() == 1){
    //   eta += r0f.replicate(1,p);
    if(xr.rows()==n){
      eta += (xr*r0f).replicate(1,p);
    }
    
    if((random(1)>0) | (random(3)>0)){
      if(random(1)>0){
        // random slopes in TMBtrait.R
        eta += xb*Br;  
      }else if(random(3)>0){
        //random slopes in gllvm.TMB
        eta += xb*(Br.colwise()+B.col(0));
      }
      
      matrix <Type> Spr(xb.cols(),xb.cols());
      matrix <Type> SprI(xb.cols(),xb.cols());
      Type logdetSpr = 0;
      
      int l = xb.cols();
      if(random(1)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB);
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL(l,l);
        SprL.fill(0.0);
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaij(i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        matrix <Type> SprIL(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      if(random(3)>0){
        // Eigen::DiagonalMatrix<Type, Eigen::Dynamic>sds(l);
        matrix<Type> sds = Eigen::MatrixXd::Zero(l,l);
        sds.diagonal() =  exp(sigmaB.segment(0,xb.cols()));
        
        vector<Type>sigmaSPij((l*l-l)/2);
        sigmaSPij.fill(0.0);
        //covariances of random effects
        matrix<Type> SprL;
        if(cs.cols()>1){
          //need a vector with covariances and zeros in the right places
          for(int i=0; i<cs.rows(); i++){
            sigmaSPij((cs(i,0) - 1) * (cs(i,0) - 2) / 2 + cs(i,1)-1) = sigmaB(xb.cols()+i);
          }
          SprL = sds*gllvmutils::constructL(sigmaSPij);
        }else{
          SprL = sds;
        }
        matrix <Type> Ir = Eigen::MatrixXd::Identity(xb.cols(),xb.cols());
        matrix <Type> SprIL(xb.cols(),xb.cols());
        SprIL = SprL.template triangularView<Eigen::Lower>().solve(Ir);
        SprI = SprIL.transpose()*SprIL;
        Spr=SprL*SprL.transpose();
        logdetSpr = 2*SprL.diagonal().array().log().sum();
      }
      
      if(colMatBlocksI(0)(0,0)==p){
        Type logdetColCorMat = 0;
        vector <Type> rhoSP(sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1));
        rhoSP.fill(1.0);
        if(sigmaB.size()>(xb.cols()+cs.rows()*(cs.cols()>1))){
          rhoSP = exp(-exp(sigmaB.segment(xb.cols()+cs.rows()*(cs.cols()>1),sigmaB.size()-xb.cols()-cs.rows()*(cs.cols()>1))));
          if(nncolMat.rows()<p){
            //need to cap this for numerical stability
            for(int re=0; re<rhoSP.size(); re++){
              rhoSP(re) = CppAD::CondExpLt(rhoSP(re), Type(1e-12), Type(1e-12), rhoSP(re));
            }
          }
        }
        
        int sp = 0;
        
        if(nncolMat.rows()<p && rhoSP.size()==1){
          //only go here if rhoSP.size()==1. Other case we need a cholesky.
          logdetColCorMat = colMatBlocksI(0).col(1).segment(1,colMatBlocksI.size()-1).sum();
          for(int cb=1; cb<colMatBlocksI.size(); cb++){
            //efficiently update inverse and determinant using rank 1 updates
            matrix<Type> colCorMatI(colMatBlocksI(cb).cols(), colMatBlocksI(cb).cols());
            gllvmutils::rank1inv(colCorMatI, colMatBlocksI(cb), logdetColCorMat, rhoSP(0));
            //formulate MVNORM manually because we have all the components
            nll += 0.5*(Br.middleCols(sp, colCorMatI.cols())*colCorMatI*Br.middleCols(sp, colCorMatI.cols()).transpose()*SprI).trace();
            sp += colCorMatI.cols();
          }
          //add cheap normalizing constant
          nll -= 0.5*(p*xb.cols()-p*logdetSpr-xb.cols()*logdetColCorMat);
        }else if(nncolMat.rows()==p){
          if(rhoSP.size()==1){//p(block) sized matrix
            int sp = 0;
            for(int cb=1; cb<colMatBlocksI.size(); cb++){
              Eigen::SparseMatrix<Type> AL(colMatBlocksI(cb).cols(),colMatBlocksI(cb).cols());
              gllvmutils::nngp(AL, colMatBlocksI(cb), logdetColCorMat, rhoSP(0), nncolMat.middleCols(sp, colMatBlocksI(cb).cols()));
              nll += 0.5*(Br.middleCols(sp, AL.cols())*AL*AL.transpose()*Br.middleCols(sp, AL.cols()).transpose()*SprI).trace();
              sp += colMatBlocksI(cb).cols();
            }
            //add cheap normalizing constant
            nll -= 0.5*(p*xb.cols()-p*logdetSpr-xb.cols()*logdetColCorMat);
            
          }else{//now colCorMatI is updated for every covariate
            int sp = 0;
            for(int cb=1; cb<colMatBlocksI.size(); cb++){
              Eigen::SparseMatrix<Type> AL(colMatBlocksI(cb).cols(), colMatBlocksI(cb).cols());
              for(int re=0; re<rhoSP.size(); re++){
                gllvmutils::nngp(AL, colMatBlocksI(cb), logdetColCorMat, rhoSP(re), nncolMat.middleCols(sp, colMatBlocksI(cb).cols()));
                Br.row(re).middleCols(sp,colMatBlocksI(cb).cols()) *= AL;
              }
              sp += colMatBlocksI(cb).cols();
            }
            nll += 0.5*(Br*Br.transpose()*SprI).trace();
            nll -= 0.5*(p*xb.cols()-p*logdetSpr-logdetColCorMat);
          }
          
        }
      }else{
        //independence across species
        MVNORM_t<Type> mvnorm(Spr);
        for (int j=0; j<p;j++){
          nll += mvnorm(Br.col(j));
        }
      }
    }
    
    //latent variables
    if(nlvr>0){
      for (int i=0; i<n; i++) {
        for(int q=0; q<u.cols(); q++){
          nll -= dnorm(u(i,q), Type(0), Type(1), true);
        }
      }
      //variances of LVs
      u *= Delta;
      if(num_lv_c>0){
        matrix<Type> b_lv2(x_lv.cols(),nlvr);
        
        b_lv2.leftCols(num_lv_c) = b_lv.leftCols(num_lv_c);
        eta += x_lv*b_lv2*newlam;
      }
      // add LV term to lin. predictor 
      lam += u*newlam;
      eta += lam;
      // if(family==10){
      //   // etaH += lam;
      //   etaH += ucopy*thetaH;
      // }
    }
    
    
    // Row/Site effects
    if((random(0)>0)){
      vector<Type> sigma = exp(log_sigma);
      eta += (dr0*r0r).replicate(1,p);//matrix<Type>(Eigen::MatrixXd::Ones(1,p));
      
      int dccounter = 0; // tracking used dc entries
      int sigmacounter = 0; // tracking used sigma entries
      for(int re=0; re<nr.size();re++){
        matrix<Type> Sr(nr(re),nr(re));Sr.setZero();
        
        // diagonal row effect
        if(cstruc(re) == 0){
          Sr.diagonal().array() = pow(sigma(sigmacounter), 2);
          sigmacounter++;
        }else if(cstruc(re) == 1){ // corAR1
          Sr = gllvm::corAR1(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
          sigmacounter+=2;
        }else if(cstruc(re) == 3){ // corCS
          Sr = gllvm::corCS(sigma(sigmacounter), log_sigma(sigmacounter+1), nr(re));
          sigmacounter += 2;
        }else if((cstruc(re) == 4) || (cstruc(re) == 2)){ // corMatern, corExp
          // Distance matrix calculated from the coordinates for rows
          matrix<Type> DiSc(dc(dccounter).cols(),dc(dccounter).cols()); DiSc.fill(0.0);
          matrix<Type> dc_scaled(dc(dccounter).rows(),dc(dccounter).cols()); dc_scaled.fill(0.0);
          DiSc.setZero();
          DiSc.diagonal().array() += 1/sigma(sigmacounter);
          sigmacounter++;
          dc_scaled = dc(dccounter)*DiSc;
          if(cstruc(re)==2){ // corExp
            Sr = gllvm::corExp(sigma(sigmacounter), Type(0), nr(re), dc_scaled);
            sigmacounter++;
          } else if(cstruc(re)==4) { // corMatern
            Sr = gllvm::corMatern(sigma(sigmacounter), Type(1), sigma(sigmacounter+1), nr(re), dc_scaled);
            sigmacounter += 2;
          }
          dccounter++;
        }
        
        if(cstruc(re)==0){
          //independence of REs
          vector<Type> r0s = r0r.col(0).segment(0,nr(re));
          for(int ir=0; ir<nr(re); ir++){
            nll -= dnorm(r0s(ir), Type(0), sigma(sigmacounter-1), true);
          }
        }else{
          nll += MVNORM(Sr)(r0r.col(0).segment(0,nr(re)));
        }
        
      }
    }
    
    // Correlated LVs
    if(num_corlv>0) {
      int i;
      if(ucopy.rows() == nu){
        eta += (dLV*ucopy)*newlamCor;
        // if(family==10){ // betaH
        //   etaH += (dLV*ucopy)*thetaH;
        // }
        
        // group specific lvs
        if(cstruclv==0){// no covariance
          matrix<Type> Slv(num_corlv,num_corlv);
          Slv.setZero();
          Slv.diagonal().fill(1.0);
          MVNORM_t<Type> mvnorm(Slv);
          for (int i=0; i<nu; i++) {
            nll += mvnorm(ucopy.row(i));
          }
          // REPORT(Slv);
        } else {
          
          matrix<Type> Slv(nu,nu);
          for(int q=0; q<num_corlv; q++){
            // site specific LVs, which are correlated between groups
            Slv.setZero();
            
            if(cstruclv==1){// AR1 covariance
              Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), nu);
            } else if(cstruclv==3) {// Compound Symm  if(cstruclv==3)
              Slv = gllvm::corCS(Type(1), rho_lvc(q,0), nu);
            } else {
              DiSc_lv.setZero();
              for(int j=0; j<dc_lv.cols(); j++){
                DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
              }
              dc_scaled_lv = dc_lv*DiSc_lv;
              if(cstruclv==2){// exp decaying
                Slv = gllvm::corExp(Type(1), Type(0), nu, dc_scaled_lv);
                // Slv = gllvm::corExp(Type(1), (rho_lvc(q,0)), nu, DistM);
              } else if(cstruclv==4) {// matern
                Slv = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,dc_lv.cols())), nu, dc_scaled_lv);
              }
            }
            
            MVNORM_t<Type> mvnormS1(Slv);
            nll += mvnormS1(ucopy.col(q));
          }
          // REPORT(Slv);
        }
      } else {
        
        matrix<Type> Slv(times,times);
        eta += ucopy*newlamCor;
        // if(family==10){// betaH
        //   etaH += ucopy*thetaH;
        // }
        for(int q=0; q<num_corlv; q++){
          // site specific LVs, which are correlated within groups
          Slv.setZero();
          // Define covariance matrix
          if(cstruclv==1){// AR1 covariance
            Slv = gllvm::corAR1(Type(1), rho_lvc(q,0), times);
          } else if(cstruclv==3) {// Compound Symm  if(cstruclv==3)
            Slv = gllvm::corCS(Type(1), rho_lvc(q,0), times);
          } else {
            DiSc_lv.setZero();
            for(int j=0; j<dc_lv.cols(); j++){
              DiSc_lv(j,j) += 1/exp(rho_lvc(q,j));
            }
            dc_scaled_lv = dc_lv*DiSc_lv;
            if(cstruclv==2){// exp decaying
              Slv = gllvm::corExp(Type(1), Type(0), times, dc_scaled_lv);
            } else if(cstruclv==4) {// matern
              Slv = gllvm::corMatern(Type(1), Type(1), exp(rho_lvc(q,dc_lv.cols())), times, dc_scaled_lv);
            }
          }
          
          MVNORM_t<Type> mvnormS2(Slv);
          
          for (i=0; i<nu; i++) {
            nll += mvnormS2(ucopy.block(i*times,q,times,1));
          }
        }
        // REPORT(Slv);
      }
      // REPORT(nu);
    }
    
    if(model<1){
      // gllvm.TMB.R
      // if(family==10){
      //   etaH += x*bH;
      // }
      eta += x*b;
      for (int j=0; j<p; j++){
        for(int i=0; i<n; i++){
          mu(i,j) = exp(eta(i,j));
        }
      }
      
    } else {
      // Fourth corner model, TMBtrait.R
      // if(family==10){
      //   matrix<Type> eta1h=x*bH;
      //   eta1h.resize(n, p);
      //   etaH += eta1h;
      // }
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
          if(!gllvmutils::isNA(y(i,j)))nll -= dpois(y(i,j), exp(eta(i,j)), true);
        }
      }
    } else if(family==1){//negative.binomial family
      if((num_RR>0) && (nlvr == 0) && (random(2)<1)){
        //use dnbinom_robust in this case - below code does not function well
        //for constrained ordination without any random-effects
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j)))nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1);
          }
        }
      }else{
        for (int j=0; j<p;j++){
          for (int i=0; i<n; i++) {
            if(!gllvmutils::isNA(y(i,j)))nll -= y(i,j)*(eta(i,j)) - y(i,j)*log(iphi(j)+mu(i,j))-iphi(j)*log(1+mu(i,j)/iphi(j)) + lgamma(y(i,j)+iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
          }
        } 
      }
      // } else if(family==2) {//binomial family
      //   for (int j=0; j<p;j++){
      //     for (int i=0; i<n; i++) {
      //       if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
      //       } else {mu(i,j) = pnorm(eta(i,j));}
      //       nll -= log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));
      //     }
      //   }
      
    } else if(family==2) {//binomial family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
          } else {mu(i,j) = pnorm(eta(i,j));}
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(1), mu(i,j)-Type(1e-12), mu(i,j)));//check if on the boundary
          mu(i,j) = Type(CppAD::CondExpEq(mu(i,j), Type(0), mu(i,j)+Type(1e-12), mu(i,j)));//check if on the boundary
          if(!gllvmutils::isNA(y(i,j))){
            nll -= y(i,j)*log(mu(i,j))+log(1-mu(i,j))*(Ntrials(j)-y(i,j));
            if(Ntrials(j)>1 && (Ntrials(j)>y(i,j))){
              nll -= lgamma(Ntrials(j)+1.) - lgamma(y(i,j)+1.) - lgamma(Ntrials(j)-y(i,j)+1.);//norm.const.
            }
          }
        }
      }
    } else if(family==3){//gaussian family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))) nll -= dnorm(y(i,j), eta(i,j), iphi(j), true); 
        }
      }
    } else if(family==4){//gamma family
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= dgamma(y(i,j), iphi(j), exp(eta(i,j))/iphi(j), true); 
        }
      }
    } else if(family==5){//tweedie familyF
      ePower = invlogit(ePower) + Type(1);
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))) nll -= dtweedie(y(i,j), exp(eta(i,j)),iphi(j),ePower, true); 
        }
      }
    } else if(family==6) {//zero-infl-poisson
      iphi=iphi/(1+iphi);
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j)))nll -= dzipois(y(i,j), exp(eta(i,j)),iphi(j), true); 
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
          if(!gllvmutils::isNA(y(i,j))){
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
        }
        // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0r(i)/sigma,2))*random(0);
      }
    } else if(family==8) {// exponential family
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(!gllvmutils::isNA(y(i,j)))nll -= dexp(y(i,j), exp(-eta(i,j)), true);  // (-eta(i,j) - exp(-eta(i,j))*y(i,j) );
        }
      }
    } else if(family==9) {// beta family
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          if(extra(0)<1) {mu(i,j) = mu(i,j)/(mu(i,j)+1);
          } else {mu(i,j) = pnorm(eta(i,j));}
          if(!gllvmutils::isNA(y(i,j)))nll -= dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
        }
      }
    } else if(family==10) {// beta hurdle family
      int truep = (p/2);
      for (int i=0; i<n; i++) {
        for (int j=0; j<truep; j++){
          if(extra(0)<1) {
            // etaH(i,j) = exp(etaH(i,j))/(exp(etaH(i,j))+1);
            mu(i,j) = mu(i,j)/(mu(i,j)+1);
            mu(i,truep+j) = mu(i,truep+j)/(mu(i,truep+j)+1);
          } else {
            // etaH(i,j) = pnorm(etaH(i,j));
            mu(i,j) = pnorm(eta(i,j));
            mu(i,truep+j) = pnorm(eta(i,truep+j));
          }
          if(!gllvmutils::isNA(y(i,j))){
            if (y(i,j) == 0) {
              // nll -= log(1-mu(i,j));
              nll -= log(1-mu(i,truep+j));
            } else{
              // nll -= log(mu(i,j)) + dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
              nll -= log(mu(i,truep +j)) + dbeta(squeeze(y(i,j)), Type(mu(i,j)*iphi(j)), Type((1-mu(i,j))*iphi(j)), 1);
            }
          }
        }
      }
      // REPORT(mu);
      // REPORT(etaH);
    } else if(family==11) {//zero-infl-NB
      iphi=iphi/(1+iphi);
      // vector<Type> iphiZINB = exp(lg_phiZINB);
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          if(!gllvmutils::isNA(y(i,j))){
            if(y(i,j)>0){
              nll -= log(1-iphi(j)) + dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phiZINB(j), 1);
            }else{
              nll -= log(iphi(j) + (Type(1)-iphi(j))*dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phiZINB(j), 0)); 
            }
          }
        }
      }
    }
  }
  return nll;
}
