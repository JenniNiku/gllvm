// some functions
namespace gllvmutils{
// 
// TMB_ATOMIC_VECTOR_FUNCTION(
//   Hypo,
//   1,
//   // ATOMIC_DOUBLE
//   {
//     // tx, ty are double*
//     double x = tx[0];
//     double y = tx[1];
//     double abs_x = std::abs(x);
//     double abs_y = std::abs(y);
//     double max_val = (abs_x > abs_y) ? abs_x : abs_y;
//     double r = (abs_x > abs_y) ? (abs_y / abs_x) : (abs_x / abs_y);
//     if(max_val == 0) {
//       ty[0] = 0;
//     } else {
//       ty[0] = max_val * std::sqrt(1 + r * r);
//     }
//   },
//   // ATOMIC_REVERSE
//   {
//     // In here, use "Type" as the templated scalar
//     Type x = tx[0];
//     Type y = tx[1];
//     Type abs_x = fabs(x);
//     Type abs_y = fabs(y);
//     Type max_val = (abs_x > abs_y) ? abs_x : abs_y;
//     Type r = (abs_x > abs_y) ? (abs_y / abs_x) : (abs_x / abs_y);
//     Type denom;
//     if(max_val == Type(0)) {
//       denom = Type(0);
//     } else {
//       denom = max_val * sqrt(Type(1) + r * r);
//     }
//     if(denom == Type(0)) {
//       px[0] = Type(0);
//       px[1] = Type(0);
//     } else {
//       px[0] = py[0] * x / denom;
//       px[1] = py[0] * y / denom;
//     }
//   }
// );
// 
// template<class Type>
// Type Hypo(const Type &x, const Type &y) {
//   vector<Type> tx(2);
//   tx[0] = x;
//   tx[1] = y;
//   vector<Type> ty(1);
//   Hypo(tx, ty);
//   return ty[0];
// }
// 
// template<class Type>
// Type hypo(const Type &x, const Type &y) {
//   CppAD::vector<Type> tx(2);
//   tx[0] = x;
//   tx[1] = y;
//   CppAD::vector<Type> ty = Hypo(tx);
//   return ty[0];
// }

//structure to pass list from R as list in TMB
template<class Type>
struct dclist : vector<matrix<Type>> {
  dclist(){}
  
  dclist(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

//nested list from R into TMB
template<class Type>
struct nesteddclist : vector<dclist<Type>> {
  nesteddclist(SEXP x) {  // x = List of Lists from R
    this->resize(LENGTH(x));
    for(int i = 0; i < LENGTH(x); i++) {
      SEXP inner_list = VECTOR_ELT(x, i);
      (*this)[i] = dclist<Type>(inner_list);  // use existing constructor
    }
  }
};

//covariance to correlation
template<class Type>
matrix<Type> cov2cor(matrix<Type> x){
  return x.diagonal().cwiseInverse().cwiseSqrt().asDiagonal()*x*x.diagonal().cwiseInverse().cwiseSqrt().asDiagonal();
}

//check for NA
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Constructs lower cholesky of correlation matrix from off-diagonal entries
// need this because unstructured_corr does not return cholesky
// which we need for easy access to the determinant and for phylo effects
template<class Type>
matrix<Type> constructL(const vector<Type> &theta){
  int n = (1 + sqrt(1 + 8 *  theta.size())) / 2;
  
  matrix<Type> L = Eigen::MatrixXd::Zero(n, n);  // Initialize Cholesky factor
  L.diagonal().fill(1.0);
  int covscounter = 0;
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      L(i, j) = theta[covscounter];
      covscounter++;
    }
  }
  for (int i = 1; i < n; ++i) {
    L.row(i) /= L.row(i).norm();
  }
  return L;
}

//inverts a matrix of form C*rho+(1-rho)I using repeated Sherman-Morrison formula
template<class Type>
void rank1inv(matrix<Type> &covMatI, const matrix<Type> matI, Type &logdet, const Type rho){
  covMatI = matI/rho;
  logdet += matI.cols()*log(rho);//second component for log-determinant of colMat*rho
  matrix<Type> temp(matI.cols(),1);
  for(int j=0; j<matI.cols(); j++){
    logdet += log1p((1-rho)*covMatI(j,j));
    temp.col(0) = covMatI.col(j);//need to evaluate this in a temporary to prevent aliasing issues in the next line
    covMatI.template selfadjointView<Eigen::Lower>().rankUpdate(temp, -(1-rho)/(1+(1-rho)*temp(j,0))); //only update lower half of matrix
    covMatI = covMatI.template selfadjointView<Eigen::Lower>(); //return as symmetric matrix
  }
}

//approximately inverts a matrix of form C*rho+(1-rho)I using nearest neighbours
//construct the sparse (approximate) inverse of upper triangular factor of the covariance matrix
//constructing upper triangular instead of lower for efficient memory allocation
//each column then has #nearest neighbour number of non-zero entries
//also means fewer transposes in the code.
template<class Type>
void nngp(Eigen::SparseMatrix<Type> &covMatLI, const matrix<Type> matI, Type &logdet, const Type rho, const matrix<int> neighbours){
  if(matI.cols()>1){
    vector<int>k = (neighbours.array()>0).cast<int>().colwise().sum();
    covMatLI.reserve(k+1);//+1 for diagonal
    covMatLI.insert(0,0) = 1;
    //starting loop at 1 to skip first entry, which is 1
    //as this code is for a correlation matrix
    Type d = 0;
    for(int j=1; j<matI.cols(); j++){
      matrix <Type> C = matI(neighbours.col(j).head(k(j)).array()-1,neighbours.col(j).head(k(j)).array()-1);
      C *= rho;
      C.diagonal().fill(1.0);
      matrix<Type>C2(k(j),1);
      C2.col(0) = matI(neighbours.col(j).head(k(j)).array()-1, j)*rho;
      matrix<Type>b = C.ldlt().solve(C2);
      d = sqrt(1.0/(1.0-(C2.transpose()*b).value()));
      for(int i=0; i<k(j); i++){
        covMatLI.insert(neighbours(i,j)-1, j) = -b(i,0)*d;
      }
      covMatLI.insert(j,j) = d;
      logdet -= 2*log(d);
    }
    covMatLI.makeCompressed();
  }else{ //no neighbours, must be only 1 species in this "block", logdeterminant is zero
    covMatLI.insert(0,0) = 1/matI(0,0);//=1
  }
  // logdet -= 2*covMatLI.diagonal().array().log().sum();
}

// template<class Type>
// vector<Type> assistNNGP(vector<Type>input){
//   Type rho = input(0);
//   int k = CppAD::Integer(input(input.size()-1));
//   matrix<Type> C(k,k);
//   for(int l=1; l<(k*k+1); l++){
//     C(l-1) = input(l);
//   }
//   vector<Type>res(k);
//   C *= rho;
//   C.diagonal().fill(1.0);
//   
//   matrix <Type> b(k,k);
//   if(k<=4)//Eigen can probably do an efficient/analytic inverse here
//   {
//     b = C.inverse();
//   }else{//invert larger matrices via cholesky
//     CppAD::vector<Type> res1 = atomic::invpd(atomic::mat2vec(C));
//     b = atomic::vec2mat(res1, k, k, 1);
//   }
//   res = b.col(k-1)/sqrt(b(k-1,k-1));
//   return(res);
// }
// REGISTER_ATOMIC(assistNNGP)
//   
//   
//   //approximately inverts a matrix of form C*rho+(1-rho)I using nearest neighbours
//   //construct the sparse (approximate) inverse of upper triangular factor of the covariance matrix
//   //constructing upper triangular instead of lower for efficient memory allocation
//   //each column then has #nearest neighbour number of non-zero entries
//   //also means fewer transposes in the code.
//   template<class Type>
//   void nngp(Eigen::SparseMatrix<Type> &covMatLI, const matrix<Type> matI, Type &logdet, const Type rho, const matrix<int> neighbours){
//     // typedef Eigen::Triplet<Type> T;
//     // std::vector<T> AL;
//     if(matI.cols()>1){
//       // AL.push_back(T(0,0,1));
//       vector<int>k = (neighbours.array()>0).cast<int>().colwise().sum();
//       covMatLI.reserve(k);
//       covMatLI.insert(0,0) = 1;
//       //starting loop at 1 to skip first entry of AL, which is 1
//       //as this code is for a correlation matirx
//       for(int j=1; j<matI.cols(); j++){
//         matrix <Type> C = matI(neighbours.col(j).head(k(j)).array()-1,neighbours.col(j).head(k(j)).array()-1);
//         // C *= rho;
//         // C.diagonal().fill(1.0);
//         // 
//         // matrix <Type> b(k(j),k(j));
//         // if(k(j)<=4)//Eigen can probably do an efficient/analytic inverse here
//         // {
//         //   b = C.inverse();
//         // }else{//invert larger matrices via cholesky
//         //   CppAD::vector<Type> res = atomic::invpd(atomic::mat2vec(C));
//         //   b = atomic::vec2mat(res,k(j),k(j),1);
//         // }
//         //AL(j,(neighbours.col(j)).head(k).array()-1) = b.col(k-1)/sqrt(b(k-1,k-1));
//         vector<Type> input(2+k(j)*k(j));
//         input(0) = rho;
//         for(int i=1; i<(k(j)*k(j)+1); i++){
//           input(i)= C(i-1);
//         }
//         input(input.size()-1) = k(j);
//         // input << rho, atomic::mat2vec(C), k(j);
//         vector<Type>b = assistNNGP(input);
//         for(int i=0; i<k(j); i++){
//           // AL.push_back(T(j,neighbours(i,j)-1,b(i, k(j)-1)/sqrt(b(k(j)-1,k(j)-1))));
//           covMatLI.insert(neighbours(i,j)-1, j) = b(i);
//         }
//       }
//       covMatLI.makeCompressed();
//     }else{ //no neighbours, must be only 1 species in this "block", logdeterminant is zero
//       // AL.push_back(T(0,0,1/matI(0,0)));
//       covMatLI.insert(0,0) = 1/matI(0,0);
//     }
//     // covMatLI.setFromTriplets(AL.begin(),AL.end());
//     
//     logdet -= 2*covMatLI.diagonal().array().log().sum();
//   }

// robust exponential function, similar to the one in ADMB
template<class Type>
Type mfexp(Type x) {
  Type one_true = exp(x);
  Type one_false = exp(Type(-30.))*(1.-x-30.)/(1.+2.*(-x-30.));
  Type two_false = exp(Type(30.))*(1.+2.*(x-30.))/(1.+x-30.);
  
  Type two_true = CppAD::CondExpLt(Type(-30.), x, one_true, one_false);
  Type two_value = CppAD::CondExpLt(x, Type(30.), two_true, two_false);
  return two_value;
}

double log1plus(double x) {
  return log1p(x);  
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  log1plus
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = log1plus(tx[0])
  ,
  // ATOMIC_REVERSE
  px[0] = py[0]/(tx[0]+1);
)

// more stable way to compute log(1+x)
template<class Type>
Type log1plus(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return log1plus(tx)[0];
}

double expminus1(double x) {
  return expm1(x);  
}
  
TMB_ATOMIC_VECTOR_FUNCTION(
    // ATOMIC_NAME
    expminus1
    ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty[0] = expminus1(tx[0])
    ,
    // ATOMIC_REVERSE
    px[0] = py[0]*exp(tx[0]);
)

// robust way to compute exp(x)-1      
template<class Type>
Type expminus1(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return expminus1(tx)[0];
}
    
VECTORIZE1_t(expminus1)

template<class Type>
Type logsinh(Type x) {
  return x - Type(log(2.0)) + log(-expminus1(-2*x));
}

template<class Type>
Type pgtilt(Type m, Type v, int K, int sgn, Type tDD, Type hDh, Type hDADh) {
  Type T = m*m + v;
  Type s2 = v - 2*tDD;
  Type out = 0;
  for(int k=1; k<=K; k++){
    Type d = Type((2*k-1)*(2*k-1))*Type(M_PI*M_PI);
    Type kap = 1/(d + T);
    Type den = 1 + 2*kap*s2;
    Type t = 2*kap*m/den;
    Type c = 2*kap/den;
    Type sh = s2 + 2*Type(sgn)*t*hDh;
    Type mr = m - t*s2 + Type(sgn)*(c - t*t)*hDh;
    Type vr = s2 + 4*Type(sgn)*t*hDh + 4*t*t*hDADh - c*sh*sh + 2*tDD - 4*c*hDADh + 2*c*c*hDh*hDh;
    Type delta = kap*(m*m + s2) - 0.5*log(den) - kap*m*m/den + log1plus((mr*mr + vr)/d) - log1plus(T/d);
    out += CppAD::CondExpLt(delta, Type(0), delta, Type(0));
  }
  return out;
}

extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  void Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
}

double logit_pnorm(double x) {
  double log_p_lower, log_p_upper;
  Rf_pnorm_both(x, &log_p_lower, &log_p_upper, 2 /* both tails */, 1 /* log_p */);
  return log_p_lower - log_p_upper;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  logit_pnorm
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = logit_pnorm(tx[0])
  ,
  // ATOMIC_REVERSE
  Type tmp3 = CppAD::CondExpLe(ty[0], Type(18.0), -ty[0] + 2*log1plus(exp(ty[0])), 2*exp(-ty[0]));
  Type tmp4 = dnorm(tx[0], Type(0), Type(1), true) + tmp3;
  px[0] = exp( tmp4 ) * py[0];
)

// used for more robust version of probit
template<class Type>
Type logit_pnorm(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return logit_pnorm(tx)[0];
}

VECTORIZE1_t(logit_pnorm)

template<class Type>
Type probitscale(Type m, Type v, int sgn, Type tDD, Type hDh, int bnd) {
  Type s2 = v - 2*tDD;
  Type kap = 1/sqrt(6 + m*m + v);
  Type nu = (m + kap*s2 - Type(sgn)*kap*kap*hDh)*exp(-0.5*kap*kap*s2);
  Type delta = logspace_add(Type(0), -logit_pnorm(m)) - logspace_add(Type(0), -logit_pnorm(nu)) + 0.5*(nu*nu - m*m) - kap*kap*s2;
  if(bnd>1){
    return CppAD::CondExpGt(delta, Type(0), delta, Type(0));
  }else{
    return delta;
  }
}

template<class Type>
Type pnormlogup(Type z) {
  Type za = CppAD::CondExpLt(z, Type(-30), z, Type(-30));
  Type zb = CppAD::CondExpLt(z, Type(-30), Type(-30), z);
  Type asym = -0.5*za*za - 0.5*log(2*M_PI) - log(-za) + log(1 - 1/(za*za) + 3/(za*za*za*za));
  return CppAD::CondExpLt(z, Type(-30), asym, log(pnorm(zb)));
}

template<class Type>
Type cloglogBt(Type th, Type ta, Type t) {
  return exp(-t*th + 0.5*t*t*ta*ta + pnormlogup(Type(th/ta - t*ta)));
}

template<class Type>
Type cloglogshape3q(Type ma, Type s2, Type vs, Type ld, Type Ee) {
  Type th = log(Type(3)) - Type(1.0/6.0);
  Type t2 = Type(1.0/3.0);
  Type kap = Type(1);
  for(int it=0; it<2; it++){
    Type tht = th - kap*ma;
    Type tat = sqrt(t2 + kap*kap*s2);
    Type B1 = cloglogBt(tht, tat, Type(1));
    Type B2 = cloglogBt(tht, tat, Type(2));
    Type B3 = cloglogBt(tht, tat, Type(3));
    Type Ew = 2 + pnorm(Type(-tht/tat)) + 0.5*B1 + 0.25*B2 + 0.25*B3;
    t2 = 1/(2 + pnorm(Type(-tht/tat)) + B1 + 0.75*B2 + B3);
    kap = t2*Ew*(vs + Type(1e-12))/(s2 + Type(1e-12) - t2*Ew*(s2 - vs));
    th = log(Ew) - 0.5*t2 - (1 - kap)*ma + 0.5*ld - 0.5*(1 - kap)*(1 - kap)*vs;
  }
  Type tht = th - kap*ma;
  Type tat = sqrt(t2 + kap*kap*s2);
  Type P0 = pnorm(Type(-tht/tat));
  Type out = 3*Ee + 2*tht + tht*P0 - tat*dnorm(Type(tht/tat), Type(0), Type(1), false) - P0*log(Type(2)) + 0.5*log(2*M_PI*t2) + 0.5;
  out -= exp(th + 0.5*t2 + (1 - kap)*ma - 0.5*ld + 0.5*(1 - kap)*(1 - kap)*vs);
  Type lz = pnormlogup(Type(tht/tat));
  Type sk = 0;
  Type e = Type(1);
  for(int k=0; k<7; k++){
    sk += log(1 + exp(-e*log(Type(2)) - e*tht + 0.5*e*e*tat*tat + pnormlogup(Type(tht/tat - e*tat)) - lz));
    e = 2*e;
  }
  out -= pnorm(Type(tht/tat))*sk;
  return out;
}

template<class Type>
Type cloglogshape3(Type m, Type v) {
  return cloglogshape3q(m, v, v, Type(0), m);
}

}


// needed to use ldlt see tmb issue #398
// been addressed in TMB
// namespace std {
// template<>
// struct
//   numeric_limits<TMBad::ad_aug> : numeric_limits<double> {};
// }