// some functions
namespace gllvmutils{
// 
// //numerically safe coshl with known derivative tanh
// // Double version (for ATOMIC_DOUBLE)
// // assumes x is non-negative
// double logcosh(double x) {
//   //log(cosh(sqrt(x^2+a))) is a pain numerically
//   //I don't think anything can be done for the function inside the cosh
//   //but this version of log(cosh(x)) performs best in my tests
//   //in comparison to e.g., abs(x)+log1p(exp(-2*abs(x)))-log(2.0);
//   //or something like if(|x|<1e-4){0.5x^2}else if(|x|<17){log(cosh(x))}else{|x|-log(2)}
//   double x_abs = ((x>0)-(x<0))*x;
//   if(x_abs<1e-4){
//     return 0.5*x*x;
//   }else if(x_abs<10){//not technically wrong, but this should be x_abs..
//     return x_abs + log1p(exp(-2*x_abs))-log(2.0);//std::log(std::coshl(x));
//   }else{
//     return x_abs- log(2.0);
//   }
// }
// 
// template<class Type>
// Type dlogcosh(Type x){
//   //log(cosh(sqrt(x^2+a))) is a pain numerically
//   //I don't think anything can be done for the function inside the cosh
//   //but this version of log(cosh(x)) performs best in my tests
//   //in comparison to e.g., abs(x)+log1p(exp(-2*abs(x)))-log(2.0);
//   //or something like if(|x|<1e-4){0.5x^2}else if(|x|<17){log(cosh(x))}else{|x|-log(2)}
//   double x_abs = ((x>0)-(x<0))*x;
//   Type y = 0;
//   y = CppAD::CondExpLt(x_abs, 1e-4, x, CppAD::CondExpGt(x_abs,10,((x>0)-(x<0))*1, tanh(x)));
// 
//   return y;
// }
// 
// TMB_ATOMIC_VECTOR_FUNCTION(
//   // ATOMIC_NAME
//   coshl
//   ,
//   // OUTPUT_DIM
//   1,
//   // ATOMIC_DOUBLE
//   ty[0] = logcosh(tx[0]);
// ,
// // ATOMIC_REVERSE
// Type d_result = dlogcosh(tx[0]);//dlogcosh(tx[0]);//tanh(0.5*sqrt(tx[0]))*1/(4*sqrt(tx[0]));
// px[0] = d_result * py[0];  // Reverse mode chain rule
// )
// 
// // Scalar version
// template<class Type>
// Type logcosh(Type x){
//   CppAD::vector<Type> tx(1);
//   tx[0] = x;
//   return coshl(tx)[0];
// }  

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
  
//structure to pass list from R as list in TMB
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

}


// needed to use ldlt see tmb issue #398
// been addressed in TMB
// namespace std {
// template<>
// struct
//   numeric_limits<TMBad::ad_aug> : numeric_limits<double> {};
// }