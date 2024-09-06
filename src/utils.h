// some functions
namespace gllvmutils{

template<class Type>
matrix<Type> cov2cor(matrix<Type> x){
  return x.diagonal().cwiseInverse().cwiseSqrt().asDiagonal()*x*x.diagonal().cwiseInverse().cwiseSqrt().asDiagonal();
}

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
namespace std {
template<>
struct
  numeric_limits<TMBad::ad_aug> : numeric_limits<double> {};
}