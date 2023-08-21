// some distributions
namespace gllvm{

/* , int give_log=0 */
// template<class Type>
// Type dnegbinva(Type y, Type eta, Type phi, Type cQ)
// {
//   Type logres = y*(eta-cQ) - (y+phi)*log(phi+exp(eta-cQ)) + lgamma(y+phi) - phi*cQ + phi*log(phi) - lgamma(phi) -lfactorial(y);
//   // if(!give_log) return exp(logres);
//   return logres;
// }
// 
// /* Hess for nb EVA */
// template<class Type>
// Type nb_Hess(Type y, Type eta, Type phi)
// {
//   Type res = (((phi+y) / (phi+exp(eta))) * exp(eta) - ((phi+y)*pow(phi+exp(eta),-2))*pow(exp(eta),2));
//   // if(!give_log) return exp(logres);
//   return res;
// }
// 
// /* Hess for bin logit EVA */
// template<class Type>
// Type dbinom_logit_eva(Type y, Type eta, Type cQ)
// {
//   Type mu = 0.0;
//   Type mu_prime = 0.0;
//   
//   CppAD::vector<Type> z(4);
//   z[0] = eta;
//   z[1] = 0;
//   z[2] = 1/(1+exp(-z[0]));
//   z[3] = exp(z[0])/(exp(z[0])+1);
//   
//   mu = Type(CppAD::CondExpGe(z[0], z[1], z[2], z[3]));
//   mu_prime = mu * (1-mu);
//   Type logres = y * eta + log(1-mu) -  mu_prime*cQ;
//   return logres;
// }

/* Hess for bin probit EVA */
// template<class Type>
// Type bin_probit_Hess(Type y, Type eta, Type etaP, Type etaD)
// {
//   // Type etaD =  dnorm(Type(eta), Type(0), Type(1), true); // log normal density evaluated at eta(i,j)
//   Type res = ((y*(etaP*exp(etaD)*(-eta)-pow(exp(etaD),2))*pow(1-etaP,2) + (1-y)*((1-etaP)*exp(etaD)*eta-pow(exp(etaD),2))*pow(etaP,2) )/(etaP*etaP*(etaP*etaP-2*etaP+1))); 
//   return res;
// }

// AR(1) correlation matrix
template <class Type>
matrix<Type> corAR1(Type s0, Type s1, int nr)
  {
  matrix<Type> S(nr,nr);
  Type rho = s1 / sqrt(1.0 + pow(s1, 2));
  for (int d=0;d<nr;d++) {
    S(d,d)=s0*s0;
    for (int j=0;j<d;j++){
      S(d,j)=s0*pow(rho,(d-j))*s0;       // ar1 correlation
      S(j,d)=S(d,j);
    }
  }
  return S;
}

// compound symmetry correlation matrix
template <class Type>
matrix<Type> corCS(Type s0, Type s1, int nr)
{
  matrix<Type> S(nr,nr);
  Type rho = s1 / sqrt(1.0 + pow(s1, 2));
  S.fill(s0*rho*s0);
  S.diagonal().fill(s0*s0);
  // for (int d=0;d<nr;d++) {
    // S(d,d)=s0*s0;
    // for (int j=0;j<d;j++){
    //   S(d,j)=s0*rho*s0;
    //   S(j,d)=S(d,j);
    // }
  // }
  return S;
}

// Exp decaying correlation matrix
template <class Type>
matrix<Type> corExp(Type s0, Type s1, int nr, matrix<Type> dc)
{  //matrix<Type> corExp(Type s0, CppAD::vector<Type> s1, int nr, matrix<Type> dc)
  matrix<Type> S(nr,nr);
  // matrix<Type> alf(dc.cols(),dc.cols());
  // for(int i=0; i<dc.cols(); i++){//old
  //     alf(i,i) = 1/(exp(s1)*exp(s1));
  // }
  // matrix<Type> alf = 1/exp(s1);
  Type alf = 1/exp(s1);
  
  
  for (int d=0;d<nr;d++) {
    S(d,d)=s0*s0;
    for (int j=0;j<d;j++){
      S(d,j)=s0*exp(-sqrt( (((dc.row(d)-dc.row(j))*alf)*(dc.row(d)-dc.row(j)).transpose()).sum() ) )*s0;//
      // // S(d,j)=s0*exp(-sqrt(((dc.row(d)-dc.row(j))*(dc.row(d)-dc.row(j)).transpose()).sum())/alf)*s0;
      // S(d,j)=s0*exp(-dc(d,j)*alf )*s0;
      S(j,d)=S(d,j);
    }
  }
  return S;
}


// Matern correlation matrix
template <class Type>
matrix<Type> corMatern(Type s0, Type phi, Type kappa, int nr, matrix<Type> dc)
{
//   // s0 covariance
//   //phi range
//   //kappa smoothness
  matrix<Type> S(nr,nr);
  for (int d=0;d<nr;d++) {
    S(d,d)=s0*s0;
    for (int j=0;j<d;j++){
      // S(d,j)=s0*matern(dc(d,j), ph, kappa)*s0;
      S(d,j)=s0*matern(sqrt(((dc.row(d)-dc.row(j))*(dc.row(d)-dc.row(j)).transpose()).sum()), phi, kappa)*s0; //old
      S(j,d)=S(d,j);
    }
  }
  return S;
}

}


