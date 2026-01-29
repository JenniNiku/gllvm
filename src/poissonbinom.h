//-----------------------------------------------------------------------------
// Adapted by Bert van der Veen, with optimisation by chatGPT
// code originally from 'poisbinom.R' in the recordTest package by Jorge Castillo-Mateo
//-----------------------------------------------------------------------------

#include <unsupported/Eigen/FFT> // Use Eigen's native FFT
#include <complex>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Eigen;
using std::complex;

Eigen::MatrixXd dpoisbinom(Eigen::MatrixXd prob) {
  int m = prob.rows(); // number of distributions
  int n = prob.cols(); // number of Bernoulli trials
  int N = n + 1;       // FFT size
  
  // Handle single-trial case
  if (n == 1) {
    Eigen::MatrixXd res(m, 2);
    res.col(0) = Eigen::VectorXd::Ones(m) - prob.col(0);
    res.col(1) = prob.col(0);
    return res;
  }
  
  // FFT parameters
  int L = N / 2; 
  double w_base = 2.0 * M_PI / N;
  
  std::vector<complex<double>> factors(L);
  for(int i = 0; i < L; i++){
    double angle = w_base * (i + 1);
    factors[i] = complex<double>(cos(angle) - 1.0, sin(angle));
  }
  
  Eigen::MatrixXcd fft_data(m, N);
  
  fft_data.col(0) = Eigen::VectorXcd::Ones(m);
  
#pragma omp parallel for schedule(static)
  for (int l = 0; l < L; ++l) {
    
    std::vector<double> real_acc(m, 0.0);
    std::vector<double> imag_acc(m, 0.0);
    
    complex<double> f = factors[l];
    double f_re = f.real();
    double f_im = f.imag();
    
    for(int t = 0; t < n; ++t) {
      
      for(int d = 0; d < m; ++d) {
        double p = prob(d, t);
        
        double z_re = 1.0 + p * f_re;
        double z_im = p * f_im;
        
        double norm_sq = z_re*z_re + z_im*z_im;
        real_acc[d] += 0.5 * std::log(norm_sq);
        imag_acc[d] += std::atan2(z_im, z_re);
      }
    }
    
    for(int d = 0; d < m; ++d) {
      double mag = std::exp(real_acc[d]);
      double theta = imag_acc[d];
      
      complex<double> val(mag * std::cos(theta), mag * std::sin(theta));
      
      int idx_left = l + 1;
      int idx_right = N - 1 - l;
      
      fft_data(d, idx_left) = val;
      if(idx_left != idx_right) {
        fft_data(d, idx_right) = std::conj(val);
      }
    }
  }
  
  Eigen::MatrixXd result(m, N);
  
#pragma omp parallel
{
  Eigen::FFT<double> fft;
  
#pragma omp for
  for(int d = 0; d < m; ++d) {

    Eigen::VectorXcd row = fft_data.row(d) / (double)N;
    Eigen::VectorXcd out_freq(N);
    
    fft.fwd(out_freq, row); 
    
    for(int i = 0; i < N; ++i) {
      result(d, i) = std::abs(out_freq[i]); 
    }
  }
}

return result;
}
