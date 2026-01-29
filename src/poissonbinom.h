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
  int L = N / 2; // integer division (floor) is fine here given the logic
  double w_base = 2.0 * M_PI / N;
  
  // Precompute Twiddle Factors (e^{iw} - 1)
  // This vector is small, so we just compute it once.
  std::vector<complex<double>> factors(L);
  for(int i = 0; i < L; i++){
    double angle = w_base * (i + 1);
    factors[i] = complex<double>(cos(angle) - 1.0, sin(angle));
  }
  
  // We will construct the frequency domain signal directly.
  // We need to store complex values for the FFT.
  // Dimensions: m rows (distributions), N cols (frequencies)
  // We store row-wise for distributions to make FFT loop easy later.
  Eigen::MatrixXcd fft_data(m, N);
  
  // set DC component (k=0) to 1.0 (prob sum is 1)
  fft_data.col(0) = Eigen::VectorXcd::Ones(m);
  
  // 1. Calculate Characteristic Function (Log Domain)
  // Parallelize over frequencies (L)
#pragma omp parallel for schedule(static)
  for (int l = 0; l < L; ++l) {
    
    // Temporary accumulators for this frequency
    // We use vectors to allow the compiler to vectorize the inner loop (SIMD)
    std::vector<double> real_acc(m, 0.0);
    std::vector<double> imag_acc(m, 0.0);
    
    complex<double> f = factors[l];
    double f_re = f.real();
    double f_im = f.imag();
    
    // Loop over trials
    for(int t = 0; t < n; ++t) {
      
      // INNER LOOP: Distributions
      // Iterating 'd' here is cache-friendly (contiguous memory in prob)
      // This is the most critical optimization.
      for(int d = 0; d < m; ++d) {
        double p = prob(d, t);
        
        // z = 1 + p * factor
        double z_re = 1.0 + p * f_re;
        double z_im = p * f_im;
        
        // Accumulate log(z)
        // log(|z|) = 0.5 * log(re^2 + im^2) -> Avoids sqrt()
        double norm_sq = z_re*z_re + z_im*z_im;
        real_acc[d] += 0.5 * std::log(norm_sq);
        imag_acc[d] += std::atan2(z_im, z_re);
      }
    }
    
    // Write back to the FFT matrix with symmetry
    for(int d = 0; d < m; ++d) {
      // Exponentiate back to get CF value
      double mag = std::exp(real_acc[d]);
      double theta = imag_acc[d];
      
      complex<double> val(mag * std::cos(theta), mag * std::sin(theta));
      
      // Fill symmetric parts of the spectrum
      // P_l and P_{N-l} are conjugates for real output
      int idx_left = l + 1;
      int idx_right = N - 1 - l;
      
      fft_data(d, idx_left) = val;
      if(idx_left != idx_right) {
        fft_data(d, idx_right) = std::conj(val);
      }
    }
  }
  
  // 2. Perform FFT
  // Parallelize over distributions
  Eigen::MatrixXd result(m, N);
  
#pragma omp parallel
{
  // Create a thread-local FFT plan (setup can be expensive)
  Eigen::FFT<double> fft;
  
#pragma omp for
  for(int d = 0; d < m; ++d) {
    // Extract row, FFT, and store magnitude
    // R's fft() computes unnormalized DFT.
    // Poisson-Binomial formula usually requires IDFT scaling (1/N).
    // However, the original code divided input by (n+1) and used R's fft (Forward).
    // We replicate the exact logic: Input / N -> Forward FFT -> Abs.
    
    Eigen::VectorXcd row = fft_data.row(d) / (double)N;
    Eigen::VectorXcd out_freq(N);
    
    // Eigen's fwd matches R's fft behavior (Forward DFT)
    fft.fwd(out_freq, row); 
    
    // Store real part (magnitude of real signal is abs(real))
    // Theoretically imaginary part is 0 due to symmetry.
    for(int i = 0; i < N; ++i) {
      result(d, i) = std::abs(out_freq[i]); 
    }
  }
}

return result;
}
