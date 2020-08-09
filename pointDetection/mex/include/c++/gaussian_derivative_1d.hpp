#ifndef	GAUSSIAN_DERIVATIVE_1D_HPP
# define GAUSSIAN_DERIVATIVE_1D_HPP

# include <cassert>

template <int N>
struct hermite_polynomial
{
  static double res(double x, double inv_v)
  {
    return inv_v * (x * hermite_polynomial<N - 1>::res(x, inv_v) -
		    hermite_polynomial<N - 1>::derivate(x, inv_v));
  }

  static double derivate(double x, double inv_v)
  {
    return N * hermite_polynomial<N - 1>::res(x, inv_v);
  }
};

template <>
struct hermite_polynomial<0>
{
  static double res(double x, double inv_v) { return 1; }

  static double derivate(double x, double inv_v) { return 0; }
};

template <int N>
class gaussian_derivative_1d
{
public:
  gaussian_derivative_1d(double sigma) :
    sigma_(sigma),
    size_(0),
    data_(0)
  {
    // Note: sigma should be > 0. However, sigma < 1 yields to a very
    // degrated discretization of the Gaussian filter.

    assert(sigma >= 1);

    int hside = ((int) ceil(6 * sigma));

    size_ = 2 * hside + 1;

    data_ = new double[size_];

    double inv_v = 1 / (sigma * sigma);

    double c = (N & 1 ? -1 : 1) / (sqrt(2 * M_PI) * sigma);
    
    for (int x = -hside; x <= hside; ++x)
      data_[x + hside] = c * hermite_polynomial<N>::res(x, inv_v) *
	exp(- 0.5 * x * x * inv_v);
  }
  
  ~gaussian_derivative_1d()
  {
    delete[] data_;
  }

public:
  double sigma() const { return sigma_; }

  int size() const { return size_; }

  double operator[](int x) const
  {
    assert(x >= 0 && x < size_);

    return data_[x];
  }

  double norm_l1() const
  {
    double res = 0;

    for (int x = 0; x < size_; ++x)
      res += data_[x];

    return res;
  }

  double norm_l2() const
  {
    double res = 0;

    for (int x = 0; x < size_; ++x)
      res += data_[x] * data_[x];

    return sqrt(res);
  }

private:
  double sigma_;

  int size_;

  double* data_;
};

typedef gaussian_derivative_1d<0> gaussian;

#endif /* !GAUSSIAN_DERIVATIVE_1D_HPP */
