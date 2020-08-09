#ifndef	NON_SEPARABLE_FILTER_2D_HPP
# define NON_SEPARABLE_FILTER_2D_HPP

# include <cassert>

class non_separable_filter_2d
{
public:
  non_separable_filter_2d(int width, int height) :
    width_(width),
    height_(height)
  {
    data_ = new double*[width_];

    for (int x = 0; x < width_; ++x)
      data_[x] = new double[height_];
  }
  
  ~non_separable_filter_2d()
  {
    for (int x = 0; x < width_; ++x)
      delete[] data_[x];
    delete[] data_;
  }

public:
  int width() const { return width_; }
  int height() const { return height_; }

  double operator()(int x, int y) const
  {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);

    return data_[x][y];
  }

  double & operator()(int x, int y)
  {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);

    return data_[x][y];
  }

  double norm_l1() const
  {
    double res = 0;

    for (int x = 0; x < width_; ++x)
      for (int y = 0; y < height_; ++y)
	res += data_[x][y];

    return res;
  }

  double norm_l2() const
  {
    double res = 0;

    for (int x = 0; x < width_; ++x)
      for (int y = 0; y < height_; ++y)
	res += data_[x][y] * data_[x][y];

    return sqrt(res);
  }

private:
  int width_, height_;

  double** data_;
};

#endif /* !NON_SEPARABLE_FILTER_2D_HPP */
