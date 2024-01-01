#ifndef   	CONVOLUTION_HPP
# define   	CONVOLUTION_HPP

# include <vector>
# include <cassert>

# include <gaussian_derivative_1d.hpp>
# include <non_separable_filter_2d.hpp>
# include <image.hpp>

//////////////////////////////
// 2D separable convolution //
//////////////////////////////

template <typename T, int N, int M>
inline void convolve(const image<2, T> & src,
		     const gaussian_derivative_1d<N> & win1,
		     const gaussian_derivative_1d<M> & win2,
		     image<2, double> & dst,
		     bool mirror = false)
{
  assert(dst.width() == src.width() && dst.height() == src.height());

  int hside1 = win1.size() >> 1;
  int hside2 = win2.size() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.border_mirror(new_margin);
  else
    src.border_replicate(new_margin);

  image<2, double> tmp(src.size(), src.margin());

  double sum;

  // for all rows
  for (int j = 0; j < src.height(); ++j)
    for (int i = 0; i < src.width(); ++i)
      {
	sum = 0.0;
	for (int k = -hside1; k <= hside1; ++k)
	  sum += win1[k + hside1] * src(i - k, j);
	tmp(i, j) = sum;
      }

  if (mirror)
    tmp.border_mirror(new_margin);
  else
    tmp.border_replicate(new_margin);

  // for all columns
  for (int i = 0; i < src.width(); ++i)
    for (int j = 0; j < src.height(); ++j)
      {
	sum = 0.0;
	for (int k = -hside2; k <= hside2; ++k)
	  sum += win2[k + hside2] * tmp(i, j - k);
	dst(i, j) = sum;
      }
}

//////////////////////////////
// 3D separable convolution //
//////////////////////////////

template <typename T, int N, int M, int L>
inline void convolve(const image<3, T> & src,
		     const gaussian_derivative_1d<N> & win1,
		     const gaussian_derivative_1d<M> & win2,
		     const gaussian_derivative_1d<L> & win3,
		     image<3, double> & dst,
		     bool mirror = false)
{
  assert(dst.width() == src.width() &&
	 dst.height() == src.height() &&
	 dst.depth() == src.depth());

  int hside1 = win1.size() >> 1;
  int hside2 = win2.size() >> 1;
  int hside3 = win3.size() >> 1;

  int new_margin = std::max(std::max(hside1, hside2), hside3);

  if (mirror)
    src.border_mirror(new_margin);
  else
    src.border_replicate(new_margin);

  std::vector<double> tmp_1d(src.depth());
  image<2, double> tmp_2d(src.size());

  double sum;

  for (int k = 0; k < src.depth(); ++k)
    {
      // for all rows
      for (int j = 0; j < src.height(); ++j)
	for (int i = 0; i < src.width(); ++i)
	  {
	    sum = 0.0;
	    for (int l = -hside1; l <= hside1; ++l)
	      sum += win1[l + hside1] * src(i - l, j, k);
	    tmp_2d(i, j) = sum;
	  }
      
      if (mirror)
	tmp_2d.border_mirror(new_margin);
      else
	tmp_2d.border_replicate(new_margin);
      
      // for all columns
      for (int i = 0; i < src.width(); ++i)
	for (int j = 0; j < src.height(); ++j)
	  {
	    sum = 0.0;
	    for (int l = -hside2; l <= hside2; ++l)
	      sum += win2[l + hside2] * tmp_2d(i, j - l);
	    dst(i, j, k) = sum;
	  }
    }

  if (mirror)
    dst.border_mirror(new_margin);
  else
    dst.border_replicate(new_margin);

  // for all z
  for (int j = 0; j < src.height(); ++j)
    for (int i = 0; i < src.width(); ++i)
      {
	for (int k = 0; k < src.depth(); ++k)
	  {
	    sum = 0.0;
	    for (int l = -hside3; l <= hside3; ++l)
	      sum += win3[l + hside3] * dst(i, j, k - l);
	    tmp_1d[k] = sum;
	  }

	for (int k = 0; k < src.depth(); ++k)
	  dst(i, j, k) = tmp_1d[k];
      }
}

//////////////////////////////////
// 2D non separable convolution //
//////////////////////////////////

template <typename T>
inline void convolve(const image<2, T> & src,
		     const non_separable_filter_2d & win,
		     image<2, double> & dst,
		     bool mirror = false)
{
  assert(dst.width() == src.width() && dst.height() == src.height());

  int hside1 = win.width() >> 1;
  int hside2 = win.height() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.border_mirror(new_margin);
  else
    src.border_replicate(new_margin);

  double sum;
  
  for (int i = 0; i < src.width(); ++i)
    for (int j = 0; j < src.height(); ++j)
      {
	sum = 0.0;
	for (int x = -hside1; x <= hside1; ++x)
	  for (int y = -hside2; y <= hside2; ++y)
	    sum += win(x + hside1, y + hside2) * src(i - x, j - y);
	dst(i, j) = sum;
      }  
}

#endif	    /* !CONVOLUTION_HPP */
