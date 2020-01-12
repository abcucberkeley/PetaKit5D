#ifndef WINDOW_HPP
# define WINDOW_HPP

# include <algorithm>

# include <vector.hpp>

template <unsigned N, unsigned S, typename T>
class window
{
public:
  window() {}

  window(const vector<N,T> points[])
  {
    std::copy(points, points + S, points_);
  }
  
  const vector<N,T> & point(unsigned i) const { return points_[i]; }

  vector<N,T> & point(unsigned i) { return points_[i]; }

  static unsigned size() { return S; }
  
private:
  vector<N,T> points_[S];
};

template <typename T>
inline
const window<2,4,T> & neighb_c4()
{
  static window<2,4,T> win;
  static bool first = true;
  
  if (first)
    {
      win.point(0)[0] = 1;	win.point(0)[1] = 0;
      win.point(1)[0] = -1;	win.point(1)[1] = 0;
      win.point(2)[0] = 0;	win.point(2)[1] = 1;
      win.point(3)[0] = 0;	win.point(3)[1] = -1;
      first = false;
    }
  return win;
}

template <typename T>
inline
const window<2,8,T> & neighb_c8()
{
  static window<2,8,T> win;
  static bool first = true;

  if (first)
    {
      win.point(0)[0] = -1;	win.point(0)[1] = -1;
      win.point(1)[0] = 0;	win.point(1)[1] = -1;
      win.point(2)[0] = 1;	win.point(2)[1] = -1;
      win.point(3)[0] = -1;	win.point(3)[1] = 0;
      win.point(4)[0] = 1;	win.point(4)[1] = 0;
      win.point(5)[0] = -1;	win.point(5)[1] = 1;
      win.point(6)[0] = 0;	win.point(6)[1] = 1;
      win.point(7)[0] = 1;	win.point(7)[1] = 1;
      first = false;
    }
  
  return win;
}

template <typename T>
inline
const window<3,6,T> & neighb_c6()
{
  static window<3,6,T> win;
  static bool first = true;

  if (first)
    {
      win.point(0)[0] = -1; win.point(0)[1] = 0;  win.point(0)[2] = 0;
      win.point(1)[0] = 0;  win.point(1)[1] = -1; win.point(1)[2] = 0;
      win.point(2)[0] = 1;  win.point(2)[1] = 0;  win.point(2)[2] = 0;
      win.point(3)[0] = 0;  win.point(3)[1] = 1;  win.point(3)[2] = 0;
      win.point(4)[0] = 0;  win.point(4)[1] = 0;  win.point(4)[2] = 1;
      win.point(5)[0] = 0;  win.point(5)[1] = 0;  win.point(5)[2] = -1;
      first = false;
    }

  return win;
}

template <typename T>
inline
const window<3,26,T> & neighb_c26()
{
  static window<3,26,T> win;
  static bool first = true;

  if (first)
    {
      win.point(0)[0]  = -1; win.point(0)[1]  = 0;  win.point(0)[2] = 0;
      win.point(1)[0]  = 0;  win.point(1)[1]  = -1; win.point(1)[2] = 0;
      win.point(2)[0]  = 1;  win.point(2)[1]  = 0;  win.point(2)[2] = 0;
      win.point(3)[0]  = 0;  win.point(3)[1]  = 1;  win.point(3)[2] = 0;      
      win.point(4)[0]  = -1; win.point(4)[1]  = -1; win.point(4)[2] = 0;
      win.point(5)[0]  = -1; win.point(5)[1]  = 1;  win.point(5)[2] = 0;
      win.point(6)[0]  = 1;  win.point(6)[1]  = -1; win.point(6)[2] = 0;
      win.point(7)[0]  = 1;  win.point(7)[1]  = 1;  win.point(7)[2] = 0;

      win.point(8)[0]  = -1; win.point(8)[1]  = -1; win.point(8)[2] = 1;
      win.point(9)[0]  = 0;  win.point(9)[1]  = -1; win.point(9)[2] = 1;
      win.point(10)[0] = 1;  win.point(10)[1] = -1; win.point(10)[2] = 1;
      win.point(11)[0] = -1; win.point(11)[1] = 0;  win.point(11)[2] = 1;
      win.point(12)[0] = 1;  win.point(12)[1] = 0;  win.point(12)[2] = 1;
      win.point(13)[0] = -1; win.point(13)[1] = 1;  win.point(13)[2] = 1;
      win.point(14)[0] = 0;  win.point(14)[1] = 1;  win.point(14)[2] = 1;
      win.point(15)[0] = 1;  win.point(15)[1] = 1;  win.point(15)[2] = 1;
      win.point(16)[0] = 0;  win.point(16)[1] = 0;  win.point(16)[2] = 1;

      win.point(17)[0] = -1; win.point(17)[1] = -1; win.point(17)[2] = -1;
      win.point(18)[0] = 0;  win.point(18)[1] = -1; win.point(18)[2] = -1;
      win.point(19)[0] = 1;  win.point(19)[1] = -1; win.point(19)[2] = -1;
      win.point(20)[0] = -1; win.point(20)[1] = 0;  win.point(20)[2] = -1;
      win.point(21)[0] = 1;  win.point(21)[1] = 0;  win.point(21)[2] = -1;
      win.point(22)[0] = -1; win.point(22)[1] = 1;  win.point(22)[2] = -1;
      win.point(23)[0] = 0;  win.point(23)[1] = 1;  win.point(23)[2] = -1;
      win.point(24)[0] = 1;  win.point(24)[1] = 1;  win.point(24)[2] = -1;
      win.point(25)[0] = 0;  win.point(25)[1] = 0;  win.point(25)[2] = -1;      

      first = false;
    }

  return win;
}

#endif /* WINDOW_HPP */
