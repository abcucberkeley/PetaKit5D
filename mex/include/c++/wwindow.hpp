#ifndef WWINDOW_HPP
# define WWINDOW_HPP

# include <window.hpp>

template <unsigned N, unsigned S, typename T, typename W>
class wwindow : public window<N,S,T>
{
public:
  wwindow() {}

  wwindow(const vector<N,T> points[], const W weights[]) :
    window<N,S,T>(points)
  {
    std::copy(weights, weights + S, weights_);
  }
  
  const W & weight(unsigned i) const { return weights_[i]; }

  W & weight(unsigned i) { return weights_[i]; }

private:
  W weights_[S];
};

template <typename T, typename W>
inline
const wwindow<2,2,T,W> & chanfrein11_fwd()
{
  static wwindow<2,2,T,W> win;
  static bool first = true;
  
  if (first)
    {
      win.point(0)[0] = -1;	win.point(0)[1] = 0;	win.weight(0) = 1;
      win.point(1)[0] = 0;	win.point(1)[1] = -1;	win.weight(1) = 1;
      first = false;
    }
  return win;
}

template <typename T, typename W>
inline
const wwindow<2,2,T,W> & chanfrein11_bkd()
{
  static wwindow<2,2,T,W> win;
  static bool first = true;
  
  if (first)
    {
      win.point(0)[0] = 1;	win.point(0)[1] = 0;	win.weight(0) = 1;
      win.point(1)[0] = 0;	win.point(1)[1] = 1;	win.weight(1) = 1;
      first = false;
    }
  return win;
}

template <typename T, typename W>
inline
const wwindow<3,3,T,W> & chanfrein111_fwd()
{
	static wwindow<3,3,T,W> win;
	static bool first = true;
	
	if (first)
    {
		win.point(0)[0] = -1;	win.point(0)[1] = 0;	win.point(0)[2] = 0;    win.weight(0) = 1;
		win.point(1)[0] = 0;	win.point(1)[1] = -1;	win.point(1)[2] = 0;    win.weight(1) = 1;
		win.point(2)[0] = 0;	win.point(2)[1] = 0;	win.point(2)[2] = -1;    win.weight(2) = 1;
		first = false;
    }
	return win;
}

template <typename T, typename W>
inline
const wwindow<3,3,T,W> & chanfrein111_bkd()
{
	static wwindow<3,3,T,W> win;
	static bool first = true;
	
	if (first)
    {
		win.point(0)[0] = 1;	win.point(0)[1] = 0;	win.point(0)[2] = 0;    win.weight(0) = 1;
		win.point(1)[0] = 0;	win.point(1)[1] = 1;	win.point(1)[2] = 0;    win.weight(1) = 1;
		win.point(2)[0] = 0;	win.point(2)[1] = 0;	win.point(2)[2] = 1;    win.weight(2) = 1;
		first = false;
    }
	return win;
}

#endif /* WWINDOW_HPP */
