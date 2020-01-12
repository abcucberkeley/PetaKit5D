#ifndef   	SOLVER_HPP
# define   	SOLVER_HPP

//# include <cmath>
# include <math.h>
# include <cassert>
# include <limits>

class solver
{
public:
  double prec() const { return prec_; }

  int nroots() const { return nroots_; }

  double root(int k) const { return roots_[k]; }

public:
  solver(double prec = std::numeric_limits<double>::epsilon()) :
    prec_(prec), nroots_(0)
  {
    for (int i = 0; i < 4; ++i)
      roots_[i] = 0;
  }

  // Resolve a linear equation of form: aX + b = 0 in R
  void operator()(double a, double b)
  {
    nroots_ = 0;

    if (a != 0)
      {
	roots_[0] = -b / a;
	nroots_ = 1;
      }
  }

  // Resolve a quadratic equation of form: aX2 + bX + c = 0 in R
  // ref: "Numerical Recipes in C" p. 183--184
  void operator()(double a, double b, double c)
  {
    nroots_ = 0;

    if (a == 0)
      {
	this->operator()(b, c);
	return;
      }

    double delta = b * b - 4 * a * c;

    if (fabs(delta) < prec_)
      {
	roots_[0] = - b / (2 * a);
	nroots_ = 1;
      }
    else
      if (delta > 0)
	{
	  double q = - 0.5 * (b + (b < 0 ? -1 : 1) * sqrt(delta));

	  roots_[0] = q / a;
	  roots_[1] = c / q;
	  nroots_ = 2;
	}
  }

  // Resolve a cubic equation of form:
  // aX3 + bX2 + cX + d = 0 in R
  void operator()(double a, double b, double c, double d)
  {
    nroots_ = 0;

    if (a == 0)
      {
	this->operator()(b, c, d);
	return;
      }

    cubic_(b / a, c / a, d / a);
  }

  // Resolve a 4th order equation of form: 
  // aX4 + bX3 + cX2 + dX + e = 0 in R
  void operator()(double a, double b, double c, double d, double e)
  {
    nroots_ = 0;

    if (a == 0)
      {
	this->operator()(b, c, d, e);
	return;
      }

    quartic_(b / a, c / a, d / a, e / a);
  }

private:
  // ref: "Numerical Recipes in C" p. 184--185
  void cubic_(double a, double b, double c)
  {
    double q = (a * a - 3 * b) / 9;
    double r = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

    double r2 = r * r;
    double q3 = q * q * q;
    double a_3 = a / 3;

    if (r2 < q3)
      {
	double theta = acos(r / sqrt(q3));
	  
	roots_[0] = - 2 * sqrt(q) * cos(theta / 3) - a_3;
	roots_[1] = - 2 * sqrt(q) * cos((theta + 2 * M_PI) / 3) - a_3;
	roots_[2] = - 2 * sqrt(q) * cos((theta - 2 * M_PI) / 3) - a_3;

	nroots_ = 3;
      }
    else
      {
	double a = - (r < 0 ? -1 : 1) * pow(fabs(r) + sqrt(r2 - q3), 1 / 3.0);

	double b = fabs(a) < prec_ ? 0.0 : q / a;

	roots_[0] = a + b - a_3;

	nroots_ = 1;
      }
  }

  // ref: http://mathworld.wolfram.com/QuarticEquation.html
  void quartic_(double a3, double a2, double a1, double a0)
  {
    cubic_(-a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);

    double y1 = roots_[0];

    nroots_ = 0;

    double delta_r = 0.25 * a3 * a3 - a2 + y1;

    if (delta_r < 0)
      return;

    double r = sqrt(delta_r);

    double d, e;
    double delta_d, delta_e;

    if (fabs(r) < prec_)
      {
	delta_d = 0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0);
	delta_e = 0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0);

	if (delta_d >= 0)
	  {
	    d = sqrt(delta_d);
	    roots_[0] = - 0.25 * a3 + 0.5 * d;
	    roots_[1] = - 0.25 * a3 - 0.5 * d;
	    nroots_ = 2;
	  }

	if (delta_e >= 0)
	  {
	    e = sqrt(delta_e);
	    roots_[nroots_] = - 0.25 * a3 + 0.5 * e;
	    roots_[nroots_ + 1] = - 0.25 * a3 - 0.5 * e;
	    nroots_ += 2;
	  }
      }
    else
      {
	delta_d = 0.75 * a3 * a3 - r * r - 2 * a2 +
	  (a3 * a2 - 2 * a1 - 0.25 * a3 * a3 * a3) / r;

	delta_e = 0.75 * a3 * a3 - r * r - 2 * a2 -
	  (a3 * a2 - 2 * a1 - 0.25 * a3 * a3 * a3) / r;

	if (delta_d >= 0)
	  {
	    d = sqrt(delta_d);
	    roots_[0] = - 0.25 * a3 + 0.5 * r + 0.5 * d;
	    roots_[1] = - 0.25 * a3 + 0.5 * r - 0.5 * d;
	    nroots_ = 2;
	  }

	if (delta_e >= 0)
	  {
	    e = sqrt(delta_e);
	    roots_[nroots_] = - 0.25 * a3 - 0.5 * r + 0.5 * e;
	    roots_[nroots_ + 1] = - 0.25 * a3 - 0.5 * r - 0.5 * e;
	    nroots_ += 2;
	  }
      }
  }

private:
  double prec_;

  double roots_[4];
  int nroots_;
};

#endif	    /* !SOLVER_HPP */
