#ifndef VECTOR_HPP
# define VECTOR_HPP

# include <limits>
# include <cmath>
# include <cassert>
# include <iostream>

template <unsigned n, typename T>
class vector
{
public:
  vector()
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = 0;
  }
	
  template <typename U>
  vector(const vector<n,U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = rhs[i];
  }
	
  template <typename U>
  vector & operator=(const vector<n,U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = rhs[i];
    return *this;
  }
	
  template <typename U>
  vector & set_all(U rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = rhs;
    return *this;
  }

  const T & operator[](unsigned i) const
  {
    return data_[i];
  }

  T & operator[](unsigned i)
  {
    return data_[i];
  }
	
  vector operator-() const
  {
    vector tmp;
    for (unsigned i = 0; i < n; ++i)
      tmp[i] = -data_[i];
    return tmp;
  }

  vector & operator+=(const vector & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] += rhs[i];
    return *this;
  }
	
  vector & operator/=(T val)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] /= val;
    return *this;
  }
	
  unsigned size() const
  {
    return n;
  }
	
  bool is_unit() const
  {
    return fabs(this->norm() - 1.0) < std::numeric_limits<double>::epsilon();
  }
	
  bool is_null() const
  {
    return this->norm() < std::numeric_limits<double>::epsilon();
  }
	
  double norm_l2() const
  {
    double n_l2 = 0;
    for (unsigned i = 0; i < n; ++i)
      n_l2 += data_[i] * data_[i];
    return sqrt(n_l2);
  }
	
  const vector & normalize()
  {
    double n_l2 = this->norm_l2();
    for (unsigned i = 0; i < n; ++i)
      data_[i] = T(data_[i] / n_l2);
    return *this;
  }
	
protected:
  T data_[n];
};

template <unsigned n, typename T, typename U>
inline
bool operator==(const vector<n,T> & lhs, const vector<n,U> & rhs)
{
  for (unsigned i = 0; i < n; ++i)
    if (lhs[i] != rhs[i])
      return false;
  return true;
}

template <unsigned n, typename T>
inline
vector<n,T> operator+(const vector<n,T> & lhs, const vector<n,T> & rhs)
{
  vector<n,T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] + rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n,T> operator-(const vector<n,T> & lhs, const vector<n,T> & rhs)
{
  vector<n,T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] - rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
T operator*(const vector<n,T> & lhs, const vector<n,T> & rhs)
{
  T tmp = 0;
  for (unsigned i = 0; i < n; ++i)
    tmp += lhs[i] * rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n,T> operator*(const vector<n,T> & lhs, const T & rhs)
{
  vector<n,T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] * rhs;
  return tmp;
}

template <unsigned n, typename T, typename U>
inline
vector<n,U> operator/(const vector<n,T> & lhs, const U & rhs)
{
  vector<n,U> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] / rhs;
  return tmp;
}

template <unsigned n, typename T1, typename T2>
inline
double dist(const vector<n,T1> & lhs, const vector<n,T2> & rhs)
{
  double tmp, d = 0;
  for (unsigned i = 0; i < n; ++i)
    {
      tmp = lhs[i] - rhs[i];
      d += tmp * tmp;
    }
  return sqrt(d);
}

template <typename T>
inline
vector<3, T> vprod(const vector<3,T> & lhs, const vector<3,T> & rhs)
{
  vector<3,T> tmp;
  tmp[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
  tmp[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
  tmp[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
  return tmp;
}

template <unsigned n, typename T>
inline
std::ostream & operator<<(std::ostream & ostr, const vector<n,T> & v)
{
  ostr << '(';
  for (unsigned i = 0; i < n; ++i)
    ostr << v[i] << (i == n - 1 ? ")" : ", ");
  return ostr;
}

#endif
