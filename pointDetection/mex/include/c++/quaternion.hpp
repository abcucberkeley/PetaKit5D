#ifndef QUATERNION_HPP
# define QUATERNION_HPP

# include <vector.hpp>
# include <matrix.hpp>

class quaternion : public vector<4,double>
{
public:
  typedef vector<4, double> super_type;

public:
  quaternion() {}
	
  quaternion(double s, double x, double y, double z)
  {
    data_[0] = s;
    data_[1] = x;
    data_[2] = y;
    data_[3] = z;
  }
	
  quaternion(double s, const vector<3,double> & v)
  {
    data_[0] = s;
    data_[1] = v[0];
    data_[2] = v[1];
    data_[3] = v[2];
  }
	
  quaternion(const vector<4,double> & v) : super_type(v)
  {
  }
  
  quaternion & operator=(const vector<4,double> & v)
  {
    static_cast<super_type*>(this)->operator=(v);

    return *this;
  }

  const vector<4,double> & to_vector() const
  {
    return *this;
  }
	
  matrix<3,3,double> to_matrix() const
  {
    matrix<3,3,double> mat;
		
    double q0 = data_[0];
    double q1 = data_[1];
    double q2 = data_[2];
    double q3 = data_[3];
		
    mat(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    mat(1, 0) = 2 * (q1 * q2 + q0 * q3);
    mat(2, 0) = 2 * (q1 * q3 - q0 * q2);
    mat(0, 1) = 2 * (q1 * q2 - q0 * q3);
    mat(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    mat(2, 1) = 2 * (q2 * q3 + q0 * q1);
    mat(0, 2) = 2 * (q1 * q3 + q0 * q2);
    mat(1, 2) = 2 * (q2 * q3 - q0 * q1);
    mat(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
		
    return mat;
  }
	
  const double & s() const { return data_[0]; }
  double & s() { return data_[0]; }
		
  const vector<3,double> & v() const
  {
    return *(const vector<3,double>*)(const void*)(& this->data_[1]);
  }
	
  vector<3,double> & v()
  {
    return *(vector<3,double>*)(void*)(& this->data_[1]);
  }
	
  bool is_pure() const
  {
    return fabs(data_[0]) < std::numeric_limits<double>::epsilon();
  }
	
  quaternion conj() const
  {
    return quaternion(data_[0], -v());
  }
	
  quaternion inverse() const
  {
    double f = this->norm_l2();
    return conj() / (f * f);
  }
	
  vector<3,double> rotate(const vector<3,double> & v);
};

inline
quaternion operator*(const quaternion & lhs, const quaternion & rhs)
{
  return quaternion(lhs.s() * rhs.s() - lhs.v() * rhs.v(),
		    vprod(lhs.v(), rhs.v()) + lhs.v() * rhs.s() + 
		    rhs.v() * lhs.s());
}

inline
vector<3,double> quaternion::rotate(const vector<3,double> & v)
{
  return ((*this) * quaternion(0, v) * (*this).inverse()).v();
}

#endif

