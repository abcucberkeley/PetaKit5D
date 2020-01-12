/* psfmath.h contains functions required for the calculation of point spread function models.
 *
 * Copyright (C) 2005-2013 Francois Aguet
 *
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PSFMATH_H
#define PSFMATH_H

#include <complex>
#define PI 3.14159265358979311599796346854

namespace std {
    
typedef struct {
    double ti0;
    double ni0;
    double ni0_2; // notation for Ni_0^2
    double ni;
    double ni_2;
    double tg0;
    double tg;
    double ng0;
    double ng0_2;
    double ng;
    double ng_2;
    double ns;
    double ns_2;
    double lambda;
    double k0;
    double M;
    double NA;
    double NA_2;
    double alpha;
    double pixelSize;
    int sf;
    int mode;
} parameters;


// Constants for polynomial Bessel function approximation from [Abramowitz (p. 369)]
const double j0c[7] = {1, -2.2499997, 1.2656208, -0.3163866, 0.0444479, -0.0039444, 0.0002100};
const double t0c[7] = {-.78539816, -.04166397, -.00003954, 0.00262573, -.00054125, -.00029333, .00013558};
const double f0c[7] = {.79788456, -0.00000077, -.00552740, -.00009512, 0.00137237, -0.00072805, 0.00014476};
const double j1c[7] = {0.5, -0.56249985, 0.21093573, -0.03954289, 0.00443319, -0.00031761, 0.00001109};
const double f1c[7] = {0.79788456, 0.00000156, 0.01659667, 0.00017105, -0.00249511, 0.00113653, -0.00020033};
const double t1c[7] = {-2.35619449, 0.12499612, 0.00005650, -0.00637897, 0.00074348, 0.00079824, -0.00029166};


// Bessel functions J0(x) and J1(x)
// Uses the polynomial approximations on p. 369-70 of Abramowitz & Stegun (1972).
// The error in J0 is supposed to be less than or equal to 5 x 10^-8.
__inline double J0(double x) {
	double r;

	if (x < 0.0) x *= -1.0;

	if (x <= 3.0){
		double y = x*x/9.0;
		r = j0c[0] + y*(j0c[1] + y*(j0c[2] + y*(j0c[3] + y*(j0c[4] + y*(j0c[5] + y*j0c[6])))));
	} else {
		double y = 3.0/x;
		double theta0 = x + t0c[0] + y*(t0c[1] + y*(t0c[2] + y*(t0c[3] + y*(t0c[4] + y*(t0c[5] + y*t0c[6])))));
		double f0 = f0c[0] + y*(f0c[1] + y*(f0c[2] + y*(f0c[3] + y*(f0c[4] + y*(f0c[5] + y*f0c[6])))));
		r = sqrt(1.0/x)*f0*cos(theta0);
	}
	return r;
}

__inline double J1(double x) {
	double r;
    double sign = 1.0;
	if (x < 0.0) {
        x *= -1.0;
        sign *= -1.0;
    }
	if (x <= 3.0){
		double y = x*x/9.0;
		r = x*(j1c[0] + y*(j1c[1] + y*(j1c[2] + y*(j1c[3] + y*(j1c[4] + y*(j1c[5] + y*j1c[6]))))));
	} else {
		double y = 3.0/x;
		double theta1 = x + t1c[0] + y*(t1c[1] + y*(t1c[2] + y*(t1c[3] + y*(t1c[4] + y*(t1c[5] + y*t1c[6])))));
		double f1 = f1c[0] + y*(f1c[1] + y*(f1c[2] + y*(f1c[3] + y*(f1c[4] + y*(f1c[5] + y*f1c[6])))));
		r = sqrt(1.0/x)*f1*cos(theta1);
	}
	return sign*r;
}


// Evaluates the optical path difference, with derivative d/d_theta in theta, the angle between 0 and alpha
__inline void L_theta(complex<double>* L, double theta, parameters p, double ci, double z, double z_p) {
    double ni2sin2theta = p.ni_2*sin(theta)*sin(theta);
    complex<double> sroot = sqrt(complex<double>(p.ns_2 - ni2sin2theta));
    complex<double> groot = sqrt(complex<double>(p.ng_2 - ni2sin2theta));
    complex<double> g0root = sqrt(complex<double>(p.ng0_2 - ni2sin2theta));
    complex<double> i0root = sqrt(complex<double>(p.ni0_2 - ni2sin2theta));
    L[0] = p.ni*(ci - z)*cos(theta) + z_p*sroot + p.tg*groot - p.tg0*g0root - p.ti0*i0root;
    L[1] = p.ni*sin(theta) * (z-ci + p.ni*cos(theta)*(p.tg0/g0root + p.ti0/i0root - p.tg/groot - z_p/sroot));
}


// Evaluates the optical path difference, together with its partial derivative d/d_rho in rho
__inline void L_rho(complex<double>* L, double rho, parameters p, double ci, double z, double z_p) {
    double NA2rho2 = p.NA*p.NA*rho*rho;
    complex<double> iroot = sqrt(complex<double>(p.ni_2 - NA2rho2));
    complex<double> sroot = sqrt(complex<double>(p.ns_2 - NA2rho2));
    complex<double> groot = sqrt(complex<double>(p.ng_2 - NA2rho2));
    complex<double> g0root = sqrt(complex<double>(p.ng0_2 - NA2rho2));
    complex<double> i0root = sqrt(complex<double>(p.ni0_2 - NA2rho2));
    L[0] = (ci - z)*iroot + z_p*sroot + p.tg*groot - p.tg0*g0root - p.ti0*i0root;
    L[1] = 2.0*p.NA*p.NA*rho * ((z-ci)/iroot - z_p/sroot - p.tg/groot + p.tg0/g0root + p.ti0/i0root);
}


void parseParameterStruct(parameters &p, const mxArray* pm) {
    
    mxArray* ma;
    
    ma = mxGetField(pm, 0, "ti0");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ti0' (case sensitive).\n");
    p.ti0 = mxGetScalar(ma);
    ma = mxGetField(pm, 0, "ni0");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ti0' (case sensitive).\n");
    p.ni0 = mxGetScalar(ma);
    p.ni0_2 = p.ni0*p.ni0;
    ma = mxGetField(pm, 0, "ni");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ni' (case sensitive).\n");
    p.ni  = mxGetScalar(ma);
    p.ni_2 = p.ni*p.ni;
    ma = mxGetField(pm, 0, "tg0");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'tg0' (case sensitive).\n");
    p.tg0  = mxGetScalar(ma);
    ma = mxGetField(pm, 0, "tg");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'tg' (case sensitive).\n");
    p.tg  = mxGetScalar(ma);
    ma = mxGetField(pm, 0, "ng0");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ng0' (case sensitive).\n");
    p.ng0  = mxGetScalar(ma);
    p.ng0_2 = p.ng0*p.ng0;
    ma = mxGetField(pm, 0, "ng");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ng' (case sensitive).\n");
    p.ng  = mxGetScalar(ma);
    p.ng_2 = p.ng*p.ng;
    ma = mxGetField(pm, 0, "ns");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'ns' (case sensitive).\n");
    p.ns  = mxGetScalar(ma);
    p.ns_2 = p.ns*p.ns;
    ma = mxGetField(pm, 0, "lambda");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'lambda' (case sensitive).\n");
    p.lambda  = mxGetScalar(ma);
    p.k0 = 2*PI/p.lambda;
    ma = mxGetField(pm, 0, "M");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'M' (case sensitive).\n");
    p.M  = mxGetScalar(ma);
    ma = mxGetField(pm, 0, "NA");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'NA' (case sensitive).\n");
    p.NA  = mxGetScalar(ma);
    p.NA_2 = p.NA*p.NA;
    p.alpha = asin(p.NA/p.ni);
    ma = mxGetField(pm, 0, "pixelSize");
    if (ma==NULL) mexErrMsgTxt("Parameter input 'p' must contain field 'pixelSize' (case sensitive).\n");
    p.pixelSize  = mxGetScalar(ma);
    if (p.ni <= p.NA) mexErrMsgTxt("Refractive index must be strictly larger than NA.");
    ma = mxGetField(pm, 0, "sf");
    if (ma!=NULL) {
        p.sf = (int)mxGetScalar(ma);
    } else {
        p.sf = 3;
    }
    ma = mxGetField(pm, 0, "mode");
    if (ma!=NULL) {
        p.mode = (int)mxGetScalar(ma);
    } else {
        p.mode = 1;
    }
}

} // namespace std
#endif // PSFMATH_H
