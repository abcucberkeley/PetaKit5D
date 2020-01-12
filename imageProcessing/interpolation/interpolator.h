/* B-spline based interpolation. Implementation for 1D and 2D signals
 *
 * References
 * [1] Unser, IEEE Signal Proc. Mag. 16(6), pp. 22-38, 1999
 * [2] Unser et al., IEEE Trans. Signal Proc. 41(2), pp. 834-848, 1993
 *
 * Francois Aguet ? May 2012 (last modified May 11, 2012)
 */

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

namespace std {

class Interpolator {
    
 public:
    
    typedef void (*SplineFctPtr)(const double t, double v[]);
    
    // Border conditions
    static const int PERIODIC = 2;
    static const int MIRROR = 3;
    
    Interpolator(const double pixels[], const int nx, const int ny, const int borderCondition, const int degree);
    Interpolator(const double pixels[], const int nx, const int ny, const int borderCondition);
    Interpolator(const double pixels[], const int nx, const int ny);
    void init(const double pixels[], const int nx, const int ny, const int borderCondition, const int degree);
    
    ~Interpolator();
    
    double* getCoefficients() { return coefficients_; };
    // 1D interpolation
    void interp(const double x[], const int nxi, double values[], int n);
    void interp(const double x[], const int nxi, double values[]);
    // 2D interpolation
    void interp(const double x[], const double y[], const int N, double values[], int n);
    void interp(const double x[], const double y[], const int N, double values[]);
    
 private:
    int nx_, ny_, N_, dims_;
    int degree_; // degree of the spline. 1: linear, 2: quadratic, 3: cubic
    double a_, c0_;
    
    const double* pixels_;
    double* coefficients_;
    double* cp_, *cn_; // used for coefficient calculation
    double* xi_, *yi_; // interpolation coordinates    
    int nxi_, nyi_;
    
    int borderCondition_;
    
    SplineFctPtr splineFctPtr_;
    
    void getCausalInitValueX();
    void getCausalInitValueY();
    void getAntiCausalInitValueX();
    void getAntiCausalInitValueY();
    
    void computeCoefficients();
    int wrap(int i, int N);
    
    static void getCubicSpline(const double t, double v[]);
    static void getQuadraticSpline(const double t, double v[]);

};
    
    //const double Interpolator::a_ = -2.0+sqrt(3.0);
    //const double Interpolator::c0_ = 6.0;
    
    
    void Interpolator::init(const double pixels[], const int nx, const int ny, const int borderCondition, const int degree) {

        pixels_ = pixels;
        nx_ = nx;
        ny_ = ny;
        N_ = nx*ny;
        if (nx>1 && ny>1) {
            dims_ = 2;
        } else {
            dims_ = 1;
        }
        borderCondition_ = borderCondition;
        degree_ = degree;
        
        coefficients_ = new double[N_];
        cp_ = new double[N_];
        cn_ = new double[N_];
        
        splineFctPtr_ = &getCubicSpline; // default: cubic spline interpolation
        
        // set filter coefficients
        switch (degree_) {
            case 3:
                a_ = -2.0+sqrt(3.0);
                c0_ = 6.0;
                break;
            case 2:
                a_ = -3.0+2.0*sqrt(2.0);
                c0_ = 8;
                splineFctPtr_ = &getQuadraticSpline;
                break;
            default:
                break;
        }
        
        if (degree>1) {
            computeCoefficients();
        } else {
            memcpy(coefficients_, pixels_, N_*sizeof(double));
        }
    }
    
    
    Interpolator::Interpolator(const double pixels[], const int nx, const int ny, const int borderCondition, const int degree) {
        init(pixels, nx, ny, borderCondition, degree);
    }
    
    
    Interpolator::Interpolator(const double pixels[], const int nx, const int ny, const int borderCondition) {
       init(pixels, nx, ny, borderCondition, 3); // Default: cubic spline
    }
    
    
    Interpolator::Interpolator(const double pixels[], const int nx, const int ny) {
        init(pixels, nx, ny, MIRROR, 3);
    }
            
    
    Interpolator::~Interpolator() {
        delete coefficients_;
        delete cp_;
        delete cn_;
    }
    
    
    void Interpolator::computeCoefficients() {
        int i;
        
        // First filter on x-dimension
        getCausalInitValueX();
        for (int y=0;y<ny_;++y) {
            for (int x=1;x<nx_;++x) {
                i = x+y*nx_;
                cp_[i] = pixels_[i] + a_*cp_[i-1];
            }
        }
        getAntiCausalInitValueX();
        for (int y=0;y<ny_;++y) {
            for (int x=nx_-2;x>=0;--x) {
                i = x+y*nx_;
                cn_[i] = a_*(cn_[i+1]-cp_[i]);
            }
        }
        
        // If signal is 2D, compute coefficients on y-dimension
        if (dims_==2) {
            getCausalInitValueY();
            for (int x=0;x<nx_;++x) {
                for (int y=1;y<ny_;++y) {
                    i = x+y*nx_;
                    cp_[i] = cn_[i] + a_*cp_[i-nx_]; //i-1 on y becomes i-nx
                }
            }
            getAntiCausalInitValueY();
            for (int x=0;x<nx_;++x) {
                for (int y=ny_-2;y>=0;--y) {
                    i = x+y*nx_;
                    cn_[i] = a_*(cn_[i+nx_]-cp_[i]);
                }
            }
        }
            
        // Multiply by constant component of prefilter
        double c0 = c0_;
        if (dims_==2) {
            c0 *= c0;
        }
        for (int i=0;i<N_;++i) {
            coefficients_[i] = c0*cn_[i];
        }
    }
    
    
    // 1D interpolation
    // Note: degree of interpolating spline/coefficients can be different (i.e., when computing derivatives)
    void Interpolator::interp(const double x[], const int N, double v[], int n) {
        
        SplineFctPtr getSpline = &getCubicSpline; // default
        if (n==2) {
            getSpline = &getQuadraticSpline;
        }
        
        int xi;
        double dx;
        double b[4];
        if (n==3 || n==2) {
            for (int i=0;i<N;++i) {
                xi = (int)floor(x[i]);
                dx = x[i] - xi;
                getSpline(dx, b);
                v[i] = b[0]*coefficients_[wrap(xi-1,nx_)] + b[1]*coefficients_[wrap(xi,nx_)] + b[2]*coefficients_[wrap(xi+1,nx_)] + b[3]*coefficients_[wrap(xi+2,nx_)];
            }
        } else { // linear interpolation
            int x0, x1;
            for (int i=0;i<N;++i) {
                x0 = (int)floor(x[i]);
                dx = x[i] - x0;
                x1 = wrap(x0+1, nx_);
                x0 = wrap(x0, nx_);
                v[i] = dx*coefficients_[x1] + (1.0-dx)*coefficients_[x0];
            }
        }
    }
    
    
    void Interpolator::interp(const double x[], const int N, double v[]) {
        interp(x, N, v, degree_);
    }
    
    
    // 2D interpolation
    void Interpolator::interp(const double x[], const double y[], const int N, double v[], int n) {
        
        SplineFctPtr getSpline = &getCubicSpline; // default
        if (n==2) {
            getSpline = &getQuadraticSpline;
        }
        
        int xi, yi;
        double dx, dy;

        if (n==3 || n==2) {
            double bx[4], by[4];
            int row[4], col[4];
            for (int i=0;i<N;++i) {
                xi = (int)floor(x[i]);
                yi = (int)floor(y[i]);
                dx = x[i]-xi;
                dy = y[i]-yi;
                
                // row index
                row[0] = wrap(yi-1,ny_);
                row[1] = wrap(yi,ny_);
                row[2] = wrap(yi+1,ny_);
                row[3] = wrap(yi+2,ny_);
                
                col[0] = wrap(xi-1,nx_);
                col[1] = wrap(xi,nx_);
                col[2] = wrap(xi+1,nx_);
                col[3] = wrap(xi+2,nx_);
                
                // get x-interpolated value for each row
                getSpline(dx, bx); // CHANGE TO FCT POINTER
                getSpline(dy, by);
                v[i] = by[0]*(bx[0]*coefficients_[col[0]+row[0]*nx_] + bx[1]*coefficients_[col[1]+row[0]*nx_] +
                              bx[2]*coefficients_[col[2]+row[0]*nx_] + bx[3]*coefficients_[col[3]+row[0]*nx_]) + 
                       by[1]*(bx[0]*coefficients_[col[0]+row[1]*nx_] + bx[1]*coefficients_[col[1]+row[1]*nx_] + 
                              bx[2]*coefficients_[col[2]+row[1]*nx_] + bx[3]*coefficients_[col[3]+row[1]*nx_]) + 
                       by[2]*(bx[0]*coefficients_[col[0]+row[2]*nx_] + bx[1]*coefficients_[col[1]+row[2]*nx_] + 
                              bx[2]*coefficients_[col[2]+row[2]*nx_] + bx[3]*coefficients_[col[3]+row[2]*nx_]) + 
                       by[3]*(bx[0]*coefficients_[col[0]+row[3]*nx_] + bx[1]*coefficients_[col[1]+row[3]*nx_] + 
                              bx[2]*coefficients_[col[2]+row[3]*nx_] + bx[3]*coefficients_[col[3]+row[3]*nx_]);
            }
        } else { // linear interpolation
            int x0, x1, y0, y1;
            for (int i=0;i<N;++i) {
                x0 = (int)floor(x[i]);
                y0 = (int)floor(y[i]);
                dx = x[i]-x0;
                dy = y[i]-y0;
                
                x0 = wrap(x0, nx_);
                x1 = wrap(x0+1, nx_);
                y0 = wrap(y0, ny_);
                y1 = wrap(y0+1, ny_);
                
                // x00 = x0+y0*nx
                // x10 = x1+y0*nx
                // x01 = x0+y1*nx
                // x11 = x1+y1*nx
                // top row:    dx*pixels_[x10] + (1.0-dx)*pixels_[x00];
                // bottom row: dx*pixels_[x11] + (1.0-dx)*pixels_[x01];
                v[i] = (1.0-dy)*(dx*pixels_[x1+y0*nx_] + (1.0-dx)*pixels_[x0+y0*nx_])
                           + dy*(dx*pixels_[x1+y1*nx_] + (1.0-dx)*pixels_[x0+y1*nx_]);                
            }
        }
    }
    
    
    void Interpolator::interp(const double x[], const double y[], const int N, double v[]) {
        interp(x, y, N, v, degree_);
    }
        
        
    void Interpolator::getCausalInitValueX() {
        int i, offset;
        if (borderCondition_ == MIRROR) {
            for (int y=0;y<ny_;++y) { // loop through rows
                i = y*nx_;
                double ap = 1.0; // store powers of 'a'
                cp_[i] = 0.0;
                for (int k=0;k<nx_;++k) {
                    cp_[i] += pixels_[i+k]*ap;
                    ap *= a_;
                }
                for (int k=nx_-2;k>0;--k) { // mirror: loop backwards
                    cp_[i] += pixels_[i+k]*ap;
                    ap *= a_;
                }
                cp_[i] /= (1.0 - ap);
            }
        } else {
            for (int y=0;y<ny_;++y) { // loop through rows
                i = y*nx_;
                offset = i+nx_;
                double ap = a_; // store powers of 'a'
                cp_[i] = pixels_[i];
                for (int k=1;k<nx_;++k) {
                    cp_[i] += pixels_[offset-k]*ap;
                    ap *= a_;
                }
                cp_[i] /= (1.0 - ap);
            }
        }
    }
    
    
    void Interpolator::getCausalInitValueY() {
        if (borderCondition_ == MIRROR) {
            for (int x=0;x<nx_;++x) { // loop through columns
                double ap = 1.0; // store powers of 'a'
                cp_[x] = 0.0;
                for (int k=0;k<ny_;++k) {
                    cp_[x] += cn_[x+k*nx_]*ap;
                    ap *= a_;
                }
                for (int k=ny_-2;k>0;--k) { // mirror: loop backwards
                    cp_[x] += cn_[x+k*nx_]*ap;
                    ap *= a_;                    
                }
                cp_[x] /= (1.0 - ap);
            }
        } else {
            int offset;
            for (int x=0;x<nx_;++x) {
                offset = (ny_)*nx_+x;
                double ap = a_;
                cp_[x] = cn_[x];
                for (int k=1;k<ny_;++k) {
                    cp_[x] += cn_[offset-k*nx_]*ap;
                    ap *= a_;
                }
                cp_[x] /= (1.0 - ap);
            }
        }
    }
    
    
    void Interpolator::getAntiCausalInitValueX() {
        int i;
        if (borderCondition_ == MIRROR) {
            for (int y=0;y<ny_;++y) { // loop through rows
                i = (y+1)*nx_-1; // index of last point on row
                cn_[i] = (a_/(a_*a_-1.0)) * (cp_[i] + a_*cp_[i-1]); // [1], Box 2 (has errors)
            }
        } else { // periodic
            for (int y=0;y<ny_;++y) {
                i = (y+1)*nx_-1;
                double ap = a_;
                cn_[i] = cp_[i];
                for (int k=0;k<nx_-1;++k) {
                    cn_[i] += ap*cp_[y*nx_+k];
                    ap *= a_;
                }
                cn_[i] *= -a_/(1.0-ap);
            }
        }
    }
    
    
    void Interpolator::getAntiCausalInitValueY() {
        int i;
        int offset = (ny_-1)*nx_;
        if (borderCondition_ == MIRROR) {
            for (i=offset;i<N_;++i) { // loop through last row
                cn_[i] = (a_/(a_*a_-1.0)) * (cp_[i] + a_*cp_[i-nx_]);
            }
        } else { // periodic
            for (int x=0;x<nx_;++x) {
                i = (ny_-1)*nx_+x; // index of last point on columns (= last row)
                double ap = a_;
                cn_[i] = cp_[i];
                for (int k=0;k<ny_-1;++k) {
                    cn_[i] += ap*cp_[k*nx_+x]; // s[N-1] + a*s[0] + a^2*s[2] + ...
                    ap *= a_;
                }
                cn_[i] *= -a_/(1.0-ap);
            }            
        }
    }
    
    
    void Interpolator::getCubicSpline(const double t, double v[]) {
        double t1 = 1.0 - t;
        double t2 = t*t;
        v[0] = (t1 * t1 * t1) / 6.0;
        v[1] = (2.0 / 3.0) + 0.5 * t2 * (t-2.0);
        v[3] = (t2 * t) / 6.0;
        v[2] = 1.0 - v[3] - v[1] - v[0];
    }   
    
    
    void Interpolator::getQuadraticSpline(const double t, double v[]) {
        if (t<=0.5) {
            v[0] = (t-0.5)*(t-0.5)/2.0;
            v[1] = 0.75 - t*t;
            v[2] = 1.0-v[1]-v[0];
            v[3] = 0.0;
        } else {
            v[0] = 0.0;
            v[1] = (t-1.5)*(t-1.5)/2.0;
            v[3] = (t-0.5)*(t-0.5)/2.0;
            v[2] = 1.0-v[3]-v[1];
        }
    }
    
    
    // If rounded interpolation coordinates fall outside signal range, apply border conditions
    int Interpolator::wrap(int i, int N) {
        if (borderCondition_ == MIRROR) {
            if (i<0) {
                return -i;
            } else if (i>=N) {
                return 2*N-2 - i;
            } else {
                return i;
            }
        } else { // periodic
            if (i<0) {
                return i+N;
            } else if (i>=N) {
                return i-N;
            } else {
                return i;
            }
        }
    }
    
}
#endif // INTERPOLATOR_H
