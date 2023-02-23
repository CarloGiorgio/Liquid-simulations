#include <iostream>
#include <algorithm>
#include <iomanip>
#include<cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <list>

using namespace std;

//useful functions for making natural power
double MyPow(double a,int n){
    double r=1;
    for(int i=0;i<n;i++){
        r*=a;
    }
    return r;
}

int MyPow(int a,int n){
    int i=1;
    for(int j=0;j<n;j++){
        i*=a;
    }
    return i;
}

//definition of pseudo-random numbers
//generator
mt19937 mt;
//type of distributions
uniform_real_distribution<double> xsi(0.0,1.0);
normal_distribution<double> gaus(0.0,1);

//Class defining vector algebra,it is in 3d 
//It is based on doubles
class MyVec{
    public:

    double *r;

    //empty vector
    MyVec(){
        r = new double[3];
        for(int i=0;i<3;i++){
            r[i]=0.;
        }
    }

    //isotropic vector with strength l
    MyVec(double l){
        r = new double[3];
        for(int i=0;i<3;i++){
            r[i]=1.*l;
        }
    }

    //spcified vector
    MyVec(double x,double y,double z){
        r = new double[3];
        r[0]=x;
        r[1]=y;
        r[2]=z;
    };

    //initialize vector with precise coordinates
    void set_initial(double x,double y,double z){
        r[0] = x;
        r[1] = y;
        r[2] = z;
    };

    void set_unit(){
        for(int i=0;i<3;i++){
            r[i] = 1.;
        }
    }

    //initialize vector with random between 0 and 1
    //defined in the programe according to the generator
    void set_random(){
        for(int i=0;i<3;i++){
            r[i] = xsi(mt);
        }
    };

    void set_random(double l){
        for(int i=0;i<3;i++){
            r[i] = l*xsi(mt);
        }
    };

    //random vector in uniform volume between -l/2 and l/2
    void set_random_2(){
        for(int i=0;i<3;i++){
            r[i] = (xsi(mt)-0.5);
        }
    };

    void set_random_2(double l){
        for(int i=0;i<3;i++){
            r[i]=(xsi(mt)-0.5)*l;
        }
    };

    //uniform vector on a sphere of radius 1
    void set_random_sphere(){
        double  xisq, xi1, xi2, xi;
        xisq = 1.0;
        while (xisq >= 1.0)
            {
            xi1  = 1.0 - 2.0*xsi(mt);
            xi2  = 1.0 - 2.0*xsi(mt);
            xisq = xi1 * xi1 + xi2 * xi2;
            }

        xi = sqrt(fabs(1.0 - xisq));
        r[0] = 2.0 * xi1 * xi;
        r[1] = 2.0 * xi2 * xi;
        r[2] = 1.0 - 2.0 * xisq;
    };
    // Periodicaly Boundary conditions
    void set_PB(double l){
        for (int i=0;i<3;i++){
            r[i] -= rint(r[i]/l)*l;
        }    
    }
    //copy a vector
    void copy(MyVec&v){
        for(int i=0;i<3;i++){
            r[i] = v.r[i];
        }
    }
    //norm of the vector
    double norm(){
        double R = 0.;
        for(int i=0;i<3;i++){
            R += r[i]*r[i];
        }
        return sqrt(R);
    };
    //norm with PB conditions
    double norm(double l){
        double R=0.;
        double _r=0;
        for(int i=0;i<3;i++){
            _r = r[i]-l*rint(r[i]/l);
            R += _r*_r;
        }
        return sqrt(R);
    };
    //sum between two vector
    MyVec operator+(MyVec&v){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i] = r[i] + v.r[i];
        }
        return vr;
    };

    //difference between vectors
    MyVec operator-(MyVec&v){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i] = r[i] - v.r[i];
        }
        return vr;
    };

    MyVec operator-(){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i] = r[i]*(-1.);
        }
        return vr;
    };
    
    //addition by a scalar
    MyVec operator &&(double l){
        MyVec vr(l);
        for(int i =0;i<3;i++){
            vr.r[i] += r[i];
        }
        return vr;
    };

    //product with a scalar
    MyVec operator *(double l){
        MyVec vr;
        for(int i =0;i<3;i++){
            vr.r[i] = l*r[i];
        }
        return vr;
    };

    MyVec operator /(double l){
        MyVec vr;
        for(int i =0;i<3;i++){
            vr.r[i] = r[i]/l;
        }
        return vr;
    };

    //scalar product 
    double operator &(MyVec &v){
        double d=0.;
        for(int i=0;i<3;i++){
            d+=r[i]*v.r[i];
        }
        return d;
    };

    
    MyVec operator ^(MyVec &v){
        //external product
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i]=r[(i+1)%3]*v.r[(i+2)%3]-(r[(i+2)%3]*v.r[(i+1)%3]);
        }
        return vr;
    };
    
    
    double dist(MyVec & v){
        //distance between vector
        MyVec vr;
        double d;
        vr = (*this) - v;
        d = vr.norm();
        return d;
    }

    //minimum image convention
    MyVec MIC(MyVec &v,double l){
        MyVec w;
        double rij;
        for(int i=0;i<3;i++){
            rij=r[i]-v.r[i];
            w.r[i]=rij-l*rint(rij/l);
        }
        return w;
    }
    //distance between vector in PBC
    double dist(MyVec & v,double l){
        MyVec vr;
        double d;
        for(int i=0;i<3;i++){
            d= r[i] - v.r[i];
            vr.r[i] = d - rint(d/l)*l;
        }
        d = vr.norm();
        return d;
    }

    MyVec project(MyVec &v){
        //perform the projection of a vector on r
        MyVec vr;
        double d;
        d = (*this)&v;
        for(int i = 0;i<3; i++){
            vr.r[i] = r[i] * d;
        }
        return vr;
    };

    //print vector on stdout
    void print_vec(){
        for(int i=0;i<3;i++){
            cout<<r[i]<<' ';
        }
        cout<<endl;
    };

    //print on a file
    void print_vec(ostream&f){
          for(int i=0;i<3;i++){
            f<<r[i]<<' ';
        }
        f<<endl;
    };
};