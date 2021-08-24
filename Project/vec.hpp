#include <iostream>
#include <iomanip>
#include<cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <list>

using namespace std;

mt19937 mt;
uniform_real_distribution<double> xsi(0.0,1.0);
normal_distribution<double> gaus(0.0,1);
//Class defining the vector properties, it will work in 3 dimension 
// All the values will be considered as double 
class MyVec{
    public:

    double *r;

    //empty vector
    MyVec(){
        r=new double[3];
        for(int i=0;i<3;i++){
            r[i]=0.;
        }
    }

    //isotropic vector with strength l
    MyVec(double l){
        r=new double[3];
        for(int i=0;i<3;i++){
            r[i]=1.*l;
        }
    }

    //spcified vector
    MyVec(double x,double y,double z){
        r=new double[3];
        r[0]=x;
        r[1]=y;
        r[2]=z;
    };

    //initialize vector with precise coordinates
    void set_initial(double x,double y,double z){
        r[0]=x;
        r[1]=y;
        r[2]=z;
    };

    void set_unit(){
        for(int i=0;i<3;i++){
            r[i]=1.;
        }
    }

    //initialize vector with random between 0 and 1
    //defined in the programe according to the generator
    void set_random(){
        for(int i=0;i<3;i++){
            r[i]=xsi(mt);
        }
    };
    void set_random(double l){
        for(int i=0;i<3;i++){
            r[i]=l*xsi(mt);
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
            r[i]-=rint(r[i]/l)*l;
        }    
    }
    //copy a vector
    void copy(MyVec&v){
        for(int i=0;i<3;i++){
            r[i]=v.r[i];
        }
    }
    //norm of the vector
    double norm(){
        double R=0;
        for(int i=0;i<3;i++){
            R+=r[i]*r[i];
        }
        return sqrt(R);
    };
    //norm with PB conditions
    double norm(double l){
        double R=0.;
        double _r=0;
        for(int i=0;i<3;i++){
            _r=r[i]-l*rint(r[i]/l);
            R+=_r*_r;
        }
        return sqrt(R);
    };
    //sum between two vector
    MyVec operator+(MyVec&v){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i]=r[i]+v.r[i];
        }
        return vr;
    };

    //difference between vectors
    MyVec operator-(MyVec&v){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i]=r[i]-v.r[i];
        }
        return vr;
    };

    MyVec operator-(){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i]=r[i]*(-1.);
        }
        return vr;
    };
    
    //scalar addition 
    MyVec operator &&(double l){
        MyVec vr(l);
        for(int i =0;i<3;i++){
            vr.r[i]+=r[i];
        }
        return vr;
    };

    //product with a scalar
    MyVec operator *(double l){
        MyVec vr;
        for(int i =0;i<3;i++){
            vr.r[i]=l*r[i];
        }
        return vr;
    };

    //dot product 
    double operator &(MyVec &v){
        double d=0.;
        for(int i=0;i<3;i++){
            d+=r[i]*v.r[i];
        }
        return d;
    };

    //external product
    MyVec operator ^(MyVec &v){
        MyVec vr;
        for(int i=0;i<3;i++){
            vr.r[i]=r[(i+1)%3]*v.r[(i+2)%3]-(r[(i+2)%3]*v.r[(i+1)%3]);
        }
        return vr;
    };
    //distance 
    double dist(MyVec & v){
        MyVec vr;
        double d;
        vr=(*this)-v;
        d=vr.norm();
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
    double dist(MyVec & v,double l){
        MyVec vr;
        double d;
        for(int i=0;i<3;i++){
            d=r[i]-v.r[i];
            vr.r[i]=d-rint(d/l)*l;
        }
        d=vr.norm();
        return d;
    }
    //print vec
    void print_vec(){
        for(int i=0;i<3;i++){
            cout<<r[i]<<' ';
        }
        cout<<endl;
    };
};