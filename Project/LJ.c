# include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
double L;
double e;
int N;
double p;
double** pos;
double delta;
double rc;
double u_shift;
double V;
double rho;
double T;
double p_tail;
double u_tail;

double MyPow(double a,int n){
    double r=1;
    for(int i=0;i<n;i++){
        r*=a;
    }
    return r;
}

void tail(){
    double r6,r3;
    r3=1./MyPow(rc,3);
    r6=r3*r3;
    p_tail=32./9.*M_PI*rho*rho*r3*(r6-1.5);
    u_tail=8./9.*M_PI*rho*r3*(r6-3.);
}
void set_random(){
    for(int i=0;i<N;i++){
        for(int j=0;j<3;j++){
            pos[i][j]=drand48()*L;
        }
    }
}

double dist_PB (double * v,double *w){
    double d=0.;
    double r;
    for(int i=0;i<3;i++){
        r=v[i]-w[i];
        r=r-rint(r/L)*L;
        d+=MyPow(r,2);
    }
    return sqrt(d);
}

void displace(int n,double *p){
    for(int i=0;i<3;i++){
        p[i]=pos[n][i]+(drand48()-0.5)*delta;
    }
}

double potential(double *v,double *w){
    double d,d12,d6;
    d=dist_PB(v,w);
    d6=1./MyPow(d,6);
    d12=d6*d6;
    if(d<rc){
        return 4*(d12-d6)- u_shift;
    }else {
        return 0;
    }
}

double energy_displace(double *v,int n){
    double _e=0.;
    for(int i=0;i<N;i++){
        if(i!=n){
            _e+=potential(v,pos[i]);
        }
    }
    return _e;
}

void energy_pressure(){
    e=u_tail;
    p=p_tail;
    double d,d6,d12;
    for(int i=0;i<N-1;i++){
        for(int j=i+1;j<N;j++){
            d=dist_PB(pos[i],pos[j]);

            if(d<rc){
                d6=1./MyPow(d,6);
                d12=d6*d6;
                p+=48*(d12-0.5*d6)/V/3.;
                e+=(4*(d12-d6)-u_shift)/N;
            }

        }
    }
}

void copy(double *v,double*w){
    for(int i=0;i<3;i++){
        v[i]=w[i];
    }
}
void print(double *v){
    for(int i=0;i<3;i++){
        printf("%e ",v[i]);
    }
    printf("\n");
}
void mc_step(){
    double *p;
    p=(double*)malloc(3*sizeof(double));
    double d;
    int n=0;
    for(int i=0;i<N;i++){
        n=floor(drand48()*(N-1));
        displace(n,p);
        d=energy_displace(p,n)-energy_displace(pos[n],n);
        if(d<=0||drand48()<exp(-d/T)){
            copy(pos[n],p);
        }
    }
}

int main(){
    N=200;
    rho=0.5;
    V=N/rho;
    L=pow(V,1./3.);
    rc=2.5;
    T=0.9;
    delta=0.2;
    u_shift=4.0*(1.0/MyPow(rc,12)-1.0/MyPow(rc,6));
    tail();
    unsigned long int seed;
    seed=time(NULL);
    srand48(seed);

    pos=(double**)malloc(N*sizeof(double *));
    for(int i=0;i<N;i++){
        pos[i]=(double *)malloc(3*sizeof(double));
    }
    set_random();
    energy_pressure();
    //printf("%.2e %.2e\n",e,p);
    for(int i=1;i<3000;i++){
        mc_step();
        if(i >1000){
            energy_pressure();
            printf("%d %e %e\n",i,e,p+rho*T);
        }
    }
    return 0;
}