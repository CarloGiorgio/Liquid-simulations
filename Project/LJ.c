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
double beta;
double Pv;
double dv;

double MyPow(double a,int n){
    double r=1;
    for(int i=0;i<n;i++){
        r*=a;
    }
    return r;
}

//tail correction
void tail(int i){
    if (i==0){
        p_tail=0.;
        u_tail=0.;
    }else{
        double r6,r3;
        r3=1./MyPow(rc,3);
        r6=r3*r3;
        p_tail=32./9.*M_PI*rho*rho*r3*(r6-1.5);
        u_tail=8./9.*M_PI*rho*r3*(r6-3.);
    }
}

//random configuration
void set_random(){
    for(int i=0;i<N;i++){
        for(int j=0;j<3;j++){
            pos[i][j]=drand48()*L;
        }
    }
}

//minimum image convention for the distance
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

//random vector between -delta/2 and delta/2
double* rand_displace(){
    double *p;
    p=(double*)malloc(3*sizeof(double));
    for(int i=0;i<3;i++){
        p[i]=(drand48()-0.5)*delta;
    }
    return p;
}
double LJ(double r){
    return 0*4*(1.0/MyPow(r,12)-1.0/(MyPow(r,6)));
}
double potential(double *v,double *w){
    double d;
    d=dist_PB(v,w);
    if(d<rc){
        return LJ(d)- u_shift;
    }else {
        return 0.0;
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

double energy_tot(){
    double _e;
    _e=0.;
    for(int i=0;i<N-1;i++){
        for(int j=i+1;j<N;j++){
            _e+=potential(pos[i],pos[j]);
        }
    }
    return _e;
}

void energy_pressure(){
    e=u_tail;
    p=p_tail+rho*T;
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

void add (double *v,double*w){
    for(int i=0;i<3;i++){
        v[i]+=w[i];
    }
}
void sub(double *v,double*w){
    for(int i=0;i<3;i++){
        v[i]-=w[i];
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
void set_PB(double* v){
    for(int j=0;j<3;j++){
        v[j]-=rint(v[j]/L)*L;
    }
}

// Function for NPT ensemble
double gibbs_energy(){
    return energy_tot()+Pv*V-N*T*log(V);
}

void scale(double s){
    for(int i=0;i<N;i++){
        for(int j=0;j<3;j++){
            pos[i][j]*=s;
        }
    }
}

void mc_step(){
    double _r=drand48()*(N+1);
    if(_r<N){
        double *p;
        p=(double*)malloc(3*sizeof(double));
        double d;
        int n;
        for(int i=0;i<N;i++){
            n=rint(drand48()*(N-1));
            d=energy_displace(pos[n],n);
            p=rand_displace();
            add(p,pos[n]);
            set_PB(p);
            d=energy_displace(p,n)-d;
            if(d<=0||drand48()<exp(-d*beta)){
                copy(pos[n],p);
                e+=d;
            }
        }
    }else{
        double ran_v;
        double V0,L0,G;
        ran_v=(drand48()-0.5)*dv;
        V0=V;
        G=energy_tot()+Pv*V0-N*log(V0)*T;
        L0=L;
        V=V0*exp(ran_v);
        L=L0*exp(ran_v/3.0);
        scale(L/L0);
        G=energy_tot()+Pv*V-N*log(V)*T-G-log(V/V0)*T;
        if(G>0&&drand48()>exp(-beta*G)){
            scale(L0/L);
            V=V0;
            L=L0;
        }
        rho=N/V;
    }
}

int main(){
    N=400;
    rho=0.2;
    V=N/rho;
    L=pow(V,1./3.);
    rc=-2.5;
    T=1.4;
    beta=1.0/T;
    delta=0.2;
    dv=0.5;
    Pv=2.0;
    u_shift=LJ(rc);
    tail(0);
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
    for(int i=1;i<10001;i++){
        mc_step();
        //printf("%.2lf\n",L);

        if(i >7000){
            energy_pressure();
            printf("%d %g %g %g\n",i,e,p+rho*T,rho);
        }
    }
    return 0;
}