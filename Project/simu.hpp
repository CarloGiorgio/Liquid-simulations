#include "vec.hpp"

uniform_int_distribution<> n_mc;

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
//Class that contains the initialization of the vector and 
class MySet{
    public:
    //configuration variables
    int N;
    double rho;
    double V;
    double L;
    vector<MyVec> pos;

    //energy variables
    double e;
    double rc;
    double u_shift;
    double u_tail;
    double p_tail;

    //radial distribution variables
    vector<double> g_r;
    vector <double> _r;
    int n_bin;
    double dg;
    double n_count;
    
    //LL variables
    int M;
    double rn;
    vector <list <int> > LL;
    vector<vector<int> > NN_list;
    vector <vector <int> > pos_cell;
    vector <int> pos_cell_1;

    //Verlet List
    double rv;
    vector <list <int> > VL;
    vector <MyVec> pos_v;

    void set_value(double _rho,int n){
        double rc6,rc3;
        rho=_rho;
        N=n;
        V=N/rho;
        L=pow(V,1./3.);
        pos.resize(N);
        e=0.;
        rc=2.5;
        rc3=1.0/MyPow(rc,3);
        rc6=rc3*rc3;
        u_shift=4.*rc6*(rc6-1.);
        u_tail=8./9.*M_PI*rho*rc3*(rc6-3.0);
        p_tail=32.0/9.*M_PI*rho*rho*rc3*(rc6-3./2.);
        n_mc=uniform_int_distribution<>(0,N-1);
    };

    void set_value(double _rho,int n,double r_c){
        double rc6,rc3;
        rho=_rho;
        N=n;
        V=N/rho;
        L=pow(V,1./3.);
        pos.resize(N);
        e=0.;
        rc=r_c;
        rc3=1.0/MyPow(rc,3);
        rc6=1.0/MyPow(rc,6);
        u_shift=4.*rc6*(rc6-1.0);
        u_tail=8./9.*M_PI*rho*rc3*(rc6-3.0);
        p_tail=32.0/9.*M_PI*rho*rho*rc3*(rc6-3./2.);
        n_mc=uniform_int_distribution<>(0,N-1);
    };

    //g(r) setting; see Frenkel
    //creation of variables for g(r)
    void set_g_r(int _n){
        n_bin=_n;
        g_r.resize(n_bin);
        _r.resize(n_bin);
        dg=L/2./n_bin;
        n_count=0.;
    };

    //sampling for g(r)
    void add_g_r(){
        double r_g;
        int i_g;
        for(int i=0;i<N-1;i++){
            for(int j=i+1;j<N;j++){
                r_g=pos[i].dist(pos[j],L);
                if(r_g<L/2.){
                    i_g=rint(r_g/dg);
                    g_r[i_g]+=2;
                }
            }
        }
        n_count++;
    };

    //evaluation of g(r)
    void eval_g_r(){
        double nd=4./3.*M_PI*rho;
        double vb;
        for(int i=0;i<g_r.size();i++){
            _r[i]=dg*(i+0.5);
            vb=(MyPow(i+1,3)-MyPow(i,3))*MyPow(dg,3);
            g_r[i]/=(n_count*N*vb*nd);
        }
    };


    //initial position
    //random
    void pos_random(){
        for(int i=0;i<pos.size();i++){
            pos[i].set_random(L);
        }
    };

    //simple cubic
    void pos_SC(){
        double a=pow(V/N,1./3.);
        int j=0,k=0,l=0;
        for(int i=0;i<pos.size();i++){
            if(j*a>=L){
                j=0;
                k++;
            }

            if(k*a>=L){
                k=0;
                l++;
            }

            pos[i].set_initial(j*a,k*a,l*a);
            j++;
        }
    };

    //Body centered cubic
    void pos_BCC(){
        double a=pow(2.*V/N,1./3.);
        int j=0,k=0,l=0;
        for(int i=0;i<pos.size();i++){
            if(j*a+a/2.*(l%2)>=L){
                j=0;
                k++;
            }

            if(k*a+a/2.*(l%2)>=L){
                k=0;
                l++;
            }

            pos[i].set_initial(j*a+a/2.*(l%2),k*a+a/2.*(l%2),l*a/2.);
            j++;
        }
    };

    //Face centered cubic
    void pos_FCC(){
        double a=pow(4.*V/N,1./3.);
        int j=0,k=0,l=0;
        for(int i=0;i<pos.size();i++){
            if(j*a+a/2.*(k%2+l%2-2*(l%2)*(k%2))>=L){
                j=0;
                k++;
            }
            if(k*a/2.>=L){
                k=0;
                l++;
            }
            pos[i].set_initial(a*j+a/2.*(k%2+l%2-2*(l%2)*(k%2)),a/2.*k,a/2.*l);
            j++;
        }
    };

    //Interaction Potential
    double potential(MyVec &v,MyVec &w){
        double d,d6,d12;
        MyVec _r;
        _r=v.MIC(w,L);
        d=_r.norm();
        if (d<rc){
            d6=1.0/MyPow(d,6);
            d12=d6*d6;
            return 4.*(d12-d6)-u_shift;
        }else{
            return 0.;
        }
    };

    //force
    MyVec force(MyVec &vij){
        MyVec f(0.);
        double d,d2,d6,d12,r_f;
        d=vij.norm();
        if (d<rc){
            d2=1.0/(d*d);
            d6=d2*d2*d2;
            d12=d6*d6;
            r_f=48.*d2*(d12-0.5*d6);
            for(int i=0;i<3;i++){
                f.r[i]=r_f*vij.r[i];
            }
        }
        return f;
    };

    void energy_brute_tot(){
        e=u_tail;
        for(int i=0;i<N-1;i++){
            for(int j=i+1;j<N;j++){
                e+=potential(pos[i],pos[j])/N;
            }
        }
    };

    int find_cell(MyVec v){
        int n=0;
        int n_p;
        for(int i=0;i<3;i++){
            n_p=int(rint(v.r[i]/rn)+M)%M;
            n+=n_p*MyPow(M,i);
        }
        return n;
    };

    //Linked List codes
    void create_LL(){
        rn=rc;
        //rn=rc;
        M=rint(L/rn);
        NN_list.resize(M*M*M);
        LL.resize(M*M*M);
        pos_cell_1.resize(N);

        for(int i=0;i<N;i++){
            pos_cell_1[i]=find_cell(pos[i]);
            LL[pos_cell_1[i]].push_front(i);
        }

        int NN=0;
        for(int n=0;n<NN_list.size();n++){
            for(int i=-1;i<2;i++){
                for(int j=-1;j<2;j++){
                    for(int k=-1;k<2;k++){
                        NN=(n+k+M)%M\
                        +((n/M+j+M)%M)*M\
                        +(((n/M)/M+i+M)%M)*M*M;
                        NN_list[n].push_back(NN);
                    }
                }
            }
        }
    };

    double energy_LL_one(MyVec v,int n){
        double _e=0;
        for(int i : NN_list[pos_cell_1[n]]){
            for(int j:LL[i]){
                if(j!=n){
                    _e+=potential(v,pos[j]);
                }
            }
        }
        return _e;
    };

    void energy_LL(){
        //find a way to evaluate the neighbours and the 
        //loop over the list 
        e=u_tail;
        for(int i=0;i<N;i++){
            e+=energy_LL_one(pos[i],i)/N/2.;
        }
    };

    void create_VL(){
        rv=rc+0.2;
        VL.resize(N);
        double rv_t;
        pos_v.resize(N);

        for(int i=0;i<N;i++){
            pos_v[i].copy(pos[i]);
        }
        for(int i=0;i<N-1;i++){
            
            for(int j=i+1;j<N;j++){
                rv_t=pos_v[i].dist(pos_v[j],L);
                if(rv_t<rv){
                    VL[i].push_front(j);
                    VL[j].push_front(i);
                }
            }
        }  
    };

    void update_VL(){
        for(int i=0;i<VL.size();i++){
            VL[i].clear();
            pos_v[i].copy(pos[i]);
        }

        double rv_t;
        for(int i=0;i<N-1;i++){
            for(int j=i+1;j<N;j++){
                rv_t=pos_v[i].dist(pos_v[j],L);
                if(rv_t<rv){
                    VL[i].push_front(j);
                    VL[j].push_front(i);
                }
            }
        }  
    };
    double energy_VL_one(MyVec v,int n){
        double _e=0.;
        for(int j:VL[n]){
            _e+=potential(v,pos[j]);
        }
        return _e;
    }

    void energy_VL(){
        e=u_tail;
        for(int i=0;i<N;i++){
            e+=energy_VL_one(pos[i],i)/N/2.;
        }
    }

    
};

/*
************************************************************************
*/
// Definition of a class for MC simulation in NTV ensemble, Derived by MySet
class NTV:public MySet{
    public:
    double p;
    double T;
    int time;
    double delta;

    void set_value_NTV(double t,int ti,double d){
        T=t;
        time=ti;
        delta=d;
        p=0.;
    };

    void print_para(){
        cout<<"#N:"<<N<<" rho:"<<rho<<" L:"<<L<<" r_c:"<<rc<<" T:"<<T<<" Delta:" <<delta;
    }

    void print_para(fstream &f){
        f<<"#N:"<<N<<" rho:"<<N/V<<" r_c:"<<rc<<" T:"<<T<<" Delta:" <<delta<<endl;
    }

    void energy_pressure_brute(){
        double d,d6,d12;
        MyVec rij;
        MyVec fij;
        p=p_tail;
        e=u_tail;
        for(int i=0;i<N-1;i++){
            for(int j=i+1;j<N;j++){
                rij=pos[i].MIC(pos[j],L);
                d=pos[i].dist(pos[j],L);

                if (d<rc){
                    d6=MyPow(d,6);
                    d6=1./d6;
                    d12=d6*d6;
                    fij=force(rij);
                    p+=(fij&rij)/V/3.;
                    e+=potential(pos[i],pos[j])/N;
                }
            }
        }
    };

    double energy_brute_one(MyVec v, int n){
        double _e=0.;
        for(int j=0;j<N;j++){
            if(j!=n){
                _e+=potential(v,pos[j]);
            }
        }
        return _e;
    }

    void mc_step_brute(){
        MyVec P;
        int n=0;
        double d=0.;
        double dis;
        //cout<<L<<endl<<endl;
        for(int i=0;i<N;i++){
            n=rint(xsi(mt)*(N-1));
            //cout<<n<<" I'm in"<<endl;
            for(int j=0;j<3;j++){
                dis=(xsi(mt)-0.5)*delta;
                P.r[j]=(pos[n].r[j]+dis);
            }

            P.set_PB(L);
            d=energy_brute_one(P,n)-energy_brute_one(pos[n],n);
            //cout<<d<<" ";
            if(d<=0||xsi(mt)<exp(-d/T)){
                pos[n].copy(P);
            }
        }
    };

    void simu_brute(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,pm=0.,pm2=0.;
        
        f.open(s,ios::out | ios::trunc);
        /*
        if(!f.is_open()){
            f.clear();
            f.open(s, ios::out | ios::trunc); //Create file.
            f.close();
            f.open(s);
        }*/
        print_para(f);
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2 5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        //f<<0<<" "<<e<<" "<<em<<" "<<em2-em*em<<" "<<p<<" "<<pm<<" "<<pm2-MyPow(pm,2)<<endl;
        for(int i=1;i<time+1;i++){
            mc_step_brute();
            if(i>time_start){
                energy_pressure_brute();
                //cout<<e<<" "<<p<<endl;
                em+=e;
                em2+=e*e;
                pm+=p;
                pm2+=p*p;
                if(i/time_take>k){
                    f<<i<<" "<<e<<" "<<em/(i-time_start)<<" "<<sqrt(em2/(i-time_start)-MyPow(em/(i-time_start),2))<<" "<<p<<" "<<pm/(i-time_start)<<" "<<sqrt(pm2/(i-time_start)-MyPow(pm/(i-time_start),2))<<endl;
                    k=i/time_take;
                }
            }
        }
        f.close();
    };

    void energy_pressure_VL(){
        double d,d6,d12;
        p=p_tail;
        e=u_tail;
        for(int i=0;i<N;i++){
            for(int j:VL[i]){
                d=pos[i].dist(pos[j],L);
                d6=MyPow(d,6);
                d6=1./d6;
                d12=d6*d6;
                p+=48.*(d12-0.5*d6)/V/6.;
                e+=(4.*(d12-d6)-u_shift)/2./N;
            }
        }
    };

    void mc_step_VL(){
        MyVec P;
        MyVec _P;
        int n=0;
        double d=0.;
        double dis;
        for(int i=0;i<N;i++){
            n=rint(xsi(mt)*(N-1));
            //cout<<n<<" I'm in"<<endl;
            for(int j=0;j<3;j++){
                dis=(xsi(mt)-0.5)*delta;
                P.r[j]=(pos[n].r[j]+dis);
            }
            _P.copy(pos[n]);
            if(pos[n].dist(pos_v[n],L)<((rv-rc)*0.5)){
                update_VL();
            }
            d=(-1)*energy_VL_one(pos[n],n);
            pos[n].copy(P);

            if(pos[n].dist(pos_v[n],L)<((rv-rc)*0.5)){
                update_VL();
            }

            d+=energy_VL_one(pos[n],n);
            if(d>0&&xsi(mt)>exp(-d/T)){
                pos[n].copy(_P);
            }
        }
    };

    void simu_VL(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,pm=0.,pm2=0.;
        
        f.open(s,ios::out | ios::trunc);
        print_para(f);
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2 5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        
        create_VL();
        for(int i=1;i<time+1;i++){
            mc_step_VL();
            if(i>time_start){
                energy_pressure_VL();
                //cout<<e<<" "<<p<<endl;
                em+=e;
                em2+=e*e;
                pm+=p;
                pm2+=p*p;
                if(i/time_take>k){
                    f<<i<<" "<<e<<" "<<em/(i-time_start)<<" "<<sqrt(em2/(i-time_start)-MyPow(em/(i-time_start),2))<<" "<<p<<" "<<pm/(i-time_start)<<" "<<sqrt(pm2/(i-time_start)-MyPow(pm/(i-time_start),2))<<endl;
                    k=i/time_take;
                }
            }
        }
        f.close();
    };
    void energy_pressure_LL(){
        double d,d12,d6;
        p=p_tail+rho*T;
        e=u_tail;
        for(int i=0;i<N;i++){
            for(int l : NN_list[pos_cell_1[i]]){
                for(int j:LL[l]){
                    if(i!=j){
                        d=pos[i].dist(pos[j],L);
                        if(d<rc){
                            d6=MyPow(d,6);
                            d6=1./d6;
                            d12=d6*d6;
                            p+=48.*(d12-0.5*d6)/V/6.;
                            e+=(4.*(d12-d6)-u_shift)/2./N;
                        }
                    }
                }
            }
        }
    };
    void mc_step_LL(){
        MyVec P;
        int n=0;
        double d=0.;
        double dis;
        for(int i=0;i<N;i++){
            n=n_mc(mt);
            //cout<<n<<" I'm in"<<endl;
            for(int j=0;j<3;j++){
                dis=(xsi(mt)-0.5)*delta;
                P.r[j]=(pos[n].r[j]+dis);
            }
            P.set_PB(L);
            d=energy_LL_one(P,n)-energy_LL_one(pos[n],n);

            if(d<=0||xsi(mt)<exp(-d/T)){
                pos[n].copy(P);
                LL[pos_cell_1[n]].remove(n);
                pos_cell_1[n]=find_cell(P);
                LL[pos_cell_1[n]].push_front(n);
            }
        }
    };

    void simu_LL(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,pm=0.,pm2=0.;
        f.open(s,ios::out | ios::trunc);
        print_para(f);
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2 "\
        <<"5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        
        create_LL();
        set_g_r(50);
        for(int i=1;i<time+1;i++){
            mc_step_LL();
            if(i>time_start){
                energy_pressure_LL();
                em+=e;
                em2+=e*e;
                pm+=p;
                pm2+=p*p;
                add_g_r();
                if(i/time_take>k){
                    f<<i<<" "<<e<<" "<<em/(i-time_start)<<" "<<\
                    sqrt(em2/(i-time_start)-MyPow(em/(i-time_start),2))<<" "<<\
                    p<<" "<<pm/(i-time_start)<<" "<<\
                    sqrt(pm2/(i-time_start)-MyPow(pm/(i-time_start),2))<<endl;

                    k=i/time_take;
                }
            }
        }
        eval_g_r();
        f.close();
        char sg[100];
        sprintf(sg,"gr_%.2f_%.2f.txt",rho,T);
        f.open(sg,ios::out | ios::trunc);
        f<<"#g(r) using samples:"<<n_count<<" and bins:"<<n_bin<<endl;
        print_para(f);
        for(int i=0;i<g_r.size();i++){
            f<<_r[i]<<" "<<g_r[i]<<endl;
        }
        f.close();
    };
};

class NTP:public NTV{
    public:
    double delta_v;

    void set_value_NTP(double t,double _p,int ti,double d,double dv){
        T=t;
        time=ti;
        delta=d;
        p=_p;
        delta_v=dv;
    };

    void mc_step_volume(){
        double Vo,Lo,_d;
        double _e=e;
        Vo=V;
        V=V*exp((xsi(mt)-0.5)*delta_v);
        Lo=L;
        L=pow(V,1./3.);

        for(int i=0;i<N;i++){
            pos[i]=pos[i]*(L/Lo);
        }
        create_LL();
        energy_LL();
        _d=e-_e;
        _d=-1.*(_d+p*(V-Vo)-(N+1)*log(V/Vo)*T)/T;
        if(xsi(mt)>exp(_d)){
            for(int i=0;i<N;i++){
                pos[i]=pos[i]*(Lo/L);
            }
            L=Lo;
            V=Vo;
            create_LL();
        }
        
    };

    void mc_step_NTP(){
        double r=xsi(mt)*(N+1)+1;
        if(r<N){
            mc_step_LL();
        }
        else{
            mc_step_volume();
        }
    };

    void print_para_NTP(fstream &f){
        f<<"#N:"<<N<<" rho:"<<rho<<\
        " P:"<<p<<" r_c:"<<rc<<" T:"<<\
        T<<" Delta:" <<delta<<" Delta V:"<<delta_v<<endl;
    }
    void simu_LL_NTP(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,rhom=0.,rhom2=0.;
        
        f.open(s,ios::out | ios::trunc);
        print_para_NTP(f);
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2"<< \
        " 5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        
        create_LL();
        for(int i=1;i<time+1;i++){
            mc_step_NTP();
            if(i>time_start){
                energy_LL();
                //cout<<e<<" "<<p<<endl;
                em+=e;
                em2+=e*e;
                rho=N/V;
                rhom+=rho;
                rhom2+=rho*rho;
                if(i/time_take>k){
                    f<<i<<" "<<e<<" "<<em/(i-time_start)<<" "\
                    <<sqrt(em2/(i-time_start)-MyPow(em/(i-time_start),2))<<" "<<\
                    rho<<" "<<rhom/(i-time_start)<<" "<<\
                    sqrt(rhom2/(i-time_start)-MyPow(rhom/(i-time_start),2))<<endl;
                    k=i/time_take;
                }
            }
        }
        f.close();
    };  
};