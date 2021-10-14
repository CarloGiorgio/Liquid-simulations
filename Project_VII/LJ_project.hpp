#include "vec.hpp"

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
    vector<list<int> > NN_list;
    vector <int> pos_cell_1;

    void set_value(double _rho,int n){
        rho=_rho;
        N=n;
        V=N/rho;
        L=pow(V,1./3.);
        pos.resize(N);
        e=0.;
        rc=2.5;
        u_shift=LJ(rc);
        tail();
    };

    void set_value(double _rho,int n,double r_c){
        rho=_rho;
        N=n;
        V=N/rho;
        L=pow(V,1./3.);
        pos.resize(N);
        e=0.;
        rc=r_c;
        u_shift=LJ(rc);
        if(rc<0||N<2){
            u_tail=0.;
            p_tail=0.;
        }else{
            tail();
        }
    };

    void tail(){
        double rc3,rc6;
        rc3=1.0/MyPow(rc,3);
        rc6=rc3*rc3;
        u_tail=8./9.*M_PI*rho*rc3*(rc6-3.0);
        p_tail=32.0/9.*M_PI*rho*rho*rc3*(rc6-3./2.);
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

    //Interaction Potential
    double LJ(double r){
        return 4*(1.0/MyPow(r,12)-1.0/MyPow(r,6));
    }

    double potential(MyVec &v,MyVec &w){
        double d;
        MyVec _r;
        _r=v.MIC(w,L);
        d=_r.norm();
        if (d<rc){
            return LJ(d)-u_shift;
        }else{
            return 0.0;
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

/*
Function for Linked list
*/
//function for the cell individuation
    int find_cell(MyVec v){
        int n=0;
        int n_p;
        for(int i=0;i<3;i++){
            n_p=(int)(rint(v.r[i]/rn)+M)%M;
            n+=n_p*MyPow(M,i);
        }
        return n;
    };

    //Linked List creation
    //It's created both the list and the neighbours
    void create_LL(){
        if(rc<0.){
            M=4;
            rn=L/M;
        }else{
            rn=rc+0.2;
            //rn=rc;
            M=rint(L/rn);
        }
        
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

    void create_LL(double _rc){
        rn=_rc;
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
        double _e=0.;
        for(int i : NN_list[pos_cell_1[n]]){
            for(int j:LL[i]){
                if(j!=n){
                    _e+=potential(v,pos[j]);
                }
            }
        }
        return _e;
    };

    double energy_LL_one_new(MyVec v,int n){
        double _e=0.;
        for(int i : NN_list[n]){
            for(int j:LL[i]){
                _e+=potential(v,pos[j]);
                }
            }
        return _e;
    };

    double energy_LL(){
        //find a way to evaluate the neighbours and the 
        //loop over the list 
        double _e=0.;
        for(int l=0;l<N;l++){
            _e+=energy_LL_one(pos[l],l)/2.0;
        }
        return _e;
    };
    
    void update_LL(){
        for(int i=0;i<LL.size();i++){
            LL[i].clear();
        }
        for(int i=0;i<N;i++){
            pos_cell_1[i]=find_cell(pos[i]);
            LL[pos_cell_1[i]].push_front(i);
        }        
    };
    
};

/*
************************************************************************
*/
// Definition of a class for MC simulation in NTV ensemble, Derived by MySet
class NTV:public MySet{

    public:
    double p;
    double T;
    double beta;
    int time;
    double delta;

    void set_value_NTV(double t,int ti,double d){
        T=t;
        beta=1.0/T;
        time=ti;
        delta=d;
        p=0.;
    };

    /*
    void print_para(){
        cout<<"#N:"<<N<<" rho:"<<rho<<" L:"<<L<<" r_c:"<<rc<<" T:"<<T<<" Delta:" <<delta;
    }*/

    void print_para(fstream &f){
        f<<"#N:"<<N<<" rho:"<<N/V<<" r_c:"<<rc<<" T:"<<T<<" Delta:" <<delta<<endl;
    }

    void mean_print_file(fstream&f,int t,int freq,int start,int*k,\
        double o,double a,double *o1,double *o2,double *a1,double*a2){
        (*o1)+=o;
        (*o2)+=o*o;
        (*a1)+=a;
        (*a2)+=a*a;
        if(t/freq>(*k)){
            f<<t<<" "<<o<<" "<<(*o1)/(t-start)<<" "\
            <<sqrt((*o2)/(t-start)-MyPow((*o1)/(t-start),2))<<" "<<\
            a<<" "<<*a1/(t-start)<<" "<<\
            sqrt((*a2)/(t-start)-MyPow((*a1)/(t-start),2))<<endl;
            (*k)+=1;
            f.flush();
        }
    };  
    void mean_print_file(fstream&f,int t,int freq,int start,int*k,\
        double o,double a,double *o1,double *o2,double *a1,double*a2,vector<double> v){
        (*o1)+=o;
        (*o2)+=o*o;
        (*a1)+=a;
        (*a2)+=a*a;
        if(t/freq>(*k)){
            f<<t<<" "<<o<<" "<<(*o1)/(t-start)<<" "\
            <<sqrt((*o2)/(t-start)-MyPow((*o1)/(t-start),2))<<" "<<\
            a<<" "<<*a1/(t-start)<<" "<<\
            sqrt((*a2)/(t-start)-MyPow((*a1)/(t-start),2));
            for(double _p:v){
                f<<" "<<_p;
            }
            f<<endl;
            (*k)+=1;
            f.flush();
        }
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
        int n;
        int nn,no;

        double d;
        for(int i=0;i<N;i++){
            n=floor(xsi(mt)*N);
            d=energy_LL_one(pos[n],n);
            
            //cout<<n<<" I'm in"<<endl;
            P.set_random_2(delta);
            P=pos[n]+P;
            P.set_PB(L);
            nn=find_cell(P);
            no=pos_cell_1[n];

            if(nn!=no){
                LL[no].remove(n);
                pos_cell_1[n]=nn;
                LL[nn].push_front(n);
            }

            d=energy_LL_one(P,n)-d;
            if(d<0||xsi(mt)<exp(-d*beta)){
                pos[n].copy(P);
            }else{
                if(no!=nn){
                    LL[nn].remove(n);
                    pos_cell_1[n]=no;
                    LL[no].push_front(n);
                }
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
        //set_g_r(50);
        for(int i=1;i<time+1;i++){
            mc_step_LL();
            if(i>time_start){
                energy_pressure_LL();
                mean_print_file(f,i,time_take,time_start,\
                &k,e,p,&em,&em2,&pm,&pm2);
                //add_g_r();
            }
        }
        f.close();
    };
};

class NTP:public NTV{

    public:
    double delta_v;
    double G;
    double pv;

    void set_value_NTP(double t,double _p,int ti,double d,double dv){
        T=t;
        beta=1.0/T;
        time=ti;
        delta=d;
        pv=_p;
        delta_v=dv;
        p=0.;
        G=0.0;
    };

    void print_para_NTP(fstream &f){
        f<<"#N:"<<N<<\
        " P:"<<pv<<" r_c:"<<rc<<" T:"<<\
        T<<" Delta:" <<delta<<" Delta V:"<<delta_v<<endl;
    };

    void scale(double s){
        for(int i=0;i<N;i++){
            pos[i]=pos[i]*s;
        }
    };

    double energy_gibbs_LL(){
        return (energy_LL()+pv*V-N*T*log(V));
    }

    void mc_step_volume_LL(){
        double Vo,ran,Lo,_d;
        _d=energy_gibbs_LL();
        Vo=V;
        ran=log(Vo)+(xsi(mt)-0.5)*delta_v;
        V=exp(ran);
        Lo=L;
        L=exp(ran/3.0);
        scale(L/Lo);

        if(rint(L/rn)!=M){
            create_LL();
        }else{
            update_LL();
        }

        _d=energy_gibbs_LL()-_d-log(V/Vo)/beta;
        if(_d>0&&xsi(mt)>exp(-beta*_d)){
            scale(Lo/L);
            L=Lo;
            V=Vo;

            if(rint(L/rn)!=M){
                create_LL();
            }else{
                update_LL();
            }
        }else{
            rho*=Vo/V;
            tail();
        }
    };

    void mc_step_NTP_LL(){
        double r=xsi(mt)*(N+1);
        if(r<N){
            mc_step_LL();
        }
        else{
            mc_step_volume_LL();
            //cerr<<L<<endl;
        }
        
    };


    void simu_LL_NTP(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,rhom=0.,rhom2=0.;
        vector<double> can;
        can.resize(2);
        f.open(s,ios::out | ios::trunc);
        print_para_NTP(f);
        f<<"# 1: Time 2: Gibbs energy 3:<G> 4:<G^2>-<G>^2"<< \
        " 5:density 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        double corr=0.;
        if(rc>0){
            corr=8.0/3.0*M_PI*(1.0/MyPow(rc,9)-1.0/MyPow(rc,3));
        }
        create_LL();
        for(int i=1;i<time+1;i++){
            mc_step_NTP_LL();
            if(i>time_start){
                energy_pressure_LL();
                can[0]=e;
                can[1]=p+corr*rho*rho;
                G=energy_gibbs_LL()/N+u_tail;
                mean_print_file(f,i,time_take,time_start,&k,G,rho,&em,\
                &em2,&rhom,&rhom2,can);
                //cout<<e<<" "<<p<<endl;
            }
        }
        f.close();
    };  
};