#include "vec.hpp"


//Class that contains parameters and generic methods for the simulation


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
    double u_tail;
    double p_tail;

    //radial distribution variables
    vector <double> g_r;
    vector <double> _r;
    int n_bin;
    double dg;
    double n_count;

    //LL variables
    int M;
    vector <list <int> > LL;
    vector<list<int> > NN_list;
    vector <int> pos_cell;

    //initialization variables
    void set_value(double _rho,int n){
        rho=_rho;
        N=n;
        V=(double) N/rho;
        L=pow(V,1./3.);
        pos.resize(N);
        e=0.;
        rc=2.5;
        tail();
    };
    //initialization variables with cutoff radius;
    //used for perfect gas case
    void set_value(double _rho,int n,double r_c){
        rho=_rho;
        N=n;
        V=N/rho;
        L=pow(V,1./3.);
        if(pos.size()!=0){
            pos.clear();
        }
        pos.resize(N);
        e=0.;
        rc=r_c;
        tail();
    };

    //tai correction to the truncated potential
    //The energy is the intensive one

    void tail(){
        if(rc>0||N<2){
            double rc3,rc6;
            rc3=1.0/MyPow(rc,3);
            rc6=rc3*rc3;
            u_tail=8./9.*M_PI*rho*rc3*(rc6-3.0);
            p_tail=32.0/9.*M_PI*rho*rho*rc3*(rc6-3./2.);
        }else{
            u_tail=0.;
            p_tail=0.;
        }
    };

    //rescaling for new density
    void tail(double _rho){
        u_tail=u_tail*_rho/rho;
        p_tail=p_tail*_rho*_rho/rho/rho;
    };


    //************* Begin g(r) block **************//
    //g(r) setting; see Frenkel
    //creation of variables for g(r). It is given the binning
    void set_g_r(int _n){
        n_bin=_n;
        if(g_r.size()!=0){
            g_r.clear();
            _r.clear();
        }
        g_r.resize(n_bin);
        _r.resize(n_bin);
        dg=L/2./n_bin;
        n_count=0.;
    };

    //sampling for g(r)
    //run over all particles and find the ones in the bin
    // It is a O(N^2) algorithm...
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
    //************* End g(r) block **************//

    //********** initial position *********//
    //uniformly in the box
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

    //************** Energy methods ********//
    //Interaction Potential
    double LJ(double r){
        return 4*(1.0/MyPow(r,12)-1.0/MyPow(r,6));
    }

    //evaulation of 2 particles contribution
    double potential(MyVec &v,MyVec &w){
        double d;
        MyVec _r;
        _r=v.MIC(w,L);
        d=_r.norm();
        if (d<rc){
            return LJ(d);
            //return LJ(d)-u_shift;
        }else{
            return 0.0;
        }
    };

    //force between two particles
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

    //********************* Begin Linked list method *******************//
    //function for cell individuation
    int find_cell(MyVec v){
        int n=0;
        int n_p;
        for(int i=0;i<3;i++){
            n_p=(int)(rint(v.r[i]/rc)+M)%M;
            n+=n_p*MyPow(M,i);
        }
        return n;
    };

    //Linked List creation
    //It's created both the list and the neighbours

    void create_LL(){
        if(rc<0.){
            M=1;
        }else{
            M=rint(L/rc);
        }
        //forcing an algorithm brute force
        if(M<3){
            M=1;
        }
        //when called the second time clear everything!
        if(LL.size()!=0){
            for(int i=0;i<LL.size();i++){
                LL[i].clear();
                NN_list[i].clear();
            }
            LL.clear();
            NN_list.clear();
            pos_cell.clear();
        }
        
        NN_list.resize(M*M*M);
        LL.resize(M*M*M);
        pos_cell.resize(N);

        for(int i=0;i<N;i++){
            pos_cell[i]=find_cell(pos[i]);
            LL[pos_cell[i]].push_front(i);
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

    //update of the list
    //used after a volume step
    void update_LL(){
        for(int i=0;i<LL.size();i++){
            LL[i].clear();
        }
        for(int i=0;i<pos.size();i++){
            pos_cell[i]=find_cell(pos[i]);
            LL[pos_cell[i]].push_front(i);
        }
    };


    //evaluation of energy when there is a translation.
    //It avoids an infinite contribution
    double energy_one(MyVec v,int n){
        double _e=0.;
        if (rc<0){
            return _e;
        }

        //cycle over neighbours list
        for(int i : NN_list[pos_cell[n]]){
            //cycle over particles in list
            for(int j:LL[i]){
                if(j!=n){
                    _e+=potential(v,pos[j]);
                }
            }
        }
        return _e;
    };

    //evaluation of energy when a new particle is inserted
    //same as befor but we cuse all the particles in the list

    double energy_one_new(MyVec v,int n){
        double _e=0.;
        if (rc<0) {
            return _e;
        }
        for(int i : NN_list[n]){
            for(int j:LL[i]){
                _e+=potential(v,pos[j]);
                }
            }
        return _e;
    };

    //total energy evaluation. It is extensive. NOT correction to the tail
    double energy(){
        double _e=0.;
        if(rc<0.) {
            return _e;
        }
        for(int l=0;l<N;l++){
            _e+=energy_one(pos[l],l)/2.0;
        }
        return _e;
    };

    
    
    //********************* End Linked list method *******************//
    /*
    void create(double _rc){
        rn=_rc;
        //rn=rc;
        M=rint(L/rn);
        if(M<3){
            M=3;
            rn=(double) L/M;
        }
        NN_list.resize(M*M*M);
        LL.resize(M*M*M);
        pos_cell.resize(N);

        for(int i=0;i<N;i++){
            pos_cell[i]=find_cell(pos[i]);
            LL[pos_cell[i]].push_front(i);
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
*/

};

/*
******************************
*                            *
*  Class for NTV simulation  *
*                            *
******************************
*/
class NTV:public MySet{

    public:
    //variable with thermodynamic quantities : pressure, temperature, beta, simulation time and displacement
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

    //print parameters on a file
    void print_para(fstream &f){
        f<<"#N:"<<N<<" rho:"<<(double)N/V<<" r_c:"<<rc<<" T:"<<T<<" Delta:" <<delta<<endl;
    }

    //update and print with a given frequency the value, the mean and the fluctuations of the two observables
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
            (*k)=t/freq;
            f.flush();
        }
    };
    void mean_print_file(fstream&f,int t,int freq,int start,int*k,\
        double o,double a,double *o1,double *o2,\
        double *a1,double*a2,vector<double> v){
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
            (*k)=t/freq;
            f.flush();
        }
    };

    //evaluation of the energy and the pressure using linked list
    //for the pressure we use the virial
    void energy_pressure(){
        double d,d12,d6;
        p=p_tail+rho*T;
        e=u_tail;
        if(rc<0) return;
        for(int i=0;i<N;i++){
            for(int l : NN_list[pos_cell[i]]){
                for(int j:LL[l]){
                    if(i!=j){
                        d=pos[i].dist(pos[j],L);
                        if(d<rc){
                            d6=MyPow(d,6);
                            d6=1./d6;
                            d12=d6*d6;
                            p+=48.*(d12-0.5*d6)/V/6.;
                            e+=4.*(d12-d6)/2./N;
                            //e+=(4.*(d12-d6)-u_shift)/2./N;
                        }
                    }
                }
            }
        }
    };

    //single montecarlo translation using linked list
    void mc_step_single(){
        MyVec P;
        int n;
        int nn,no;
        double d;
        n=(int)(xsi(mt)*N)%N;
        d=energy_one(pos[n],n);

        //displacement of a particle in uniform box of delta spacing
        P.set_random_2(delta);
        P=pos[n]+P;
        P.set_PB(L);
        nn=find_cell(P);
        no=pos_cell[n];

        //in the case the particle crossed a boundary: new index
        if(nn!=no){
            LL[no].remove(n);
            pos_cell[n]=nn;
            LL[nn].push_front(n);
        }
        d=energy_one(P,n)-d;

        //montecarlo acceptance rule
        if(d<0||xsi(mt)<exp(-d*beta)){
            pos[n].copy(P);
        }else{
            if(no!=nn){
                LL[nn].remove(n);
                pos_cell[n]=no;
                LL[no].push_front(n);
            }
        }
    };

    void mc_step(){
        MyVec P;
        int n;
        int nn,no;

        double d;
        for(int i=0;i<N;i++){
            n=(int)(xsi(mt)*N)%N;
            d=energy_one(pos[n],n);

            //displacement of a particle in uniform box of delta spacing
            P.set_random_2(delta);
            P=pos[n]+P;
            P.set_PB(L);
            nn=find_cell(P);
            no=pos_cell[n];

            //in the case the particle crossed a boundary: new index
            if(nn!=no){
                LL[no].remove(n);
                pos_cell[n]=nn;
                LL[nn].push_front(n);
            }
            d=energy_one(P,n)-d;

            //montecarlo acceptance rule
            if(d<0||xsi(mt)<exp(-d*beta)){
                pos[n].copy(P);
            }else{
                if(no!=nn){
                    LL[nn].remove(n);
                    pos_cell[n]=no;
                    LL[no].push_front(n);
                }
            }
        }
    };

    //method for simulation
    void simu(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,pm=0.,pm2=0.;
        f.open(s,ios::out | ios::trunc);
        print_para(f);
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2 "\
        <<"5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        create_LL();
        for(int i=1;i<time+1;i++){
            mc_step();
            if(i>time_start){
                energy_pressure();
                mean_print_file(f,i,time_take,time_start,\
                &k,e,p,&em,&em2,&pm,&pm2);
            }
        }
        f.close();
    };

    //simulation where also the g(r) is evaluated
    void simu(string s,string gs,int time_start,\
     int time_take,int bin_size){
        fstream f;
        fstream fi;
        double em=0.,em2=0.,pm=0.,pm2=0.;
        f.open(s,ios::out | ios::trunc);
        fi.open(gs,ios::out | ios::trunc);
        print_para(f);
        print_para(fi);
        //fi<<"#Evaluation of g(r)"<<endl<<"# 1:r 2:g(r)"<<endl;
        f<<"# 1: Time 2: energy 3:<e> 4:<e^2>-<e>^2 "\
        <<"5:pressure 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        create_LL();
        set_g_r(bin_size);
        for(int i=1;i<time+1;i++){
            mc_step();
            if(i>time_start){
                energy_pressure();
                if(i/time_take>k) add_g_r();
                mean_print_file(f,i,time_take,time_start,\
                &k,e,p,&em,&em2,&pm,&pm2);
                
            }
        }
        f.close();
        eval_g_r();
        for(int i=0; i<g_r.size();i++){
            fi<<_r[i]<<" "<<g_r[i]<<endl;
        }
        fi.flush();
        fi.close();
    };
};


/*
******************************
*                            *
*  Class for NPT simulation  *
*                            *
******************************
*/
class NPT:public NTV{

    public:

    //varibale for NPT: volume displacement, enthalpy and pressure (true)

    double delta_v;
    double G;
    double pv;

    void set_value_NPT(double t,double _p,int ti,double d,double dv){
        T=t;
        beta=1.0/T;
        time=ti;
        delta=d;
        pv=_p;
        delta_v=dv;
        p=0.;
        G=0.0;
    };

    void print_para_NPT(fstream &f){
        f<<"#N:"<<N<<\
        " P:"<<pv<<" r_c:"<<rc<<" T:"<<\
        T<<" Delta:" <<delta<<" Delta V:"<<delta_v<<endl;
    };

    //scaling of particles after a volume move
    void scale(double s){
        for(int i=0;i<N;i++){
            pos[i]=pos[i]*s;
        }
    };

    //Gibbs free energy
    //inserted also tail correction
    double energy_gibbs(){
        return (energy()+pv*V-N*T*log(V)+N*u_tail);
    }

    //Montecarlo volume step
    // a random walk in the log of the volume is performed

    void mc_step_volume(){
        int _n;
        double Vo,ran,Lo,_d;
        _d=energy_gibbs();
        Vo=V;
        ran=log(Vo)+(xsi(mt)-0.5)*delta_v;
        V=exp(ran);
        Lo=L;
        L=exp(ran/3.0);
        tail(N/V);
        scale(L/Lo);
        rho*=Vo/V;

        //updating list, in case creating a new one
        _n=rint(L/rc);
        if(_n==M||(_n<3&&M<3)){
            update_LL();
        }else{
            create_LL();
        }

        //must be considered the tail correction!
        _d=energy_gibbs()-_d-log(V/Vo)*T;

        //acceptance rule: actually moved accepted and seeing if must be rejected
        if(_d>0&&xsi(mt)>exp(-beta*_d)){
            scale(Lo/L);
            rho*=V/Vo;
            L=Lo;
            V=Vo;

            tail(N/V);
            _n=rint(L/rc);
            if(_n==M||(_n<3&&M<3)){
                update_LL();
            }else{
                create_LL();
            }

        }
    };

    //proposal move: 1/N probability of volume move
    void mc_step_NPT(){
        double r=xsi(mt)*(N+1)+1;
        if(r<=N){
            mc_step();
        }
        else{
            mc_step_volume();
        }

    };


    void simu_NPT(string s,int time_start, int time_take){
        fstream f;
        double em=0.,em2=0.,rhom=0.,rhom2=0.;

        //stores energy and pressure from virial
        vector<double> can;
        can.resize(2);

        f.open(s,ios::out | ios::trunc);
        print_para_NPT(f);
        f<<"# 1: Time 2: Gibbs energy 3:<G> 4:<G^2>-<G>^2"<< \
        " 5:density 6:<p> 7:<p^2>-<p>^2"<<endl;
        int k=0;
        create_LL();
        for(int i=1;i<time+1;i++){
            mc_step_NPT();
            if(i>time_start){
                energy_pressure();
                can[0]=e;
                can[1]=p;
                G=energy_gibbs()/N;
                mean_print_file(f,i,time_take,time_start,&k,G,rho,&em,\
                &em2,&rhom,&rhom2,can);
                //cout<<e<<" "<<p<<endl;
            }
        }
        f.close();
    };
};

