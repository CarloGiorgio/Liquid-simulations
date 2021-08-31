#define ense 2
#define BRUTE 0
#define PERF 
#include "simu.hpp"

int N,t_tot,t_start,t_freq;
double delta;
void error_NTV(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T filename"<<endl;
    exit(EXIT_FAILURE);
}
void error_NPT(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T p filename"<<endl;
    exit(EXIT_FAILURE);
}
void error_muTV(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T zeta filename"<<endl;
    exit(EXIT_FAILURE);
}

void error_muTV_SUS(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T zeta data_filename\
    hist_head_filename"<<endl;
    exit(EXIT_FAILURE);
}
void assign(fstream &f){
    char s[100];
    f>>s;
    N=atoi(s);
    f>>s;
    delta=atof(s);
    f>>s;
    t_tot=atoi(s);
    f>>s;
    t_start=atoi(s);
    f>>s;
    t_freq=atoi(s);
}

int main(int argc,char* argv[]){
    if(argc<3){
        cerr<<"Error"<<endl<<"At least 2 parameter"<<endl;
        exit(EXIT_FAILURE);
    }
    double rho=atof(argv[1]);
    double T=atof(argv[2]);
    fstream f;
    f.open("para.txt",ios::in);
    assign(f);
    mt.seed(time(NULL));

    muTV LJ;
    #ifdef PERF
    LJ.set_value(rho,N,-1);
    fprintf(stderr,"Perfect Gas simulation\n");
    #else
    LJ.set_value(rho,N);
    #endif
    //LJ.set_value_NTV(T,300,0.05);
    LJ.pos_BCC();
    //LJ.create_LL();

    #if ense==0
    if(argc<4){
        error_NTV(argv[0]);
    }
    LJ.set_value_NTV(T,t_tot,delta);

    //brute force algorithm
    #if BRUTE
    fprintf(stderr,"Simulation NTV using brute algorithm\n");
    LJ.simu_brute(argv[3],t_start,t_freq);

    #else
    fprintf(stderr,"Simulation NTV using Linked List\n");
    LJ.simu_LL(argv[3],t_start,t_freq);
    #endif

    #elif ense==1
    if(argc<5){
        error_NPT(argv[0]);
    }
    double P=atof(argv[3]);
    char s[100];
    f>>s;
    double dv;
    dv=atof(s);
    LJ.set_value_NTP(T,P,t_tot,delta,dv);
    #if BRUTE
    fprintf(stderr,"Simulation NTP using brute algorithm\n");
    LJ.simu_brute_NTP(argv[4],t_start,t_freq);
    #else
    fprintf(stderr,"Simulation NTP using LL algorithm\n");
    LJ.simu_LL_NTP(argv[4],t_start,t_freq);
    #endif

    #elif ense==2
    if(argc<5){
        error_muTV(argv[0]);
    }   
    double z=atof(argv[3]);
    char s[100];
    f>>s;
    int _n=atoi(s);
    LJ.set_value_muTV(T,z,t_tot,delta,_n);
    #if BRUTE
    fprintf(stderr,"Simulation muTV using brute algorithm\n");
    LJ.simu_brute_muTV(argv[4],t_start,t_freq);
    #else
    fprintf(stderr,"Simulation muTV using LL algorithm\n");
    LJ.simu_LL_muTV(argv[4],t_start,t_freq);
    #endif    
    #else
    if(argc<6){
        error_muTV_SUS(argv[0]);
    }   
    double z=atof(argv[3]);
    char s[100];
    f>>s;
    int _n=atoi(s);
    LJ.set_value_muTV(T,z,t_tot,delta,_n);
    _n=1;
    t_tot=rint(0.3*LJ.V);
    string _s;
    #if BRUTE
    fprintf(stderr,"Simulation muTV using brute algorithm\n\
    Successive Umbrella Sampling");
    while(_n<t_tot){
        _s=argv[5]+"_"+_n+".txt";
        LJ.set_value(_n/LJ.V,_n);
        LJ.simu_brute_muTV_SUS(argv[4],t_start,t_freq,_s,_n);
        _n+=1;
    }
    #else
    fprintf(stderr,"Simulation muTV using LL algorithm\n\
    Successive Umbrella Sampling");
    while(_n<t_tot){
        _s=argv[5]+"_"+_n+".txt";
        LJ.simu_LL_muTV_SUS(argv[4],t_start,t_freq,_s,_n);
        _n+=1;
    }
    #endif

    return 0;
}