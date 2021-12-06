/*
Main for simulation.
1: NTV
2: NPT
3: SUS
If defined it is possible to perform also the perfect gas simulation. It will put the cutoff radius equal to -1
*/

//#define ense 2
//#define PERF 

#include "LJ_project.hpp"

int N,t_tot,t_start,t_freq;
double delta;

//functions for printing errors in absence of arguments
void error_NTV(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T filename"<<endl;
    exit(EXIT_FAILURE);
}

void error_NTV_gr(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T bin_size \
    filename g_r_filename"<<endl;
    exit(EXIT_FAILURE);
}

void error_NPT(char*s){
    cerr<<"Error"<<endl<<"Usage:"<<s<<" rho T p filename"<<endl;
    exit(EXIT_FAILURE);
}

//assign function
//data are taken by the file para.txt and must be at least 5
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
    //common error
    if(argc<3){
        cerr<<"Error"<<endl<<"At least 2 parameter"<<endl;
        exit(EXIT_FAILURE);
    }
    double rho=atof(argv[1]);
    double T=atof(argv[2]);
    fstream f;
    f.open("para.txt",ios::in);
    assign(f);
    unsigned long long int seed=time(NULL);
    mt.seed(seed);

    NPT LJ;
    LJ.set_value(rho,N);
    LJ.pos_SC();

    //simulation in NTV
    #if ense==0
    //SImulation evaluating also g_r
    #ifdef gr

    if(argc<6){
        error_NTV_gr(argv[0]);
    }
    LJ.set_value_NTV(T,t_tot,delta);
    int n_bin;
    n_bin=atoi(argv[3]);
    cerr<<"Simulation NTV using Linked List"<<endl\
    <<"Radial distribution functionis evaluated"<<endl;
    LJ.simu(argv[4],argv[5],t_start,t_freq,n_bin);
    cerr<<"Done!"<<endl;
    #else
    if(argc<4){
        error_NTV(argv[0]);
    }
    LJ.set_value_NTV(T,t_tot,delta);
    cerr<<"Simulation NTV using Linked List"<<endl;
    LJ.simu(argv[3],t_start,t_freq);
    cerr<<"Done!"<<endl;
    #endif

    #else
    //NPT simulation: takes pressure and filename. 
    //In the para.txt there must be another 
    //value for volume displacement
    if(argc<5){
        error_NPT(argv[0]);
    }
    double P=atof(argv[3]);
    char s[100];
    f>>s;
    double dv;
    dv=atof(s);
    LJ.set_value_NPT(T,P,t_tot,delta,dv);
    fprintf(stderr,"Simulation NPT using LL algorithm\n");
    LJ.simu_NPT(argv[4],t_start,t_freq);
    #endif
    return 0;
}
