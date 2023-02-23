#include <vec.hpp>

class set_equation {
    public:
    //configuration variables
    int N;
    double rho;
    double V;
    double L;
    vector<MyVec> pos;
    vector<MyVec> vel;
    vector<MyVec> forces;


    void set_value(double _rho,int n){
        rho = _rho;
        N = n;
        V = (double) N/rho;
        L = pow(V,1./3.);
        pos.resize(N);
        vel.resize(N);
        forces.resize(N);
    };


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

    void vel_random(){
        for(int i=0;i<pos.size();i++){
            vel[i].set_random(L);
        }
    };

    MyVec force(MyVec &vij){
        //insert evaluation of the forces...
        //TODO Use r-12 and hydrodynamic forces with fluid of bacteria...
        double dij;
        dij = vij.norm();
        return (12./MyPow(dij,12))*vij

        return d;
    }

    MyVec force_vel(MyVec &vij,MyVec &r){
        //TODO look better in literature if this is working or not
        //forces acting on the velocity. 
        //Since low reynolds number, it is directly the velocity field
        MyVec w;
        double d;
        w = r.copy();
        d = r.norm();
        return 1/(6.*M_PI/d)*(vij + (w.project(vij)/(d*d)));
    };

    
    int find_cell(MyVec v){
        //function for cell individuation
        int n=0;
        int n_p;
        for(int i=0;i<3;i++){
            n_p=(int)(rint(v.r[i]/rc)+M)%M;
            n+=n_p*MyPow(M,i);
        }
        return n;
    };


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

    
    void update_LL(){
        //update of the list
        for(int i=0;i<LL.size();i++){
            LL[i].clear();
        }
        for(int i=0;i<pos.size();i++){
            pos_cell[i]=find_cell(pos[i]);
            LL[pos_cell[i]].push_front(i);
        }
    };
};


class active_matter {

    public:
    // TODO implement integration timestep; evaluation of forces and other stuffs
    double dt;
    double T;

    void set_time(double dt_,double T_){
        dt = dt_;
        T = T_
    };

    void evaluate_forces(){
        //evaluation of the forces
        MyVec d;
        for(int i=0;i<N;i++){
            for(int l : NN_list[pos_cell[i]]){
                for(int j:LL[l]){
                    if(i!=j){
                        d = pos[i] - pos[j];
                        f = force(d);
                        forces[i] += f;
                        forces[j] += -f;
                    }
                }
            }
        }
    };
    
    void perform_timestep(){
        //Euler method integration
        for (int i = 0;i < N; i++){
            pos[i] += vel[i]*dt;
            vel[i] += (forces[i] + force_vel(vel[i]))*dt;
        }
    };


};

