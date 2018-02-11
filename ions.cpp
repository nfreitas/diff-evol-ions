#include <stdio.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>

#include <vector>
#include <algorithm>

#include <random>
#include <time.h>

#include "diffevol.h"

using namespace std;

typedef vector<double> VECTOR;

double energy(VECTOR* X, VECTOR* Y, VECTOR* Z, double alpha_x, double alpha_y){

    double E_x = 0, E_y = 0, E_z = 0, E_R = 0;

    int N = X->size();

    for(int i=0; i < N; i++){
        
        double x_i = (*X)[i];
        double y_i = (*Y)[i];
        double z_i = (*Z)[i];

        E_x += pow(x_i,2);
        E_y += pow(y_i,2);
        E_z += pow(z_i,2);

        for(int j=i+1; j < N; j++){

            double x_j = (*X)[j];
            double y_j = (*Y)[j];
            double z_j = (*Z)[j];

            double R = pow( pow(x_i-x_j,2)+pow(y_i-y_j,2)+pow(z_i-z_j,2) , -.5);
            E_R += R;
        }
    }
    return pow(alpha_x,2)*E_x + pow(alpha_y,2)*E_y + E_z + 2*E_R;
}

void gradient(VECTOR* G, VECTOR* X, VECTOR* Y, VECTOR* Z, double alpha_x, double alpha_y){
    
    int N = X->size();
    
    G->resize(3*N);

    for(int i=0; i < N; i++){
        
        double x_i = (*X)[i];
        double y_i = (*Y)[i];
        double z_i = (*Z)[i];

        double S_x = 0, S_y = 0, S_z = 0;

        for(int j=0; j < N; j++){

            if(j==i) continue;

            double x_j = (*X)[j];
            double y_j = (*Y)[j];
            double z_j = (*Z)[j];

            double IR3 = pow( pow(x_i-x_j,2)+pow(y_i-y_j,2)+pow(z_i-z_j,2) , -1.5);
            S_x += (x_j-x_i)*IR3;
            S_y += (y_j-y_i)*IR3;
            S_z += (z_j-z_i)*IR3;
        }

        (*G)[i] = pow(alpha_x,2)*x_i + S_x;
        (*G)[N+i] = pow(alpha_y,2)*y_i + S_y;
        (*G)[2*N+i] = z_i + S_z;

    }

}

namespace diffevol
{

    class CustomObjetive: public Objetive{
        public:
            CustomObjetive(){}

            virtual void extract_coords(Individual& ind, VECTOR& X, VECTOR& Y, VECTOR& Z){};
            virtual Individual eval_gradient(Individual& ind){};

    };

    class Crystal: public CustomObjetive{
    
        public:

            double alpha_x;
            double alpha_y;
            
            Crystal(){}

            Crystal(double alpha_x, double alpha_y): CustomObjetive(){
                this->alpha_x = alpha_x; 
                this->alpha_y = alpha_y; 
            }

            virtual void extract_coords(Individual& ind, VECTOR& X, VECTOR& Y, VECTOR& Z){

                int N = ind.size();
                int Ni = N/3;

                X.assign(ind.begin(), ind.begin()+Ni);
                Y.assign(ind.begin()+Ni, ind.begin()+2*Ni);
                Z.assign(ind.begin()+2*Ni, ind.begin()+3*Ni);

            }
    
            virtual double eval(Individual& ind){
            
                VECTOR X,Y,Z;
            
                extract_coords(ind, X, Y, Z);

                return energy(&X,&Y,&Z,alpha_x,alpha_y);
            
            }

            virtual Individual eval_gradient(Individual& ind){

                VECTOR X,Y,Z;
 
                extract_coords(ind, X, Y, Z);

                VECTOR grad;
                gradient(&grad,&X,&Y,&Z,alpha_x,alpha_y);

                Individual corr(grad);

                return corr;
            }
    
    };


    class NoZOrderCrystal: public CustomObjetive{
    
        public:

            double alpha_x;
            double alpha_y;
            
            NoZOrderCrystal(){}

            NoZOrderCrystal(double alpha_x, double alpha_y): CustomObjetive(){
                this->alpha_x = alpha_x; 
                this->alpha_y = alpha_y; 
            }

            virtual void extract_coords(Individual& ind, VECTOR& X, VECTOR& Y, VECTOR& Z){

                int N = ind.size();
                int Ni = N/3;

                X.assign(ind.begin(), ind.begin()+Ni);
                Y.assign(ind.begin()+Ni, ind.begin()+2*Ni);
                Z.assign(ind.begin()+2*Ni, ind.begin()+3*Ni-1);
                Z.push_back(ind.at(2*Ni));

            }
    
            virtual double eval(Individual& ind){
            
                VECTOR X,Y,Z;
            
                extract_coords(ind, X, Y, Z);

                return energy(&X,&Y,&Z,alpha_x,alpha_y);
            
            }

            virtual Individual eval_gradient(Individual& ind){

                VECTOR X,Y,Z;
 
                extract_coords(ind, X, Y, Z);

                VECTOR grad;
                gradient(&grad,&X,&Y,&Z,alpha_x,alpha_y);

                Individual corr(grad);

                return corr;
            }
    
    };


    class SymCrystal: public CustomObjetive{
    
        private:
    
            double alpha_x;
            double alpha_y;
            int *syms;
            int dim;
    
        public:
            SymCrystal(){}
            
            SymCrystal(double alpha_x, double alpha_y, int *syms, int dim): CustomObjetive(){
                this->alpha_x = alpha_x; 
                this->alpha_y = alpha_y; 
                this->syms = syms; 
                this->dim = dim; 
            }

            virtual void extract_coords(Individual& ind, VECTOR& X, VECTOR& Y, VECTOR& Z){

                int N = ind.size();
                int Ni = N/this->dim;

                if(this->dim==3)
                {
                    X.assign(ind.begin(), ind.begin()+Ni);
                    Y.assign(ind.begin()+Ni, ind.begin()+2*Ni);
                    Z.assign(ind.begin()+2*Ni, ind.begin()+3*Ni);
    
                    for(int i=0; i < Ni; i++){
                        X.insert( X.end(), syms[0]*(X[i]) );
                        Y.insert( Y.end(), syms[1]*(Y[i]) );
                        Z.insert( Z.end(), syms[2]*(Z[i]) );
                    }
                }
                else if(this->dim==2)
                {

                    X.assign(ind.begin(), ind.begin()+Ni);
                    Y.assign(Ni,0);
                    Z.assign(ind.begin()+Ni, ind.begin()+2*Ni);
    
                    for(int i=0; i < Ni; i++){
                        X.insert( X.end(), syms[0]*(X[i]) );
                        Y.insert( Y.end(), syms[1]*(Y[i]) );
                        Z.insert( Z.end(), syms[2]*(Z[i]) );
                    }
                
                }
                else if(this->dim==1)
                {

                    X.assign(Ni,0);
                    Y.assign(Ni,0);
                    Z.assign(ind.begin(), ind.begin()+Ni);
    
                    for(int i=0; i < Ni; i++){
                        X.insert( X.end(), syms[0]*(X[i]) );
                        Y.insert( Y.end(), syms[1]*(Y[i]) );
                        Z.insert( Z.end(), syms[2]*(Z[i]) );
                    }
                
                }


            }
    
            virtual double eval(Individual& ind){
            
                VECTOR X,Y,Z;
            
                extract_coords(ind, X, Y, Z);

                return energy(&X,&Y,&Z,alpha_x,alpha_y);
            
            }

            virtual Individual eval_gradient(Individual& ind){

                VECTOR X,Y,Z;
            
                extract_coords(ind, X, Y, Z);
 
                VECTOR total_grad;
                gradient(&total_grad,&X,&Y,&Z,alpha_x,alpha_y);

                int N = ind.size();
                int Ni = N/this->dim; 
                
                VECTOR sym_grad(N,0);

                if(this->dim==3){
                    for(int i=0; i < Ni; i++){
                        sym_grad[i] = total_grad[i] + syms[0]*total_grad[Ni+i];
                        sym_grad[Ni+i] = total_grad[(2*Ni)+i] + syms[1]*total_grad[(3*Ni)+i];
                        sym_grad[2*Ni+i] = total_grad[(4*Ni)+i] + syms[2]*total_grad[(5*Ni)+i];
                    }
                }else if(this->dim==2)
                {
                    for(int i=0; i < Ni; i++){
                        sym_grad[i] = total_grad[i] + syms[0]*total_grad[Ni+i];
                        sym_grad[Ni+i] = total_grad[(4*Ni)+i] + syms[2]*total_grad[(5*Ni)+i];
                    }
                }else if(this->dim==1)
                {
                    for(int i=0; i < Ni; i++){
                        sym_grad[i] = total_grad[(4*Ni)+i] + syms[2]*total_grad[(5*Ni)+i];
                    }
                }


                Individual corr(sym_grad);

                return corr;
            }
    
    };

class DE_custom: public DE
{

    public:

        double gd_step;
        bool gd_last_success;
        CustomObjetive *objetive;

        DE_custom(int D,int NP,Ranges ranges,CustomObjetive* obj,double gd_step_start=.1,RNG rng=RNG())
        : DE(D, NP, ranges, obj, rng){
            gd_step = gd_step_start;
            objetive =obj;
        }
    
        void evolve(){
            
            int target_index = get_target_index();
            Individual target = population[target_index];
            Individual donor = get_donor(target_index);
            Individual trial = crossover(target, donor);
    
            trial.obj_value = objetive->eval(trial);
    
            if (trial <= target) {
                population.erase(population.begin() + target_index);
                population.add(trial);
            }else{
                trial = target - gd_step*(objetive->eval_gradient(target));
                trial.obj_value = objetive->eval(trial);

                if ( trial <= target ){
                    gd_success();
                    population.erase(population.begin() + target_index);
                    population.add(trial);
                }else{
                    gd_fail();
                }
            }

        }

        void gd_success(){
            gd_last_success = true;
            gd_step *= 1.0001;
        }

        void gd_fail(){
            gd_last_success = false;
            gd_step *= 0.999;
        }

        Individual get_donor(int target_index){
    
            // choose two mutually exclusive random integers different from
            // target_index
            int index1 = rng.randint(get_NP()-1);
            int index2 = rng.randint(get_NP()-2);
    
            if(index1 >= target_index) index1++;
    
            int min1 = std::min(target_index,index1);
            int max1 = std::max(target_index,index1);
    
            if(index2 >= min1){
                index2++;
                if(index2 >= max1) index2++;
            }
    
            double F = .5+.5*rng.rand(1);
            return (population.best_individual()+F*(population[index1]-population[index2]));

        }


        Individual crossover(Individual& target, Individual& donor){
    
            Individual trial = target;
    
            double Cr = .5+.5*rng.rand(1);
            
            for(size_t c=0; c<target.size(); c++){
                if (rng.rand(1) <= Cr) trial[c] = donor[c];
            }
    
            return trial;
        }


};
}


int main(int argc, char **argv)
{

    int code;
    int Ni=-1;
    float alpha_x=-1;
    float alpha_y=-1;

    int D=-1,dim=0;
    int NP = 40;
    char symm='i';
    int syms[] = {-1,-1,-1};

    double steady_tresh = 1e-10;
    int max_steps = 5000000;
    int max_steady_steps = 100;

    bool progress = false;

    diffevol::Ranges ranges;
    diffevol::Range range (-10,10);

    diffevol::CustomObjetive *obj;

    while((code=getopt(argc,argv, "N:x:y:I:s:d:M:m:t:p")) != -1){
        switch(code){
            case 'N':
                Ni = atoi(optarg);
                break;
            case 'x':
                alpha_x = atof(optarg);
                break;
            case 'y':
                alpha_y = atof(optarg);
                break;
            case 'I':
                NP = atoi(optarg);
                break;
            case 's':
                symm = optarg[0];
                break;
            case 'd':
                dim = atoi(optarg);
                break;
            case 'M':
                max_steps = atoi(optarg);
                break;
            case 'm':
                max_steady_steps = atoi(optarg);
                break;
            case 't':
                steady_tresh = atof(optarg);
                break;
            case 'p':
                progress = true;
        }
    }


    switch(symm){

        case 'i':
            D = dim*Ni/2;
            obj = new diffevol::SymCrystal(alpha_x,alpha_y,syms,dim);
            break;

        case 'x':
            D = dim*Ni/2;
            syms[0] = 1;
            obj = new diffevol::SymCrystal(alpha_x,alpha_y,syms,dim);
            break;

        case 'y':
            D = dim*Ni/2;
            syms[1] = 1;
            obj = new diffevol::SymCrystal(alpha_x,alpha_y,syms,dim);
            break;
        
        case 'z':
            D = dim*Ni/2;
            syms[0] = 1;
            syms[1] = 1;
            obj = new diffevol::SymCrystal(alpha_x,alpha_y,syms,dim);
            break;

        case 'n':
            D = 3*Ni;
            for(int i=0; i<3; i++) syms[i] = 0;
            obj = new diffevol::Crystal(alpha_x,alpha_y);
            break;

        case '1':
            D = 3*Ni;
            for(int i=0; i<3; i++) syms[i] = 0;
            obj = new diffevol::NoZOrderCrystal(alpha_x,alpha_y);
            break;
    
    }

    if(!(Ni>0)){ cerr << "Invalid numbers of ions.\n"; return 1;}
//    if(Ni%2!=0){ cerr << "Number of ions must be pair\n"; return 1;}
    if(!(alpha_x>0)){ cerr << "Invalid harmonic constant for x axis.\n"; return 1;}
    if(!(alpha_y>0)){ cerr << "Invalid harmonic constant for y axis.\n"; return 1;}
    if(!(NP>4)){ cerr << "Invalid number of individuals.\n"; return 1;}
    if(!(dim<=3 && dim>=1)){ cerr << "Invalid dimension (" << dim << ").\n"; return 1;}


    cout << "# Number of ions: " << Ni << '\n';
    cout << "# Crystal spacial dimensions: " << dim << '\n';
    cout << "# alpha_x: " << alpha_x << '\n';
    cout << "# alpha_y: " << alpha_y << '\n';
//    cout << "# symmetries: x=" << syms[0] << " y=" << syms[1] << " z=" << syms[2] << '\n';
    cout << "# symmetries: " << symm << '\n';
    cout << "# DE - Dimensionality: " << D << '\n';
    cout << "# DE - Number of individuals: " << NP << '\n';
    cout << "# DE - Maximum steps: " << max_steps << '\n';
    cout << "# DE - Maximum steady steps: " << max_steady_steps << '\n';
    cout << "# DE - Steady treshold: " << steady_tresh << '\n';

    for(int i=0; i<D; i++){
        ranges.push_back(range);
    }

    diffevol::DE_custom de(D,NP,ranges,obj);

    int count = 0;
    double last_best;

    int iter;
    for(iter=0; iter < max_steps; iter++){
         if(iter%100 == 0){
            
            diffevol::Population pop = de.get_population();
            diffevol::Individual best = pop.best_individual();
            diffevol::Individual worst = pop.worst_individual();

            double rel_diff = (last_best-best.obj_value)/best.obj_value;
            count = (rel_diff <= steady_tresh) ? count+1:0;
            if(count > max_steady_steps) break;

            if(progress){
            double var = (worst.obj_value - best.obj_value)/best.obj_value;
                cerr << "steps: " << iter;
                cerr << "  best: " << best.obj_value;
                cerr << "  var: " << var;
                cerr << "  gd: " << de.gd_step <<"("<<(de.gd_last_success ? "Ok!" : "Fail!") << ")";
                cerr << "  rel_diff:" << rel_diff;
                cerr << "  st_count:" << count << '\n';
            }

            last_best = best.obj_value;
        }

        de.evolve();
    }

    cout.precision(15);

    diffevol::Individual best = de.get_population().best_individual();

    cout << "# DE - Best value: " << best.obj_value << '\n';
    cout << "# DE - Steps: " << iter << '\n';

    VECTOR X,Y,Z;
    obj->extract_coords(best,X,Y,Z);

    VECTOR sort_X(X.size());
    VECTOR sort_Y(Y.size());
    VECTOR sort_Z(Z.size());

    vector<int> ind(X.size());
    for(int c=0;c<X.size();c++) ind[c] = c;

    sort(ind.begin(),ind.end(),[Z](int i, int j){return Z[i]<Z[j];});
    for(int i=0;i<ind.size();i++){
        sort_X[i] = X[ind[i]];
        sort_Y[i] = Y[ind[i]];
        sort_Z[i] = Z[ind[i]];
    }

    X = sort_X;
    Y = sort_Y;
    Z = sort_Z;

    VECTOR grad;
    gradient(&grad,&X,&Y,&Z,alpha_x,alpha_y);

    for(int k=0; k<X.size(); k++){
        cout << X[k] << "  " << Y[k] << "  " << Z[k];
        cout << " | " << grad[k] << "  " << grad[Ni+k] << "  " << grad[2*Ni+k] ;
        cout << '\n';
    }

    delete(obj);
}

