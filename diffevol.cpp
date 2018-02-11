#include <vector>
#include <cstdlib>
#include <random>
#include <time.h>

#include "diffevol.h"

namespace diffevol 
{

    // Class RNG

    RNG::RNG(){ std::srand(time(NULL)); }

    int RNG::randint(int N){
        std::uniform_int_distribution<> uniform(0,N-1);
        return uniform(rng);
    }
    
    double RNG::rand(int R){
        std::uniform_real_distribution<> uniform(0,R);
        return uniform(rng);
    }

    // Class Individual
    
    bool operator<(const Individual &ind1, const Individual &ind2){
        return (ind1.obj_value < ind2.obj_value);
    }
    bool operator<=(const Individual &ind1, const Individual &ind2){
        return (ind1.obj_value <= ind2.obj_value);
    }
    Individual operator+(const Individual &ind1, const Individual &ind2){
        Individual res;
        for(int c = 0; c < ind1.size(); c++) res.push_back(ind1[c]+ind2[c]);
        return res;
    }
    Individual operator-(const Individual &ind1, const Individual &ind2){
        Individual res;
        for(int c = 0; c < ind1.size(); c++) res.push_back(ind1[c]-ind2[c]);
        return res;
    }
    Individual operator*(const Individual &ind, const double d){
        Individual res;
        for(size_t c = 0; c < ind.size(); c++) res.push_back(d*ind[c]);
        return res;
    }
    Individual operator*(const double d, const Individual &ind){
        Individual res;
        for(size_t c = 0; c < ind.size(); c++) res.push_back(d*ind[c]);
        return res;
    }

    // Class Population

    void Population::add(Individual& ind){
        Population::iterator low;
        low = lower_bound(this->begin(),this->end(),ind);
        this->insert(low, ind);
    }

    Individual& Population::best_individual(){
        return this->begin()[0];
    }

    Individual& Population::worst_individual(){
        return this->end()[-1];
    }


    // Class DE

    DE::DE(int D, int NP, Ranges ranges, Objetive* objetive, RNG rng){
           this->D = D;
           this->NP = NP;
           this->rng = rng;
           set_ranges(ranges);
           set_objetive(objetive);
           set_population(build_population());
    }

    Individual DE::build_individual(){

        Individual ind;

        for(int i=0; i<D; i++){
            Range r = ranges[i];
            ind.push_back(r.first+rng.rand(1)*(r.second-r.first));
        }
        ind.obj_value = objetive->eval(ind);
        return ind;
    }

    Population DE::build_population(){

        Population pop;

        for(int i=0; i<NP ; i++){
            Individual ind = build_individual();
            pop.add(ind);
        }

        return pop;
    }

    void DE::evolve(){
        
        int target_index = get_target_index();
        Individual target = population[target_index];
        Individual donor = get_donor(target_index);
        Individual trial = crossover(target, donor);

        trial.obj_value = objetive->eval(trial);

        if (trial <= target) {
            population.erase(population.begin() + target_index);
            population.add(trial);
        }
                
    }

    int DE::get_target_index(){
        return rng.randint(get_NP());
    }

}
