#include <vector>
#include <cstdlib>
#include <random>
#include <time.h>

#ifndef DIFFEVOL
#define DIFFEVOL

namespace diffevol 
{

    typedef std::pair<double,double> Range;
    typedef std::vector<Range> Ranges;

    class RNG
    {
        private:
            
            std::default_random_engine rng;

        public:
            RNG();

            virtual int randint(int N);
    
            virtual double rand(int R);
    };


    class Individual: public std::vector<double>
    {
        public:
            Individual(): std::vector<double>() {};
            Individual(vector<double> vect): std::vector<double>(vect) {};

            double obj_value;
        
            friend bool operator<(const Individual &ind1, const Individual &ind2);
            friend bool operator<=(const Individual &ind1, const Individual &ind2);
            friend Individual operator+(const Individual &ind1, const Individual &ind2);
            friend Individual operator-(const Individual &ind1, const Individual &ind2);
            friend Individual operator*(const Individual &ind, const double d);
            friend Individual operator*(const double d, const Individual &ind);
    };


    class Population: public std::vector<Individual>
    {
        public:
            void add(Individual& ind);

            Individual& best_individual();
            
            Individual& worst_individual();
    };

    class Objetive
    {
        public:
            Objetive(){};
            virtual double eval(Individual& ind){return 0;};
    };

    class DE
    {
        protected:

            int D; /*! dimensions of the parameter space */
            int NP; /*! number of individuals in population*/
            Ranges ranges; /*! ranges for each parameter*/
            Population population; /*! population */
            Objetive* objetive; /*! Objetive */
        
        public:

            RNG rng;
    
            //! Constructor
            DE(){}
            DE(int D, int NP, Ranges ranges, Objetive* objetive, RNG rng=RNG());

            virtual Individual build_individual();

            virtual Population build_population();

            virtual void evolve();

            virtual int get_target_index();

            //! virtual functions that must be implemented in the derived classes
            virtual Individual get_donor(int target_index)=0;
            virtual Individual crossover(Individual& target, Individual& donor)=0;

            //! setters and getters
            void set_objetive( Objetive* obj){objetive = obj;}

            void set_population(Population pop){population = pop;}
            Population& get_population(){return population;}

            void set_ranges(Ranges ranges){this->ranges=ranges;}
            Ranges get_ranges(){return ranges;}

            int get_D(){return D;}
            int get_NP(){return NP;}
    };

}

#endif
