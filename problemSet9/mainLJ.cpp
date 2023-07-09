#include<iostream>
#include<cmath>
#include<random>
#include<ctime>

#include "functionsLJ.h"
#include "classesLJ.h"



int main(){
    const int N = 64;
    const double m = 1;
    const double epsilon = 1.;
    const double sigma = 1.5;
    const double tau = sqrt(m*sigma*sigma/epsilon);
    const double dt = 0.002*tau;
    const int n_steps = 10*tau/dt;
    const double L = 6*sigma;

    const int N_test = 50000;

    std::seed_seq seed{std::time(nullptr)};

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> uniform(-L/2, L/2);

    std::cout << "Number of particles: " << N << std::endl;
    std::cout<< "The box side length is "<< L << ", so the border is at +/- "<<L/2.0<< std::endl;
    std::cout<< "The simulation will take "<< n_steps << " steps" << std::endl;
    VerletSimulation simulation(N, n_steps, m, dt, L, sigma, epsilon, true);
    
    simulation.initialize(4);

    //Initialise the array for storing the data
    double **data = new double *[n_steps/50];
    for (int i = 0; i < n_steps/50; i++){
        data[i] = new double[N_test]; //dt Ekin Epot   
    }

    std::cout<< "There will be " << n_steps/50 << " \"measurements\" of muEx"<<std::endl;

    for(int i = 0; i < n_steps; i++){

        simulation.VerletStep(i);

        if((i+1)%50 == 0){ //Singles out the step before the momenta reset

            double test_positions[3*N_test];
            for(int m = 0; m < 3*N_test; m++){
                test_positions[m] = uniform(generator);
            }

            double *real_positions = simulation.getPositionsAndMomenta();

            for(int j = 0; j < 3*N_test; j += 3){ //For every test particle we calculate the distance and potential energy to every real particle
                double deltaU = 0;

                //std::cout << "j= " << j << std::endl;
                //std::cout << "i is " << i << " and the index is " << (i+1)/50 -1 << std::endl;
                
                for(int k = 0; k < 3*N; k += 3){
                    double distance = calculateDistance(test_positions[j], test_positions[j+1], test_positions[j+2], *(real_positions + k), *(real_positions + k+1), *(real_positions + k+2));

                    double pot = LJ_potential(distance, epsilon, sigma);
                    deltaU += pot;
                }
                if(std::isnan(deltaU)){
                    //std::cout << "Delta U is NaN at step " << (i+1)/50 - 1 << std::endl;
                    deltaU = 0;
                }

                //std::cout << "Delta U: " << deltaU << std::endl;
                //std::cout << "I insert the values at (" << (i+1)/50 - 1 <<", " <<j/3<<")"<<std::endl;
                //std::cout << std::endl;

                data[(i+1)/50 - 1][j/3] = deltaU;


            }
        }

    }

    create_csv(data, "muEx50000.csv", n_steps/50, N_test);


    return 0;
}