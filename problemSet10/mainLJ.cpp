#include<iostream>
#include<cmath>
#include<random>
#include<ctime>

#include "functionsLJ.h"
#include "classesLJ.h"



int main(){
    const int N = 100;
    const bool force = false;

    const double m = 1;
    const double epsilon = 1;
    const double tau = 1;
    const double dt = 0.01*tau;

    const int n_steps = 10000*tau/dt;

    const double gamma = 1/tau;
    const double l = 1/gamma * sqrt(epsilon/m);

    const double omega = gamma/3;

    const int t_Observe = 0.1*gamma/dt;

    const double k = m*omega*omega;
    //const int N_test = 1;


    std::cout << "Number of particles: " << N << std::endl;
    std::cout<< "The simulation will take "<< n_steps << " steps and it will save " << n_steps/t_Observe << " data points" << std::endl;
    std::cout << "Gamma is " << gamma << " and omega is " << omega << std::endl; 
    std::cout << "The exponetial factor is " << gamma*dt << std::endl;
    std::cout << "l is " << l << std::endl;
    std::cout << "t_Observe is " << t_Observe << std::endl;


    VerletSimulation simulation(N, n_steps, m, dt, epsilon, gamma, force); //Last bool is for the force
    
    simulation.zeroInit();

    int steps_to_save = n_steps/t_Observe;

    //Initialise the array for storing the data
    double **data = new double *[steps_to_save];
    for (int i = 0; i < steps_to_save; i++){
        data[i] = new double[6*N + 2*N];  
    }

    for(int i = 0; i < n_steps; i++){

        simulation.VerletStep(i);

        if(i%t_Observe == 0){
            double* PositionAndMomenta = simulation.getPositionsAndMomenta();
            double* Ekin = simulation.getKineticEnergy();
            double* Epot = simulation.getPotentialEnergy();

            for(int j = 0; j < 6*N; j++){
                data[i/t_Observe][j] = *(PositionAndMomenta + j);
            }

            for(int b = 0; b < N; b++){
                data[i/t_Observe][6*N + b] = *(Ekin + b);
                data[i/t_Observe][6*N+N+b] = *(Epot + b);
            }
        }
       
    }

    create_csv(data, "freeparticlelong.csv", steps_to_save, 6*N + 2*N);


    return 0;
}