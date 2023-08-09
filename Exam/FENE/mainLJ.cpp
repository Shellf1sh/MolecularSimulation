#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"


int main(){
    const int N = 48;
    const bool periodic = false;

    const double m = 1;
    const double K = 1;
    const double rMax = 1;
    
    const double kT = 0.1*K*rMax*rMax;
    const double sigma = 1;
    const double gamma = 0.2*sqrt(K/m);
    const double tau = 1/gamma;
    const double dt = 0.001*tau;
    const int n_steps = 200*tau/dt;

    double L = N*rMax;
    if(periodic){
        L = N*rMax/3;
    }

    bool HamiltonianDyn = false; // If true then no random influence
    const int ThermostatCutOff = 0;

    const int saveSteps = 10;

    std::cout << "Number of particles: " << N << std::endl;
    std::cout << "The simulation will take "<< n_steps << " steps, but we only save " << n_steps/saveSteps << std::endl;
    std::cout << "Gamma is " << gamma << std::endl; 
    std::cout << "The exponetial factor is " << gamma*dt << std::endl;
    std::cout << "We save every " << saveSteps << "th step" << std::endl;
    std::cout << "The box has side length " << L << std::endl;
    if(HamiltonianDyn){
        std::cout << "We turn off the thermostat off after the " << ThermostatCutOff << "th step" << std::endl;
    }

    VerletSimulation simulation(N, n_steps, m, dt, L, gamma, kT, rMax, K, periodic);
    simulation.initialize(periodic);
    //simulation.testInit();


    double **data = new double *[n_steps/saveSteps];
    for (int i = 0; i < n_steps/saveSteps; i++){
        data[i] = new double[N * 3 + N * 3 + 3 ];//Position, momenta, Ekin, Epot, time   
    }

    for (int i = 0; i < n_steps; i++) {

        simulation.VerletStep(HamiltonianDyn);  // Perform a single Verlet step
        
        if(i%saveSteps == 0){

            double* posAndMom = simulation.getPositionsAndMomenta();

            for(int j = 0; j < 2*3*N; j++){
                data[i/saveSteps][j] = *(posAndMom + j);
            }
            delete[] posAndMom;

            data[i/saveSteps][2*3*N] = simulation.getKineticEnergy();
            data[i/saveSteps][2*3*N+1] = simulation.getPotentialEnergy();

            data[i/saveSteps][2*3*N+2] = i * dt;

        }
    }

    
    create_csv(data, "test.csv", n_steps/saveSteps, N*3+N*3+3);


    return 0;
}