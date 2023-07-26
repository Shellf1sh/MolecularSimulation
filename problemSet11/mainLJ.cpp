#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"


int main(){

    const int perDim = 3;
    const int N = std::pow(perDim, 3);



    const double m = 1;
    
    const double epsilon = 1.;
    const double kT = 1.5*epsilon;
    const double sigma = 1;
    const double tau = sqrt(m*sigma*sigma/epsilon);
    const double gamma = 1/tau;
    const double dt = 0.002*tau;
    const int n_steps = 200*tau/dt;
    const double L = 6*sigma;

    bool HamiltonianDyn = true;
    const int ThermostatCutOff = n_steps/6;

    const int saveSteps = 20;

    std::cout << "Number of particles: " << N << std::endl;
    std::cout << "The simulation will take "<< n_steps << " steps, but we only save " << n_steps/saveSteps << std::endl;
    std::cout << "Gamma is " << gamma << std::endl; 
    std::cout << "The exponetial factor is " << gamma*dt << std::endl;
    std::cout << "We save every " << saveSteps << "th step" << std::endl;
    if(HamiltonianDyn){
        std::cout << "We turn off the thermostat off after the " << ThermostatCutOff << "th step" << std::endl;
    }

    VerletSimulation simulation(N, n_steps, m, dt, L, sigma, epsilon, gamma, kT);
    simulation.initialize(perDim);


    double **data = new double *[n_steps/saveSteps];
    for (int i = 0; i < n_steps/saveSteps; i++){
        data[i] = new double[N * 3 + N * 3 + 1 + 2 ];//Position, momenta, Ekin, Epot, time   
    }

    for (int i = 0; i < n_steps; i++) {
            
        HamiltonianDyn = i < ThermostatCutOff;

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

    
    create_csv(data, "Maybe27.csv", n_steps/saveSteps, N*3+N*3+3);


    return 0;
}