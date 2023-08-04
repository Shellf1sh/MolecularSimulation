#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"


int main(){

    const int perDim = 4;
    
    const int particleNumber = 2* std::pow(perDim, 3);
    const double density = 0.55; //Personal parameter

    const double m = 1;
    
    const double epsilon = 1.;
    const double kT = 1.7*epsilon; //kT/epsilon = 1.7, Personal parameter
    const double sigma = 1; 
    const double tau = 1; //sqrt(m*sigma*sigma/epsilon);
    const double gamma = 1/tau;
    const double dt = 0.001*tau;
    const int n_steps = 10*tau/dt;
    const double L = std::cbrt(particleNumber/density);
    const double particleSpace = L/perDim;

    bool HamiltonianDyn = true;
    const int ThermostatCutOff = 2*n_steps;

    const int saveSteps = 10;

    std::cout << "Number of particles: " << particleNumber << std::endl;
    std::cout << "The simulation will take "<< n_steps << " steps, but we only save " << n_steps/saveSteps << std::endl;
    std::cout << "We save every " << saveSteps << "th step" << std::endl;
    std::cout << "The box has side length " << L << " and the side length of a unit cell is " << particleSpace << std::endl;
    if(HamiltonianDyn){
        std::cout << "We turn off the thermostat off after the " << ThermostatCutOff << "th step" << std::endl;
    }

    VerletSimulation simulation(particleNumber, n_steps, m, dt, L, sigma, epsilon, gamma, kT); 

    simulation.initialize(particleNumber, particleSpace); 
    //simulation.testInit();

    double **data = new double *[n_steps/saveSteps];
    for (int i = 0; i < n_steps/saveSteps; i++){
        data[i] = new double[particleNumber*3 + particleNumber*3 + 1 + 2 ];//Position, momenta, forces, Ekin, Epot, time   
    }

    for (int i = 0; i < n_steps; i++) {
            
        HamiltonianDyn = i > ThermostatCutOff; // If true then no random influence

        simulation.VerletStep(HamiltonianDyn);  // Perform a single Verlet step

        
        if(i%saveSteps == 0){

            double* posAndMom = simulation.getPositionsAndMomenta();

            for(int j = 0; j < 2*3*particleNumber; j++){
                data[i/saveSteps][j] = *(posAndMom + j);
            }
            delete[] posAndMom;

            /*double* forces = simulation.getForces();

            for(int j = 0; j < particleNumber; j++){
                data[i/saveSteps][j] = *(forces + j + 2*3*particleNumber);
            }*/

            data[i/saveSteps][2*3*particleNumber] = simulation.getKineticEnergy();
            data[i/saveSteps][2*3*particleNumber+1] = simulation.getPotentialEnergy();

            data[i/saveSteps][2*3*particleNumber+2] = i * dt;
        }
    }

    
    create_csv(data, "position_test.csv", n_steps/saveSteps, 2*particleNumber*3+3);


    return 0;
}