#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"


int main(){

    const int perDim = 4;
    
    const int particleNumber = 2* std::pow(perDim, 3);
    const double density = 0.55; //Personal parameter

    const double m = 1;

    const int test_particles = 10000;
    const bool Widoms = false;

    const double epsilon = 1.7;//Person paramter
    const double kT = epsilon; //kT/epsilon = 1.7, Personal parameter
    const double sigma = 1; 
    const double tau = 1; //sqrt(m*sigma*sigma/epsilon);
    const double gamma = 1/tau;
    const double dt = 0.001*tau * std::pow(2, 4);
    const int n_steps = 20*tau/dt;
    const double L = std::cbrt(particleNumber/density);
    const double particleSpace = L/perDim;

    bool HamiltonianDyn = true;
    const int ThermostatCutOff = 2*n_steps;

    const int saveSteps = 1;

    std::cout << "Number of particles: " << particleNumber << std::endl;
    std::cout << "The simulation will take "<< n_steps << " steps, but we only save " << n_steps/saveSteps << std::endl;
    std::cout << "We save every " << saveSteps << "th step" << std::endl;
    std::cout << "The box has side length " << L << " and the side length of a unit cell is " << particleSpace << std::endl;
    if(HamiltonianDyn){
        std::cout << "We turn off the thermostat off after the " << ThermostatCutOff << "th step" << std::endl;
    }
    if(Widoms){
        std::cout << "We insert " << test_particles << " test particles to measure the chemical potential" << std::endl;
    }

    VerletSimulation simulation(particleNumber, n_steps, m, dt, L, sigma, epsilon, gamma, kT); 

    simulation.initialize(particleNumber, particleSpace); 
    //simulation.testInit();

    double **data = new double *[n_steps/saveSteps];
    for (int i = 0; i < n_steps/saveSteps; i++){
        data[i] = new double[particleNumber*3 + particleNumber*3 + 5 ];//Position, momenta, chemical potential, virial, Ekin, Epot, time   
    }


    for (int i = 0; i < n_steps; i++) {
            
        HamiltonianDyn = i < ThermostatCutOff; // If true then no random influence
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
            data[i/saveSteps][2*3*particleNumber+1] = simulation.getVirial();
            data[i/saveSteps][2*3*particleNumber+2] = simulation.getKineticEnergy();
            data[i/saveSteps][2*3*particleNumber+3] = simulation.getPotentialEnergy();
            data[i/saveSteps][2*3*particleNumber+4] = i * dt;

            if(Widoms){
                double avg_exp = simulation.getChemicalPotential(test_particles);
                data[i/saveSteps][2*3*particleNumber] = avg_exp;
                std::cout << i << ": " << avg_exp << std::endl;
            }else{
                data[i/saveSteps][2*3*particleNumber] = 0;
            }
        }
    }

    
    create_csv(data, "k4.csv", n_steps/saveSteps, 2*3*particleNumber+5);


    return 0;
}