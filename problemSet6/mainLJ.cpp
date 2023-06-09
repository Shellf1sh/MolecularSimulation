#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"



double LJ(double deltar, double epsilon, double sigma){
    //std::cout << "The distance is: " << deltar << std::endl;
        return 4*epsilon*(12*std::pow(12, sigma)/std::pow(13, deltar) - 6*std::pow(6, sigma)/std::pow(7, deltar));
    }

int main(){


    const int N = 64;
    const double m = 1;
    const double epsilon = 1.;
    const double sigma = 1.5;
    const double tau = sqrt(m*sigma*sigma/epsilon);
    const double dt = 0.002*tau;
    const int n_steps = 60*tau/dt;
    const double L = 6*sigma;

    std::cout << "Number of particles: " << N << std::endl;
    std::cout<< "The box side length is "<<L<<std::endl;
    std::cout<< "So the border is at +/- "<<L/2.0<< std::endl;
    std::cout<< "The simulation will take "<<n_steps<< " steps" << std::endl;
    VerletSimulation simulation(N, n_steps, m, dt, L, sigma, epsilon, true);
    
    simulation.initialize(4);
    

    //Initialise the array for storing the data
    double **data = new double *[n_steps];
    for (int i = 0; i < n_steps; i++){
        data[i] = new double[1 + 1 + 1]; //dt Ekin Epot   
    }

    double Ekin_array[n_steps];
    double Epot_array[n_steps];
    double E_array[n_steps];

    for(int i = 0; i < n_steps; i++){
        simulation.VerletStep(i);

        if(i == n_steps/2){
            double H_average = timeAverage(E_array, n_steps/2, 0);
            double actual_energy = simulation.getKineticEnergy() + simulation.getPotentialEnergy();

            double energy_diff = H_average - actual_energy;

            double scaling_const = (energy_diff + simulation.getKineticEnergy())/simulation.getKineticEnergy();

            std::cout << "The energy difference is " << energy_diff << " and the scaling constant is " << scaling_const << std::endl;


            simulation.rescaleMomenta(scaling_const);
        }

        data[i][0] = dt*i;
        double Ekin = simulation.getKineticEnergy();
        double Epot = simulation.getPotentialEnergy();

        data[i][1] = Ekin;
        data[i][2] = Epot;

        Ekin_array[i] = Ekin;
        Epot_array[i] = Epot;
        E_array[i] = Ekin + Epot;


        //std::cout << "Ekin=" << Ekin << " and Epot=" << Epot << std::endl;
    }

    create_csv(data, "Energy61cTherm.csv", n_steps, 3);


    double Ekin_time_average = timeAverage(Ekin_array, n_steps/2, n_steps/2);
    double Epot_time_average = timeAverage(Epot_array, n_steps/2, n_steps/2);
    

    std::cout << "The time averages are <Ekin>=" << Ekin_time_average << "and <Epot>=" << Epot_time_average << std::endl;

    return 0;
}