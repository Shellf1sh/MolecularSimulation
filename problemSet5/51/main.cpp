#include<iostream>
#include<cmath>
#include "functions.h"
#include "classes.h"


int main(){

    const int N = 11;
    //const int n_steps = 100;
    const double m = 1;
    const double k = 1.0;
    const double tau = sqrt(m/k);
    const double l = 1.0;
    const double L = 11*l;
    const double dt = 0.01;
    const double epsilon = 1.;
    double p0 = sqrt(m*epsilon);
    double zero = 0.;
    const int n_steps = 100*tau/dt;

    std::cout << "Initialized" << std::endl;
    
    std::vector<double> momenta = {-3./7.*p0, 6./7.* p0, -2./7.*p0};
    //std::vector<double> momenta = {-1./5.*p0, 2./5.* p0, -1./5.*p0};
    const int pyIndex = 6;

    std::cout << momenta[0] << std::endl;
    std::cout << momenta[1] << std::endl;
    std::cout << momenta[2] << std::endl;

    VerletSimulation simulation(N, n_steps, m, dt, k, l, L, false);
    simulation.initialize(zero, l, L, momenta, pyIndex);


    //double data[n_steps][N*3+N*3+1];  //3 positions, 3 momenta, 1 time and 1 energy variable for every particle
    double **data = new double *[n_steps];
    for (int i = 0; i < n_steps; i++){
        data[i] = new double[N * 3 + N * 3 + 1];    
    }

    double* InitialPosAndMom = simulation.getPositionsAndMomenta();

    for(int j = 0; j < 2*3*N; j++){
        data[0][j] = *(InitialPosAndMom + j);
    }



    for (int i = 1; i < n_steps; i++) {
        simulation.VerletStep();  // Perform a single Verlet step

        double* posAndMom = simulation.getPositionsAndMomenta();

        for(int j = 0; j < 2*3*N; j++){
           data[i][j] = *(posAndMom + j);
        }
        delete[] posAndMom;

        data[i][2*3*N] = i * dt;
    }

    create_csv(data, "51a.csv", n_steps, N*3+N*3+1);


    std::cout << "Test3" << std::endl;


    //Solving 5.1b where the chain is a loop
    VerletSimulation simulationLoop(N, n_steps, m, dt, k, l, L, true);
    simulationLoop.initialize(zero, l, L, momenta, pyIndex);

    //double dataLoop[n_steps][N*3+N*3+1];  //3 positions, 3 momenta, 1 time and 1 energy variable for every particle
    double **dataLoop = new double *[n_steps];
    for (int i = 0; i < n_steps; i++){
        dataLoop[i] = new double[N * 3 + N * 3 + 1];    
    }

    double* InitialPosAndMomLoop = simulationLoop.getPositionsAndMomenta();

    for(int j = 0; j < 2*3*N; j++){
        dataLoop[0][j] = *(InitialPosAndMomLoop + j);
    }


    for (int i = 1; i < n_steps; i++) {
        simulationLoop.VerletStep();  // Perform a single Verlet step

        double* posAndMomLoop = simulationLoop.getPositionsAndMomenta();

        //std::cout << *(posAndMom + 3*N) << std::endl;

        for(int j = 0; j < 2*3*N; j++){
           dataLoop[i][j] = *(posAndMomLoop + j);
        }
        delete[] posAndMomLoop;

        dataLoop[i][2*3*N] = i * dt;
    }

    create_csv(dataLoop, "51b.csv", n_steps, N*3+N*3+1);


    return 0;
}