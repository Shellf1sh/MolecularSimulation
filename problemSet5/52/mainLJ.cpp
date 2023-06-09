#include<iostream>
#include<cmath>
#include "functionsLJ.h"
#include "classesLJ.h"


int main(){

    std::cout << "Test1" << std::endl;

    const int k = 3;

    const int perDim = 2;
    const int N = std::pow(perDim, 3);
    //const int N = std::pow(2, 3*k);

    std::cout << "Number of particles: " << N << std::endl;

    const double m = 1;
    
    const double epsilon = 1.;
    const double sigma = 1.5;
    const double tau = sqrt(m*sigma*sigma/epsilon);
    const double dt = 0.002*tau;
    //const int n_steps = 100;
    const int n_steps = 10*tau/dt;
    const double L = std::pow(2, k)*sigma;


    VerletSimulation simulation(N, n_steps, m, dt, L, sigma, epsilon, false);
    simulation.initialize(std::pow(2, k));

    //double data[n_steps][N*3+N*3+1];  //3 positions, 3 momenta, 1 time and 1 energy variable for every particle

    //std::vector<std::vector<double>> data(n_steps, std::vector<double>(N * 3 + N * 3 + 1));
    double **data = new double *[n_steps];
    for (int i = 0; i < n_steps; i++){
        data[i] = new double[N * 3 + N * 3 + 1];    
    }

    double* InitialPosAndMom = simulation.getPositionsAndMomenta();

    for(int j = 0; j < 2*3*N; j++){
        data[0][j] = *(InitialPosAndMom + j);
    }
    data[0][2*3*N] = *(InitialPosAndMom + 2*3*N);

    for(int i = 0; i < N; i++){
        std::cout << "(" << data[0][3*i] << ", " << data[0][3*i+1] << ", " << data[0][3*i+2] << ")" << std::endl;
    }

    for(int i = 0; i < 3*N; i++){
        std::cout << "(" << data[0][i] << ") " << std::endl;
    }

    for (int i = 0; i < n_steps; i++) {
        simulation.VerletStep();  // Perform a single Verlet step

        double* posAndMom = simulation.getPositionsAndMomenta();

        for(int j = 0; j < 2*3*N; j++){
           data[i][j] = *(posAndMom + j);
        }
        delete[] posAndMom;

        data[i][2*3*N] = i * dt;
    }

    
    create_csv(data, "52aNormalDeltat.csv", n_steps, N*3+N*3+1);
    //create_csv(data, "kequal1.csv", n_steps, N * 3 + N * 3 + 1);


         
    VerletSimulation simulationdt(N, n_steps, m, 2*dt, L, sigma, epsilon, false);
    simulationdt.initialize(perDim);

    //double datadt[n_steps][N*3+N*3+1];  //3 positions, 3 momenta, 1 time and 1 energy variable for every particle
    //std::vector<std::vector<double>> datadt(n_steps, std::vector<double>(N * 3 + N * 3 + 1));

    double **datadt = new double *[n_steps];
    for (int i = 0; i < n_steps; i++){
        datadt[i] = new double[N * 3 + N * 3 + 1];    
    }

    std::cout << "Test2" << std::endl;

    double* InitialPosAndMomdt = simulationdt.getPositionsAndMomenta();

    for(int j = 0; j < 2*3*N; j++){
        datadt[0][j] = *(InitialPosAndMomdt + j);
    }
    datadt[0][2*3*N] = *(InitialPosAndMomdt + 2*3*N);

    for (int i = 0; i < n_steps; i++) {
        simulationdt.VerletStep();  // Perform a single Verlet step

        double* posAndMom = simulationdt.getPositionsAndMomenta();

        for(int j = 0; j < 2*3*N; j++){
        datadt[i][j] = *(posAndMom + j);
        }
        delete[] posAndMom;

        datadt[i][2*3*N] = i * dt;
    }

    std::cout << "Test3" << std::endl;

    create_csv(datadt, "52aDoubleDeltat.csv", n_steps, N * 3 + N * 3 + 1);
    


    std::cout << "Test4" << std::endl;

    return 0;
}