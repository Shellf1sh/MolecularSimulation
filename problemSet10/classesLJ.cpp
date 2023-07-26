#include<iostream>
#include<iomanip>
#include<cmath>
#include<random>
#include<ctime>

#include "classesLJ.h"
#include "functionsLJ.h"


std::seed_seq seed{std::time(nullptr)};

std::default_random_engine generator(seed);
std::normal_distribution<double> normal(0, 1);


VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt,  double epsilon,  double gamma, bool force)
        : N(N),n_steps(n_steps), m(m), dt(dt), epsilon(epsilon),  gamma(gamma), force(force)
    {
        beads.resize(N);
    }


void VerletSimulation::printMomenta(){
    std::cout << "Momneta = (" << beads[0].x << ", " << beads[0].y << ", " << beads[0].z << ")" << std::endl;

}


void VerletSimulation::zeroInit(){
    for (int i = 0; i < N; ++i) {
            beads[i].x = 0;
            beads[i].y = 0;
            beads[i].z = 0;
            beads[i].px = 0;
            beads[i].py = 0;
            beads[i].pz = 0;
            beads[i].fx = 0;
            beads[i].fy = 0;
            beads[i].fz = 0;
        }
    computeForces();

    if(force){
        std::cout << "The force is turned on" << std::endl;
    }else{
        std::cout << "The force is turned off" << std::endl;
    }
}

void VerletSimulation::updateMomenta(){
    for (int i = 0; i < N; i++) {
            beads[i].px += 0.5 * beads[i].fx * dt;
            beads[i].py += 0.5 * beads[i].fy * dt;
            beads[i].pz += 0.5 * beads[i].fz * dt;
        }
}

void VerletSimulation::updatePositions() {
        for (int i = 0; i < N; i++) {
            beads[i].x += 0.5 * beads[i].px/m * dt;
            beads[i].y += 0.5 * beads[i].py/m * dt;
            beads[i].z += 0.5 * beads[i].pz/m * dt;
        }
    }

void VerletSimulation::updateRandom(){
    for(int i = 0; i < N; i++){
        double Nx = normal(generator);
        double Ny = normal(generator);
        double Nz = normal(generator);

        beads[i].px = beads[i].px * exp(-gamma*dt) + sqrt(m*epsilon*(1-exp(-2*gamma*dt)))*Nx;
        beads[i].py = beads[i].py * exp(-gamma*dt) + sqrt(m*epsilon*(1-exp(-2*gamma*dt)))*Ny;
        beads[i].pz = beads[i].pz * exp(-gamma*dt) + sqrt(m*epsilon*(1-exp(-2*gamma*dt)))*Nz;

    }
}


void VerletSimulation::computeForces() {
    double distance;
    double force_magnitude;

    double prefactor;

    if(force){
        prefactor = 0.5 * m * omega*omega;
    }else{
        prefactor = 0;
    }

    for (int i = 0; i < N; ++i) {

        // Compute the forces
            double rx = (beads[i].x);
            double ry = (beads[i].y);
            double rz = (beads[i].z);

            /*distance = sqrt(rx*rx + ry*ry + rz*rz);
            
            force_magnitude = prefactor*distance*distance;*/

            double fx = 2*prefactor * rx; 
            double fy = 2*prefactor * ry; 
            double fz = 2*prefactor * rz; 

            //Maybe a sign needs to be switched
            beads[i].fx = -fx;
            beads[i].fy = -fy;
            beads[i].fz = -fz;

           // std::cout << "The force is " << fx << std::endl;

        }
    }


void VerletSimulation::VerletStep(int step) {
    updateMomenta();
    updatePositions();
    updateRandom();
    updatePositions();
    computeForces();
    updateMomenta();
}


double* VerletSimulation::getPositionsAndMomenta(){
    double* positionAndMomentum = new double[2*3*N];

    for(int i = 0; i < N; i++){
        positionAndMomentum[3*i] = beads[i].x;
        positionAndMomentum[3*i+1] = beads[i].y;
        positionAndMomentum[3*i+2] = beads[i].z;

        positionAndMomentum[3*i + 3*N] = beads[i].px;
        positionAndMomentum[3*i+1 + 3*N] = beads[i].py;
        positionAndMomentum[3*i+2 + 3*N] = beads[i].pz;
    }

    return positionAndMomentum;
}


double* VerletSimulation::getKineticEnergy(){
    double* Ekin = new double[N];

    for(int i = 0; i < N; i++){
        Ekin[i] = (beads[i].px*beads[i].px + beads[i].py*beads[i].py + beads[i].pz*beads[i].pz)/(2*m);
    }
    return Ekin;

}

double* VerletSimulation::getPotentialEnergy(){
    double* Epot = new double[N];

    for(int i = 0; i < N; i++){
        Epot[i] = 0.5* m*omega*omega*(beads[i].x*beads[i].x + beads[i].y*beads[i].y + beads[i].z*beads[i].z);
    }
    
    return Epot;
}



