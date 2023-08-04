#include<iostream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<random>
#include "classesLJ.h"
#include "functionsLJ.h"

std::seed_seq seed{std::time(nullptr)};

std::default_random_engine generator_normal(seed);
std::normal_distribution<double> normal(0, 1);


VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt, double L, double gamma, double kT, double rMax, double K, bool loop)
        : N(N),n_steps(n_steps), m(m), dt(dt), L(L), gamma(gamma), kT(kT), rMax(rMax), K(K), loop(loop)
    {
        beads.resize(N);
    }

void VerletSimulation::initialize() {
    for(int i = 0; i < N; i++){
        beads[i].px = MaxwellBoltzmann(kT, m);
        beads[i].py = MaxwellBoltzmann(kT, m);
        beads[i].pz = MaxwellBoltzmann(kT, m);
    }

    for(int i = 0; i < N; i++){
        beads[i].x = i*rMax/2 - L/2;
        beads[i].y = 0;
        beads[i].z = 0;
    }
}


void VerletSimulation::testInit(){
    beads[0].x = -0.2;
    beads[0].y = 0;
    beads[0].z = 0;
    beads[1].x = 0.2;
    beads[1].y = 0;
    beads[1].z = 0;

    beads[0].px = 0;
    beads[0].py = 0.5;
    beads[0].pz = 0;
    beads[1].px = 0;
    beads[1].py = 0;
    beads[1].pz = 0;
}

void VerletSimulation::updateMomenta(){
    for (int i = 0; i < N; ++i) {
            beads[i].px += 0.5 * beads[i].fx * dt;
            beads[i].py += 0.5 * beads[i].fy * dt;
            beads[i].pz += 0.5 * beads[i].fz * dt;
        }
    }

void VerletSimulation::updatePositions(bool hamil) {
        double factor = 0.5;
        if(hamil){
            factor = 1.0;
        }
        
        for (int i = 0; i < N; ++i) {
            beads[i].x += factor*beads[i].px/m * dt;
            beads[i].y += factor*beads[i].py/m * dt;
            beads[i].z += factor*beads[i].pz/m * dt;
        }
    }

void VerletSimulation::backintotheBox() {
    for (int i = 0; i < N; ++i) {
        beads[i].x = boundary(beads[i].x, L);
        beads[i].y = boundary(beads[i].y, L);
        beads[i].z = boundary(beads[i].z, L);
    }
}

double VerletSimulation::FENEForce(double deltar){
    return -K*deltar/(1-deltar*deltar/(rMax*rMax));
}
double VerletSimulation::FENEPotential(double deltar){
    return -0.5*K*rMax*rMax*log(1-(deltar/rMax)*(deltar/rMax));
}

void VerletSimulation::computeForces() {
    double fx;
    double fy;
    double fz;

    for (int i = 0; i < N; ++i) {
        beads[i].fx = 0.0;
        beads[i].fy = 0.0;
        beads[i].fz = 0.0;
    }

    for (int i = 0; i < N-1; ++i) {
        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        // Compute the forces
        double rx = (beads[i].x - beads[i+1].x);
        double ry = (beads[i].y - beads[i+1].y);
        double rz = (beads[i].z - beads[i+1].z);

        rx = boundary(rx, L);
        ry = boundary(ry, L);
        rz = boundary(rz, L);

        double distance = sqrt(rx*rx+ry*ry+rz*rz);
        double total_force = FENEForce(distance);

        fx = abs(total_force)*rx/distance; 
        fy = abs(total_force)*ry/distance; 
        fz = abs(total_force)*rz/distance; 

        //Maybe a sign needs to be switched
        beads[i].fx += -fx;
        beads[i].fy += -fy;
        beads[i].fz += -fz;

        beads[i+1].fx += fx;
        beads[i+1].fy += fy;
        beads[i+1].fz += fz;  
    }
    if(loop){
        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        // Compute the forces
        double rx = (beads[N-1].x - beads[0].x);
        double ry = (beads[N-1].y - beads[0].y);
        double rz = (beads[N-1].z - beads[0].z);

        rx = boundary(rx, L);
        ry = boundary(ry, L);
        rz = boundary(rz, L);

        double distance = sqrt(rx*rx+ry*ry+rz*rz);
        double total_force = FENEForce(distance);

        fx = abs(total_force)*rx/distance; 
        fy = abs(total_force)*ry/distance; 
        fz = abs(total_force)*rz/distance; 

        //Maybe a sign needs to be switched
        beads[N-1].fx += fx;
        beads[N-1].fy += fy;
        beads[N-1].fz += fz;

        beads[0].fx += -fx;
        beads[0].fy += -fy;
        beads[0].fz += -fz;  
    }
}

void VerletSimulation::updateRandom(){
    for(int i = 0; i < N; i++){
        double Nx = normal(generator_normal);
        double Ny = normal(generator_normal);
        double Nz = normal(generator_normal);

        beads[i].px = beads[i].px * exp(-gamma*dt) + sqrt(m*kT*(1-exp(-2*gamma*dt)))*Nx;
        beads[i].py = beads[i].py * exp(-gamma*dt) + sqrt(m*kT*(1-exp(-2*gamma*dt)))*Ny;
        beads[i].pz = beads[i].pz * exp(-gamma*dt) + sqrt(m*kT*(1-exp(-2*gamma*dt)))*Nz;


    }
}


void VerletSimulation::VerletStep(bool hamil) {
            updateMomenta();
            updatePositions(hamil);
            backintotheBox();
            if(!hamil){
                updateRandom();
                updatePositions(hamil);
                backintotheBox();
            }
            computeForces();
            updateMomenta();
        }


double* VerletSimulation::getPositionsAndMomenta(){
    double* positionAndMomentum = new double[2*3*N];

    for(int i = 0; i < N; i++){
        positionAndMomentum[3*i] = beads[i].x;
        positionAndMomentum[3*i+1] = beads[i].y;
        positionAndMomentum[3*i+2] = beads[i].z;

        positionAndMomentum[3*N + 3*i] = beads[i].px;
        positionAndMomentum[3*N + 3*i+1] = beads[i].py;
        positionAndMomentum[3*N + 3*i+2] = beads[i].pz;
    }

    return positionAndMomentum;
}


double VerletSimulation::getKineticEnergy(){
    double Ekin = 0.0;

    for(int i = 0; i < N; i++){
        Ekin += (beads[i].px*beads[i].px + beads[i].py*beads[i].py + beads[i].pz*beads[i].pz)/(2*m);
    }
    return Ekin;
}

double VerletSimulation::getPotentialEnergy(){
    double Epot = 0.0;
    double distance = 0.0;
    double dx = 0;
    double dy = 0;
    double dz = 0;


    for(int i = 0; i < N-1; i++){
        dx = beads[i].x-beads[i+1].x;
        dy = beads[i].y-beads[i+1].y;
        dz = beads[i].z-beads[i+1].z;

        dx = boundary(dx, L);
        dy = boundary(dy, L);
        dz = boundary(dz, L);

        double distance = sqrt(dx*dx + dy*dy + dz*dz);

        Epot += FENEPotential(distance);
    }

    if(loop){ //The potential energy between the 2 ends
        dx = beads[N-1].x-beads[0].x;
        dy = beads[N-1].y-beads[0].y;
        dz = beads[N-1].z-beads[0].z;

        dx = boundary(dx, L);
        dy = boundary(dy, L);
        dz = boundary(dz, L);

        double distance = sqrt(dx*dx + dy*dy + dz*dz);

        Epot += FENEPotential(distance);
    }

    return Epot;
}

