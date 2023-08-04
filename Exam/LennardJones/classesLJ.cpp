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


VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt, double L, double sigma, double epsilon, double gamma, double kT)
        : N(N),n_steps(n_steps), m(m), dt(dt), L(L), sigma(sigma), epsilon(epsilon), gamma(gamma), kT(kT)
    {
        beads.resize(N);
    }

void VerletSimulation::initialize(int particleNumber, double spacing) {
        for(int i = 0; i < N; i++){
            beads[i].px = MaxwellBoltzmann(kT, m);
            beads[i].py = MaxwellBoltzmann(kT, m);
            beads[i].pz = MaxwellBoltzmann(kT, m);
        }

        //int perDim = ceil(std::cbrt(particleNumber));

        const int dim = std::cbrt(N/2); //The size of a unit cell

        std::cout << "Unit cells per dimension is " << dim << std::endl;

        int index = 0;

        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                for(int k = 0; k < dim; k++){
                    beads[index].x = i * spacing - L/2;
                    beads[index].y = j * spacing - L/2;
                    beads[index].z = k * spacing - L/2;
                    index++;

                    //std::cout << "("<< i * particleSpacing - L/2<< ", "<<j * particleSpacing - L/2 << ", " << k * particleSpacing - L/2 << ")"<< std::endl;

                    beads[index].x = (i + 0.5) * spacing - L/2;
                    beads[index].y = (j + 0.5) * spacing - L/2;
                    beads[index].z = (k + 0.5) * spacing - L/2;
                    index++;

                    //std::cout << "("<< (i+0.5) * particleSpacing - L/2<< ", "<<(j+0.5) * particleSpacing - L/2 << ", " << (k+0.5) * particleSpacing - L/2 << ")"<< std::endl;

                }
            }
        }
    }


void VerletSimulation::testInit(){
    beads[0].x = -2*L/3;
    beads[0].y = 0;
    beads[0].z = 0;
    beads[1].x = 0;
    beads[1].y = 0;
    beads[1].z = 0;
    beads[2].x = 2*L/3;
    beads[2].y = 0;
    beads[2].z = 0;

    beads[0].px = 1;
    beads[0].py = 0;
    beads[0].pz = 0;
    beads[1].px = 0;
    beads[1].py = 0;
    beads[1].pz = 0;
    beads[2].px = 0;
    beads[2].py = 0;
    beads[2].pz = 0;
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

void VerletSimulation::computeForces() {
    for (int i = 0; i < N; ++i) {
        beads[i].fx = 0.0;
        beads[i].fy = 0.0;
        beads[i].fz = 0.0;
    }

    for (int i = 0; i < N; ++i) {
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;

        // Compute the forces
        for(int j = i+1; j < N; j++) {
            double rx = (beads[i].x - beads[j].x);
            double ry = (beads[i].y - beads[j].y);
            double rz = (beads[i].z - beads[j].z);

            rx = boundary(rx, L);
            ry = boundary(ry, L);
            rz = boundary(rz, L);


            double distance = sqrt(rx*rx+ry*ry+rz*rz);
            double total_force = LennardJones(distance);

            fx = abs(total_force)*rx/distance; 
            fy = abs(total_force)*ry/distance; 
            fz = abs(total_force)*rz/distance; 

            beads[i].fx += fx;
            beads[i].fy += fy;
            beads[i].fz += fz;

            beads[j].fx += -fx;
            beads[j].fy += -fy;
            beads[j].fz += -fz;  
        } 
    }
}

double VerletSimulation::LennardJones(double deltar){
    return 6*4*epsilon*(2*std::pow(sigma, 12)/std::pow(deltar, 13) - std::pow(sigma, 6)/std::pow(deltar, 7));
}

void VerletSimulation::backintotheBox() {
    for (int i = 0; i < N; ++i) {
        beads[i].x = boundary(beads[i].x, L);
        beads[i].y = boundary(beads[i].y, L);
        beads[i].z = boundary(beads[i].z, L);

        
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


double* VerletSimulation::getForces(){

    double* forces = new double[3*N];

    for(int i; i < N; i++){
        forces[3*i] = beads[i].fx;
        forces[3*i+1] = beads[i].fy;
        forces[3*i+2] = beads[i].fz;
    }

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

    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double rx = (beads[i].x - beads[j].x);
            double ry = (beads[i].y - beads[j].y);
            double rz = (beads[i].z - beads[j].z);

            rx = boundary(rx, L);
            ry = boundary(ry, L);
            rz = boundary(rz, L);
            

            double distance = sqrt(rx*rx+ry*ry+rz*rz);

            Epot += 4*epsilon*(std::pow(sigma/distance, 12) - std::pow(sigma/distance, 6));
        }
    }
    return Epot;
}

