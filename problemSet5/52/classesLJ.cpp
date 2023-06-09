#include<iostream>
#include<iomanip>
#include<cmath>
#include "classesLJ.h"
#include "functionsLJ.h"



VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt, double L, double sigma, double epsilon, bool loop)
        : N(N),n_steps(n_steps), m(m), dt(dt), L(L), sigma(sigma), epsilon(epsilon), loop(loop)
    {
        beads.resize(N);
    }

void VerletSimulation::initialize(int perDim) {
        for(int i = 0; i < N; i++){
            beads[i].px = MaxwellBoltzmann(epsilon, m);
            beads[i].py = MaxwellBoltzmann(epsilon, m);
            beads[i].pz = MaxwellBoltzmann(epsilon, m);
        }

        int index = 0;
        for(int i = 0; i < perDim; i++){
            for(int j = 0; j < perDim; j++){
                for(int k = 0; k < perDim; k++){
                    //std::cout << "(" << i << ", " << j << ", " << k << ")" << std::endl;
                    beads[index].x = i * 1.5*sigma;
                    beads[index].y = j * 1.5*sigma;
                    beads[index].z = k * 1.5*sigma;
                    index++;
                }
            }
        }


    }

void VerletSimulation::updateMomenta(){
    for (int i = 0; i < N; ++i) {
            beads[i].px += 0.5 * beads[i].fx * dt;
            beads[i].py += 0.5 * beads[i].fy * dt;
            beads[i].pz += 0.5 * beads[i].fz * dt;
        }
}

void VerletSimulation::updatePositions() {
        for (int i = 0; i < N; ++i) {
            beads[i].x += beads[i].px/m * dt;
            beads[i].y += beads[i].py/m * dt;
            beads[i].z += beads[i].pz/m * dt;
        }
    }

void VerletSimulation::backintotheBox() {
    for (int i = 0; i < N; ++i) {
        if (beads[i].x > L/2.0) {
            beads[i].x -= L;
        }
        else if (beads[i].x < -L/2.0) {
            beads[i].x += L;
        }

        if (beads[i].y > L/2.0) {
            beads[i].y -= L;
        }
        else if (beads[i].y < -L/2.0) {
            beads[i].y += L;
        }

        if (beads[i].z > L/2.0) {
            beads[i].z -= L;
        }
        else if (beads[i].z < -L/2.0) {
            beads[i].z += L;
        }
    }
}

double VerletSimulation::LennardJones(double deltar){
    return 4*epsilon*(12*std::pow(12, sigma)/std::pow(13, deltar) - 6*std::pow(6, sigma)/std::pow(7, deltar));
}

void VerletSimulation::computeForces() {
    double fx;
    double fy;
    double fz;

    for (int i = 0; i < N; ++i) {
        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        // Compute the forces
        for(int j = i+1; j < N; j++) {
            double rx = (beads[i].x - beads[j].x);

            double ry = (beads[i].y - beads[j].y);

            double rz = (beads[i].z - beads[j].z);

            if (rx > L/2.0) {
            rx -= L;
            }
            else if (rx < -L/2.0) {
            rx += L;
            }
            if (ry > L/2.0) {
            ry -= L;
            }
            else if (ry < -L/2.0) {
            ry += L;
            }
            if (rz > L/2.0) {
            rz -= L;
            }
            else if (rz < -L/2.0) {
            rz += L;
            }

            fx += LennardJones(rx); 
            fy += LennardJones(ry); 
            fz += LennardJones(rz); 

            //Maybe a sign needs to be switched
            beads[i].fx += fx;
            beads[i].fy += fy;
            beads[i].fz += fz;

            beads[i].fx += -fx;
            beads[i].fy += -fy;
            beads[i].fz += -fz;
        }
        

        
    }
}

void VerletSimulation::VerletStep() {
            updateMomenta();
            updatePositions();
            backintotheBox();
            computeForces();
            updateMomenta();
        }


double* VerletSimulation::getPositionsAndMomenta(){
    double* positionAndMomentum = new double[2*3*N];

    for(int i = 0; i < N; i++){
        positionAndMomentum[3*i] = beads[i].x;
        positionAndMomentum[3*i+1] = beads[i].y;
        positionAndMomentum[3*i+2] = beads[i].z;
    }

    for(int i = N; i < 2*N; i++){
        positionAndMomentum[3*i] = beads[i].px;
        positionAndMomentum[3*i+1] = beads[i].py;
        positionAndMomentum[3*i+2] = beads[i].pz;
    }

    return positionAndMomentum;
}


