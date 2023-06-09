#include<iostream>
#include<iomanip>
#include<cmath>
#include "classes.h"



VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt, double k, double l, double L, bool loop)
        : N(N),n_steps(n_steps), m(m), dt(dt), k(k), l(l), L(L), loop(loop)
    {
        beads.resize(N);
    }

void VerletSimulation::initialize(double zero, double l, double L, const std::vector<double>& momenta, int pyIndex) {
        for (int i = 0; i < N; ++i) {
            beads[i].x = (i+1)*l-L/2;
            beads[i].y = zero;
            beads[i].z = zero;
            if (i == pyIndex) {
                beads[i].px = momenta[0];  // x-component of momentum
                beads[i].py = momenta[1];  // y-component of momentum
                beads[i].pz = momenta[2];  // z-component of momentum
            } else {
                beads[i].px = zero;
                beads[i].py = zero;
                beads[i].pz = zero;
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

void VerletSimulation::computeForces() {

    for (int i = 0; i < N; ++i) {
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;

        // Compute the forces for the harmonic bonds
        if (i > 0) {
            double rx = (beads[i].x - beads[i - 1].x);
            double ry = (beads[i].y - beads[i - 1].y);
            double rz = (beads[i].z - beads[i - 1].z);

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


            fx += -k * rx; // x-component of harmonic force
            fy += -k * ry; // y-component of harmonic force
            fz += -k * rz; // z-component of harmonic force
        }
        if (i < N - 1) {

            double rx = (beads[i].x - beads[i + 1].x);
            double ry = (beads[i].y - beads[i + 1].y);
            double rz = (beads[i].z - beads[i + 1].z);

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
            fx += -k * rx; // x-component of harmonic force
            fy += -k * ry; // y-component of harmonic force
            fz += -k * rz; // z-component of harmonic force
        }

        if(loop && (i == N-1)){
            double rx = (beads[i].x - beads[1].x);
            double ry = (beads[i].y - beads[1].y);
            double rz = (beads[i].z - beads[1].z);

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
            fx += -k * rx; // x-component of harmonic force
            fy += -k * ry; // y-component of harmonic force
            fz += -k * rz; // z-component of harmonic force

        }

        if(loop && (i == 0)){
            double rx = (beads[i].x - beads[N-1].x);
            double ry = (beads[i].y - beads[N-1].y);
            double rz = (beads[i].z - beads[N-1].z);

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
            fx += -k * rx; // x-component of harmonic force
            fy += -k * ry; // y-component of harmonic force
            fz += -k * rz; // z-component of harmonic force
        }

        beads[i].fx = fx;
        beads[i].fy = fy;
        beads[i].fz = fz;
    }
}

void VerletSimulation::VerletStep() {
            updateMomenta();
            updatePositions();
            backintotheBox();
            computeForces();
            updateMomenta();
        }

void VerletSimulation::printBeadStates() {
    // Set the width for each field
    const int fieldWidth = 10;

    // Print the header
    std::cout << std::left << std::setw(fieldWidth) << "Bead"
              << std::right << std::setw(fieldWidth) << "Position"
              << std::right << std::setw(fieldWidth) << "Momentum" << '\n';

    // Print the bead states
    for (int i = 0; i < N; ++i) {
        std::cout << std::left << std::setw(fieldWidth) << i
                  << std::right << std::setw(fieldWidth) << "(" << beads[i].x << ", " << beads[i].y << ", " << beads[i].z << ")"
                  << std::right << std::setw(fieldWidth) << "(" << beads[i].px << ", " << beads[i].py << ", " << beads[i].pz << ")\n";
    }
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




    //std::cout << *(positionAndMomentum + 3*N) << std::endl;

    return positionAndMomentum;
}


