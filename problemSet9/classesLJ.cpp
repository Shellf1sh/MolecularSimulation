#include<iostream>
#include<iomanip>
#include<cmath>
#include "classesLJ.h"
#include "functionsLJ.h"



VerletSimulation::VerletSimulation(int N,double n_steps, double m, double dt, double L, double sigma, double epsilon, bool thermo)
        : N(N),n_steps(n_steps), m(m), dt(dt), L(L), sigma(sigma), epsilon(epsilon), thermo(thermo)
    {
        beads.resize(N);
    }

//This is the thermostat
void VerletSimulation::setRandomMomenta(){
    for(int i = 0; i < N; i++){
            beads[i].px = MaxwellBoltzmann(epsilon, m);
            beads[i].py = MaxwellBoltzmann(epsilon, m);
            beads[i].pz = MaxwellBoltzmann(epsilon, m);
        }
}

void VerletSimulation::printMomenta(){
    std::cout << "Momneta=(" << beads[0].px << ", " << beads[0].py << ", " << beads[0].pz << ")" << std::endl;

}

void VerletSimulation::initialize(int perDim) {
        //Set the particles in a lattice pattern to avoid overlap
        int index = 0;
        for(int i = 0; i < perDim; i++){
            for(int j = 0; j < perDim; j++){
                for(int k = 0; k < perDim; k++){
                    beads[index].x = i * 1.5*sigma;
                    beads[index].y = j * 1.5*sigma;
                    beads[index].z = k * 1.5*sigma;
                    index++;
                }
            }
        }
        setRandomMomenta(); //Give the particle random Maxwell-Boltzmann distributed momenta
    }

void VerletSimulation::TestInit(){
    if(N == 2){
        beads[0].x = -1.3*sigma;
        beads[0].y = 0.0;
        beads[0].z = 0.0;
        beads[1].x = 1.3*sigma;
        beads[1].y = 0.0;
        beads[1].z = 0.0;

        beads[0].px = 0.0;
        beads[0].py = 0.0;
        beads[0].pz = 0.0;
        beads[1].px = 0.0;
        beads[1].py = 0.0;
        beads[1].pz = 0.0;
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
        beads[i].x = periodicDistance(beads[i].x);
        beads[i].y = periodicDistance(beads[i].y);
        beads[i].z = periodicDistance(beads[i].z);
    }
}


double VerletSimulation::periodicDistance(double d){
    if(d > L/2.0){
        d -= L;
    }else if(d < -L/2.0){
        d += L;
    }
    return d;
}

void VerletSimulation::computeForces() {
    double distance;
    double force_magnitude;

    //Reset all the forces
    for(int i = 0; i < N; i++){
        beads[i].fx = 0.0;
        beads[i].fy = 0.0;
        beads[i].fz = 0.0;
    }

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    for (int i = 0; i < N; ++i) {

        // Compute the forces
        for(int j = i+1; j < N; j++) {
            double rx = (beads[i].x - beads[j].x);
            double ry = (beads[i].y - beads[j].y);
            double rz = (beads[i].z - beads[j].z);

            //Checks whether there's a mirror particles closer
            rx = periodicDistance(rx);
            ry = periodicDistance(ry);
            rz = periodicDistance(rz);

            //std::cout << " after mirror check" << rx << std::endl;

            distance = sqrt(rx*rx + ry*ry + rz*rz);
            
            force_magnitude = lennardJones(distance, epsilon, sigma);

            fx = -force_magnitude * rx/distance; 
            fy = -force_magnitude * ry/distance; 
            fz = -force_magnitude * rz/distance; 

            //Maybe a sign needs to be switched
            beads[i].fx += -fx;
            beads[i].fy += -fy;
            beads[i].fz += -fz;

            beads[j].fx += fx;
            beads[j].fy += fy;
            beads[j].fz += fz;
        }
    }
}

void VerletSimulation::VerletStep(int step) {
    updateMomenta();
    if(thermo && ((step % 50) == 0)){
        setRandomMomenta();
    }
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

        positionAndMomentum[3*i + 3*N] = beads[i].px;
        positionAndMomentum[3*i+1 + 3*N] = beads[i].py;
        positionAndMomentum[3*i+2 + 3*N] = beads[i].pz;
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

    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            distance = sqrt((beads[i].x-beads[j].x)*(beads[i].x-beads[j].x) + (beads[i].y-beads[j].y)*(beads[i].y-beads[j].y) + (beads[i].z-beads[j].z)*(beads[i].z-beads[j].z));

            Epot += 4*epsilon*(std::pow(sigma/distance, 12) - std::pow(sigma/distance, 6));
        }
    }
    return Epot;
}

void VerletSimulation::rescaleMomenta(double scaling){
    for( int i = 0; i < N; i++){
        beads[i].px = beads[i].px/sqrt(scaling);
        beads[i].py = beads[i].py/sqrt(scaling);
        beads[i].pz = beads[i].pz/sqrt(scaling);
    }
}


