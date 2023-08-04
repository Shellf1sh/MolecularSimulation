#ifndef MYCLASS_H  // Header guard to prevent multiple inclusion
#define MYCLASS_H

#include<vector>

struct Beads {
    double x,y,z; // position in 3 dimensions
    double px, py, pz; // momentum in 3 dimensions
    double fx, fy, fz; // forces in 3 dimensions
};


class VerletSimulation{   
    public:
        VerletSimulation(int N, double n_steps, double m, double dt, double L, double sigma, double epsilon, double gamma, double kT);

    private:
        std::vector<Beads> beads;
        const int N;
        const double n_steps;
        const double m;
        const double dt;
        const double L;
        const double sigma;
        const double epsilon;
        const double gamma;
        const double kT;

    public:
        void initialize(int particleNumber, double boxSize);
                
        void VerletStep(bool hamil);

        void printBeadStates();

        double* getPositionsAndMomenta();

        double getKineticEnergy();

        double getPotentialEnergy();

        double* getForces();

        void testInit();
    
    private:
        void updatePositions(bool hamil);

        void computeForces();

        void updateMomenta();

        void updateRandom();

        void backintotheBox();

        double LennardJones(double deltar);
};



#endif  // End of header guard
