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
        VerletSimulation(int N, double n_steps, double m, double dt, double L, double gamma, double kT, double rMax, double K, bool loop);

    private:
        std::vector<Beads> beads;
        const int N;
        const double n_steps;
        const double m;
        const double dt;
        double L;
        const double gamma;
        const double kT;
        const double rMax;
        const double K;
        const bool loop;

    public:
        void initialize(bool loop);
                
        void VerletStep(bool hamil);

        void printBeadStates();

        double* getPositionsAndMomenta();

        double getKineticEnergy();

        double getPotentialEnergy();

        void testInit();

        double FENEPotential(double deltar);

        double FENEForce(double deltar);
    
    private:
        void updatePositions(bool hamil);

        void computeForces();

        void updateMomenta();

        void updateRandom();

        void backintotheBox();

        double LennardJones(double deltar);
};



#endif  // End of header guard
