#ifndef MYCLASS_H  // Header guard to prevent multiple inclusion
#define MYCLASS_H

#include<vector>

struct Beads{
    double x,y,z; // position in 3 dimensions
    double px, py, pz; // momentum in 3 dimensions
    double fx, fy, fz; // forces in 3 dimensions
};


class VerletSimulation{   
    public:
        VerletSimulation(int N, double n_steps, double m, double dt, double epsilon, double gamma, bool force);

        bool periodic = false;

    private:
        std::vector<Beads> beads;
        const int N;
        const double n_steps;
        const double m;
        const double dt;
        const double epsilon;
        int step;
        const double gamma;

        const double l = 1/gamma * sqrt(epsilon/m);
        const double omega = gamma/3;
        const bool force;

    public:
        void initialize(int perDim);
        
        void zeroInit();

        void VerletStep(int step);

        void printBeadStates();

        double* getPositionsAndMomenta();
    
        double* getKineticEnergy();

        double* getPotentialEnergy();

        void setRandomMomenta();

        void printMomenta();

        void rescaleMomenta(double scaling);

    private:
        void updatePositions();

        void computeForces();

        void updateMomenta();

        void updateRandom();

        void backintotheBox();

        double periodicDistance(double d);

};



#endif  // End of header guard
