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
        VerletSimulation(int N, double n_steps, double m, double dt, double L, double sigma, double epsilon, bool thermo);

        bool periodic = false;

    private:
        std::vector<Beads> beads;
        const int N;
        const double n_steps;
        const double m;
        const double dt;
        const double L;
        const double sigma;
        const double epsilon;
        const bool thermo;
        int step;


    public:
        void initialize(int perDim);
        
        void TestInit();

        void VerletStep(int step);

        void printBeadStates();

        double* getPositionsAndMomenta();
    
        double getKineticEnergy();

        double getPotentialEnergy();

        void setRandomMomenta();

        void printMomenta();

        double getPositionX(int bead);

        void rescaleMomenta(double scaling);

    private:
        void updatePositions();

        void computeForces();

        void updateMomenta();

        void backintotheBox();

        double LennardJones(double deltar);

        double periodicDistance(double d);

};



#endif  // End of header guard
