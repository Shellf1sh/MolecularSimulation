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
        VerletSimulation(int N, double n_steps, double m, double dt, double k, double l, double L, bool loop);

    private:
        std::vector<Beads> beads;
        const int N;
        const double n_steps;
        const double m;
        const double dt;
        const double k;
        const double l;
        const double L;
        const bool loop;

    public:
        void initialize(double zero, double l, double L,
                    const std::vector<double>& momenta,
                    int pyIndex);
                
        void VerletStep();

        void printBeadStates();

        double* getPositionsAndMomenta();
    
    private:
        void updatePositions();

        void computeForces();

        void updateMomenta();

        void backintotheBox();
};



#endif  // End of header guard
