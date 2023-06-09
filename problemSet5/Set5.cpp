#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

struct Beads {
    double x,y,z; // position in 3 dimensions
    double px, py, pz; // momentum in 3 dimensions
    double fx, fy, fz; // forces in 3 dimensions
};

class VerletSimulation {

public:
    VerletSimulation(int N,double n_steps, double m, double dt, double k, double l, double L)
        : N(N),n_steps(n_steps), m(m), dt(dt), k(k), l(l), L(L)
    {
        beads.resize(N);
    }

    void initialize(double zero, double l, double L,
                    const std::vector<double>& momenta,
                    int pyIndex) {
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

    void VerletStep() {
            updateMomenta();
            updatePositions();
            backintotheBox();
            computeForces();
            updateMomenta();
        }



void printBeadStates() {
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

private:
    std::vector<Beads> beads;
    int N; // Number of beads
    int n_steps; //number of steps
    double dt; // time step
    double k;  // spring constant
    double m; // mass
    double l;
    double L;


    void updatePositions() {
        for (int i = 0; i < N; ++i) {
            beads[i].x += beads[i].px/m * dt;
            beads[i].y += beads[i].py/m * dt;
            beads[i].z += beads[i].pz/m * dt;
        }
    }
void computeForces() {
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
            fy += -k * (beads[i].y - beads[i - 1].y); // y-component of harmonic force
            fz += -k * (beads[i].z - beads[i - 1].z); // z-component of harmonic force
        }
        if (i < N - 1) {

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
            fx += -k * (beads[i].x - beads[i + 1].x); // x-component of harmonic force
            fy += -k * (beads[i].y - beads[i + 1].y); // y-component of harmonic force
            fz += -k * (beads[i].z - beads[i + 1].z); // z-component of harmonic force
        }

        beads[i].fx = fx;
        beads[i].fy = fy;
        beads[i].fz = fz;
    }
}


    void updateMomenta() {
        for (int i = 0; i < N; ++i) {
            beads[i].px += 0.5 * beads[i].fx * dt;
            beads[i].py += 0.5 * beads[i].fy * dt;
            beads[i].pz += 0.5 * beads[i].fz * dt;
        }
    }


    void backintotheBox() {
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



};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    const int N = 2;
    const int n_steps = 100;
    double m = 1;
    double k = 1.0;
    double tau = sqrt(m/k);
    double l = 1.0;
    const double L = 11*l;
    double dt = 0.1;
    double epsilon = 1.;
    double p0 = sqrt(m*epsilon);
    double zero = 0.;

    //std::vector<double> momenta = {-3./7.*p0, 6./7.* p0, -2./7.*p0};
    std::vector<double> momenta = {-1./5.*p0, 2./5.* p0, -1./5.*p0};
    const int pyIndex = 1-1;

    VerletSimulation simulation(N,n_steps,m,  dt, k, l, L);
    simulation.initialize(zero,l,L, momenta, pyIndex);
    simulation.printBeadStates();

    for (int i = 0; i < n_steps; ++i) {
        simulation.VerletStep();  // Perform a single Verlet step
        simulation.printBeadStates();  // Print the bead states after each step
    }
    return 0;
}
