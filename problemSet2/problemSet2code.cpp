#include<iostream>
#include<cmath>
#include<fstream>
#include<cassert>

class Dumbbell{
    public: 
        double m, k, delta_t, energy; //When creating an instance of this class remember to define these variables
        double angular_momentum[6] = {0, 0, 0, 0, 0, 0};
        
        double r[6] = {0, 0, 0, 0, 0, 0}; //The elements are [x1, y1, z1, x2, y2, z2]
        double p[6] = {0, 0, 0, 0, 0, 0}; //The elements are [px1, py1, pz1, px2, py2, pz2]
        double F[6] = {0, 0, 0, 0, 0, 0}; //The elements are [fx1, fy1, fz1, fx2, fy2, fz2]

        void setInitialConditions(double initialPosition[], double initialMomentum[], int sizePosition = 6, int sizeMomentum = 6){
            for(int i = 0; i < 6; i++){
                r[i] = initialPosition[i];
                p[i] = initialMomentum[i];
            }
        }

        void updatePosition(){
            for(int i = 0; i < 6; i++){
                r[i] = r[i] + p[i]*delta_t/m;
            }
        }

        void updateMomentum(){
            for(int i = 0; i < 6; i++){
                p[i] = p[i] + F[i]*delta_t/2;
            }
        }

        void calculateForce(){
            for(int i = 0; i < 3; i++){
                F[i] = -k*(r[i] - r[i+3]);
                F[i+3] = k*(r[i] - r[i+3]);
            }
        }

        void updateEnergy(){
            double kinetic_energy = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+p[4]*p[4]+p[5]*p[5])/m;
            double potnetial_energy = 0.5*k*((r[0]-r[3])*(r[0]-r[3]) + (r[1]-r[4])*(r[1]-r[4]) + (r[2]-r[5])*(r[2]-r[5]));

            energy = kinetic_energy + potnetial_energy;
        }

        void updateAngularMomentum(){
            angular_momentum[0] = (r[1]*p[2] - r[2]*p[1]) + (r[4]*p[5] - r[5]*p[4]); 
            angular_momentum[1] = (r[2]*p[0] - r[0]*p[2]) + (r[5]*p[3] - r[3]*p[5]);
            angular_momentum[2] = (r[0]*p[1] - r[1]*p[0]) + (r[3]*p[4] - r[4]*p[3]);

            /*angular_momentum[3] = r[4]*p[5] - r[5]*p[4]; 
            angular_momentum[4] = r[5]*p[3] - r[3]*p[5];
            angular_momentum[5] = r[3]*p[4] - r[4]*p[3];*/
        }

        void verletStep(){
            updateMomentum(); //First step of the Verlet-velocity algorithm
            updatePosition(); //Second step
            calculateForce(); //Third step
            updateMomentum(); //Fourth step
        }

};


int main(){

    Dumbbell dumbbell; // Initialise the dumbbell object

    //We set the charactaristic length scale a = 1
    dumbbell.m = 5.0; //Define constants
    dumbbell.k = 0.4;
    dumbbell.delta_t = 1e-2;

    double initPosition[6] = {1.0/3.0, 0.0, 0.0, -1.0/3.0, 0.0, 0.0}; //Setting initial conditions
    double initMomentum[6] = {-6.0/7.0*sqrt(dumbbell.m*dumbbell.k), 3.0/7.0*sqrt(dumbbell.m*dumbbell.k), -2.0/7.0*sqrt(dumbbell.m*dumbbell.k), 
                                6.0/7.0*sqrt(dumbbell.m*dumbbell.k), -3.0/7.0*sqrt(dumbbell.m*dumbbell.k), 2.0/7.0*sqrt(dumbbell.m*dumbbell.k)};
    dumbbell.setInitialConditions(initPosition, initMomentum);

    int T = 10*sqrt(dumbbell.m/dumbbell.k);
    int N = T/dumbbell.delta_t;

    double trajectories[N][17]; 
    //The columns correspond to 6 positions coordinates, 6 momentum coordinates
    //1 energy value and 3 angular momentum coordinates and 1 time variable
    //[rx1, ry1, rz1, rx2, ry2, rz2, px1, py1, pz1, px2, py2, pz2, 
    // Jx, Jy, Jz, E, t]

    dumbbell.calculateForce(); //Initial force calculation

    for(int i = 0; i < N; i++){
        dumbbell.verletStep();
        
        for(int j = 0; j < 6; j++){
            trajectories[i][j] = dumbbell.r[j];
        }
        for(int j = 6; j < 12; j++){
            trajectories[i][j] = dumbbell.p[j];
        }
        dumbbell.updateAngularMomentum();
        for(int j = 12; j < 15; j++){
            trajectories[i][j] = dumbbell.angular_momentum[j];
        }
        
        dumbbell.updateEnergy();
        trajectories[i][15] = dumbbell.energy;

        trajectories[i][16] = dumbbell.delta_t * i;

    }

    //The following block of code writes the trajectories to 
    //a .csv file so they can be plotted in python
    std::ofstream write_output("trajectory1.csv");
    assert(write_output.is_open());

    for(int i = 0; i < N; i++){
        for(int j = 0; j < 16; j++){
            write_output << trajectories[i][j] << ",";
        }
        write_output << trajectories[i][16]; //This line is here to avoid adding a , at the end of the .csv file
        write_output << "\n";
    }
    write_output.close();



    std::cout << sizeof(trajectories[0]) / sizeof(double) << std::endl;

    std::cout << "Done" << std::endl;
    return 0;
}
