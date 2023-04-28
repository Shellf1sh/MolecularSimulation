#include<iostream>
#include<cmath>

class Dumbbell{
    public: 
        double m, k, delta_t, energy; //When creating an instance of this class remember to define these variables
        double angular_momentum[3];
        
        double r[6] = {0, 0, 0, 0, 0, 0}; //The elements are [x1, y1, z1, x2, y2, z2]
        double p[6] = {0, 0, 0, 0, 0, 0}; //The elements are [px1, py1, pz1, px2, py2, pz2]
        double F[6] = {0, 0, 0, 0, 0, 0}; //The elements are [fx1, fy1, fz1, fx2, fy2, fz2]

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
            angular_momentum[0] = r[1]*p[2] - r[2]*p[1]; //TODO: Check that this is correct
            angular_momentum[1] = r[2]*p[0] - r[0]*p[2];
            angular_momentum[2] = r[0]*p[1] - r[1]*p[0];
        }

};


int main(){

    Dumbbell dumbbell; // Initialise the dumbbell object

    dumbbell.m = 5.0; //Define constants
    dumbbell.k = 0.4;
    dumbbell.delta_t = 1e-2;

    dumbbell.r[0] = 0.5; // Initial conditions
    dumbbell.r[3] = -0.5;
    dumbbell.p[0] = 0.5*sqrt(dumbbell.m*dumbbell.k);
    dumbbell.p[3] = -0.5*sqrt(dumbbell.m*dumbbell.k);


    int T = 10;
    int N = T/dumbbell.delta_t;

    double trajectories[N][6];

    dumbbell.calculateForce();

    for(int i = 0; i < N; i++){
        dumbbell.updateMomentum(); //First step of the Verlet-velocity algorithm
        dumbbell.updatePosition(); //Second step
        dumbbell.calculateForce(); //Third step
        dumbbell.updateMomentum(); //Fourth step
        
        for(int j = 0; j < 6; j++){
            trajectories[i][j] = dumbbell.r[j];
        }


        if(i%100 == 0){
            std::cout << dumbbell.r[0] << "\n"; 
        }

    }

    //TODO save the trajectories array as a .csv file for plotting in python

    std::cout << "Done" << std::endl;
    return 0;
}
