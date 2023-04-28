#include<iostream>
#include<cmath>
#include<fstream>
#include<cassert>

class Dumbbell{
    public: 
        double m, k, delta_t, energy; //When creating an instance of this class remember to define these variables
        double angular_momentum[3] = {0, 0, 0};
        
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
            angular_momentum[0] = r[1]*p[2] - r[2]*p[1]; //TODO: Check that this is correct
            angular_momentum[1] = r[2]*p[0] - r[0]*p[2];
            angular_momentum[2] = r[0]*p[1] - r[1]*p[0];
        }

};


int main(){

    Dumbbell dumbbell; // Initialise the dumbbell object

    //We set the charactaristic length scale a=1
    dumbbell.m = 5.0; //Define constants
    dumbbell.k = 0.4;
    dumbbell.delta_t = 1e-2;

    double initPosition[6] = {0.5, 0, 0, -0.5, 0, 0}; //Setting initial conditions
    double initMomentum[6] = {sqrt(dumbbell.m*dumbbell.k), 0.7*sqrt(dumbbell.m*dumbbell.k), 0, 0, -0.5*sqrt(dumbbell.m*dumbbell.k), 0};
    dumbbell.setInitialConditions(initPosition, initMomentum);



    int T = 10*sqrt(dumbbell.m/dumbbell.k);
    int N = T/dumbbell.delta_t;

    double trajectories[N][10]; 
    //The columns correspond to 6 positions coordinates 
    //1 energy value and 3 angular momentum coordinates

    dumbbell.calculateForce();

    for(int i = 0; i < N; i++){
        dumbbell.updateMomentum(); //First step of the Verlet-velocity algorithm
        dumbbell.updatePosition(); //Second step
        dumbbell.calculateForce(); //Third step
        dumbbell.updateMomentum(); //Fourth step
        
        for(int j = 0; j < 6; j++){
            trajectories[i][j] = dumbbell.r[j];
        }

        dumbbell.updateEnergy();
        trajectories[i][6] = dumbbell.energy;

        dumbbell.updateAngularMomentum();
        for(int j = 7; j < 10; j++){
            trajectories[i][j] = dumbbell.angular_momentum[j];
        }

    }

    //The following block of code writes the trajectories to 
    //a .csv file so they can be plotted in python
    std::ofstream write_output("trajectory1.csv");
    assert(write_output.is_open());

    for(int i = 0; i < N; i++){
        for(int j = 0; j < 10; j++){
            write_output << trajectories[i][j] << ",";
        }
        write_output << "\n";

        /*write_output << trajectories[i][0] << "," << trajectories[i][1] << "," << trajectories[i][2] << "," << trajectories[i][3] 
        << "," << trajectories[i][4] << "," << trajectories[i][5] << "," << trajectories[i][6] <<
        "," << trajectories[i][7] << "," << trajectories[i][8] << "," << trajectories[i][9] <<"\n";*/
    }
    write_output.close();



    std::cout << "Done" << std::endl;
    return 0;
}
