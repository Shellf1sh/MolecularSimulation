#include<iostream>
#include<cmath>

class Dumbbell{
    public: double m, k, delta_t; //When creating an instance of this class remember to define these variables
    
    public: double x1, y1, z1, x2, y2, z2; //Initialising the 6 position coordinates
    public: double px1, py1, pz1, px2, py2, pz2; //Initialising the 6 position coordinates
    public: double fx1, fy1, fz1, fx2, fy2, fz2; //Initialising the 6 force coordinates


    public: void updatePosition(){
        x1 = x1 + px1*delta_t/m;
        y1 = y1 + py1*delta_t/m;
        z1 = z1 + pz1*delta_t/m;
        x2 = x2 + px2*delta_t/m;
        y2 = y2 + py2*delta_t/m;
        z2 = z2 + pz2*delta_t/m;
    }

    public: void updateMomentum(){
        px1 = px1 + fx1*delta_t/2;
        py1 = py1 + fy1*delta_t/2;
        pz1 = pz1 + fz1*delta_t/2;
        px2 = px2 + fx2*delta_t/2;
        py2 = py2 + fy2*delta_t/2;
        pz2 = pz2 + fz2*delta_t/2;
    }

    public: void calculateForce(){
        fx1 = -k*(x1-x2);
        fy1 = -k*(y1-y2);
        fz1 = -k*(z1-z2);
        fx2 = k*(x1-x2);
        fy2 = k*(y1-y2);
        fz2 = k*(z1-z2);
    }

};


int main(){

    Dumbbell dumbbell; // Initialise the dumbbell object

    dumbbell.m = 5.0; //Define constants
    dumbbell.k = 0.4;
    dumbbell.delta_t = 1e-2;

    dumbbell.x1 = 0.5; //Initial conditions: position 
    dumbbell.y1 = 0;
    dumbbell.z1 = 0;
    dumbbell.x2 = -0.5;
    dumbbell.y2 = 0;
    dumbbell.z2 = 0;

    dumbbell.px1 = 0;
    dumbbell.py1 = 0.5*sqrt(dumbbell.m*dumbbell.k);
    dumbbell.pz1 = 0;
    dumbbell.px2 = 0;
    dumbbell.py2 = -0.5*sqrt(dumbbell.m*dumbbell.k);
    dumbbell.pz2 = 0;


    int T = 10;
    int N = T/dumbbell.delta_t;

    double trajectories[6][N];

    for(int i = 0; i < N; i++){
        dumbbell.updateMomentum(); //First step of the Verlet-velocity algorithm
        dumbbell.updatePosition(); //Second step
        dumbbell.calculateForce(); //Third step
        dumbbell.updateMomentum(); //Fourth step

        trajectories[0][i] = dumbbell.x1;
        trajectories[1][i] = dumbbell.y1;
        trajectories[2][i] = dumbbell.z1;
        trajectories[3][i] = dumbbell.x2;
        trajectories[4][i] = dumbbell.y2;
        trajectories[5][i] = dumbbell.z2;

        if(i%100 == 0){
            std::cout << dumbbell.x1 << "\n"; 
        }

    }


    std::cout << "Done" << std::endl;
    return 0;
}
