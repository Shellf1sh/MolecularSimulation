#include<iostream>
#include<cmath>
#include<fstream>
#include<cassert>
#include<string>

class Dumbbell{
    public: 
        double m, k, delta_t, energy; //When creating an instance of this class remember to define these variables

        double angular_momentum[6] = {0, 0, 0, 0, 0, 0};
        bool FENE = false;
        
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
            if(!FENE){
                for(int i = 0; i < 3; i++){
                    F[i] = -k*(r[i] - r[i+3]);
                    F[i+3] = k*(r[i] - r[i+3]);          
                }
            }else{
                for(int i = 0; i < 3; i++){
                    F[i] = -k*(r[i] - r[i+3])/(2.0*(1.0-abs(r[i]-r[i+3])*abs(r[i]-r[i+3])));
                    F[i+3] = k*(r[i] - r[i+3])/(2.0*(1.0-abs(r[i]-r[i+3])*abs(r[i]-r[i+3])));      
                }
            }
        }

        void updateEnergy(){
            double kinetic_energy = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+p[4]*p[4]+p[5]*p[5])/(2.0*m);
            double potential_energy;
            if(!FENE){
                potential_energy = 0.5*k*((r[0]-r[3])*(r[0]-r[3]) + (r[1]-r[4])*(r[1]-r[4]) + (r[2]-r[5])*(r[2]-r[5]));
            }else{
                potential_energy = -0.5*k*log(1-((r[0]-r[3])*(r[0]-r[3]) + (r[1]-r[4])*(r[1]-r[4]) + (r[2]-r[5])*(r[2]-r[5])));
            }
            energy = kinetic_energy + potential_energy;
        }

        void updateAngularMomentum(){
            angular_momentum[0] = (r[1]*p[2] - r[2]*p[1]) + (r[4]*p[5] - r[5]*p[4]); 
            angular_momentum[1] = (r[2]*p[0] - r[0]*p[2]) + (r[5]*p[3] - r[3]*p[5]);
            angular_momentum[2] = (r[0]*p[1] - r[1]*p[0]) + (r[3]*p[4] - r[4]*p[3]);
        }

        void verletStep(){
            updateMomentum(); //First step of the Verlet-velocity algorithm
            updatePosition(); //Second step
            calculateForce(); //Third step
            updateMomentum(); //Fourth step
        }
};



void runSimulation(double dataArray[][17], Dumbbell dumb, int N){
    for(int i = 0; i < N; i++){
        dumb.verletStep();
        
        for(int j = 0; j < 6; j++){
            dataArray[i][j] = dumb.r[j];
        }
        for(int j = 6; j < 12; j++){
            dataArray[i][j] = dumb.p[j];
        }

        dumb.updateAngularMomentum();
        for(int j = 12; j < 15; j++){
            dataArray[i][j] = dumb.angular_momentum[j];
        }
        
        dumb.updateEnergy();
        dataArray[i][15] = dumb.energy;

        dataArray[i][16] = dumb.delta_t * i;
    }
}

//The following block of code writes the trajectories to 
//a .csv file so they can be plotted in python
void create_csv(double data[][17], std::string file_name, int columns, int rows){

    std::ofstream write_output(file_name);
    assert(write_output.is_open());

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            write_output << data[i][j] << ",";
        }
        write_output << data[i][16]; //This line is here to avoid adding a , at the end of the .csv file and creating an extra column
        write_output << "\n";
    }
    write_output.close();
}


int main(){

    Dumbbell dumbbell; // Initialise the dumbbell object

    //We set the charactaristic length scale a = 1
    dumbbell.m = 5.0; //Define constants
    dumbbell.k = 0.4;
    dumbbell.delta_t = 1e-2 * sqrt(dumbbell.m*dumbbell.k);

    double factor = 0.4;

    double initPosition[6] = {1.0/3.0, 0.0, 0.0, -1.0/3.0, 0.0, 0.0}; //Setting initial conditions
    double initMomentum[6] = {-6.0/7.0*sqrt(dumbbell.m*dumbbell.k) *factor, 3.0/7.0*sqrt(dumbbell.m*dumbbell.k) *factor, -2.0/7.0*sqrt(dumbbell.m*dumbbell.k ) *factor, 
                                6.0/7.0*sqrt(dumbbell.m*dumbbell.k) *factor, -3.0/7.0*sqrt(dumbbell.m*dumbbell.k) *factor, 2.0/7.0*sqrt(dumbbell.m*dumbbell.k) *factor};
    dumbbell.setInitialConditions(initPosition, initMomentum);

    int T = 10*sqrt(dumbbell.m/dumbbell.k);
    int N = T/dumbbell.delta_t;

    dumbbell.calculateForce(); //Initial force calculation

    
    double trajectories[N][17]; 
    //The columns correspond to 6 positions coordinates, 6 momentum coordinates
    //1 energy value and 3 angular momentum coordinates and 1 time variable
    //[rx1, ry1, rz1, rx2, ry2, rz2, px1, py1, pz1, px2, py2, pz2, 
    // Jx, Jy, Jz, E, t]

    runSimulation(trajectories, dumbbell, N);

    create_csv(trajectories, "trajectory1.csv", 17, N);


    std::cout << "First simulation done. \n";


    //Now we run the simulation, but the time step is 4 times bigger
    Dumbbell dumbbell_big_step;

    //We set the charactaristic length scale a = 1
    dumbbell_big_step.m = 5.0; //Define constants
    dumbbell_big_step.k = 0.4;
    dumbbell_big_step.delta_t = 4 * 1e-2 * sqrt(dumbbell_big_step.m*dumbbell_big_step.k); //This is where the difference is

    int N_big_step = T/dumbbell_big_step.delta_t;//Same T as before
    double trajectories_big_step[N_big_step][17]; //Initialise array for simulation data

    dumbbell_big_step.setInitialConditions(initPosition, initMomentum); //Same initial conditions as before
    dumbbell_big_step.calculateForce(); //Initial force calculation

    runSimulation(trajectories_big_step, dumbbell_big_step, N_big_step);
    create_csv(trajectories_big_step, "trajectory_big.csv", 17, N_big_step);

    std::cout << "Second simulation done. \n";

    //End of problem a

    Dumbbell dumbbell_FENE; //Creating a new dumbbell object and switching to the FENE potential
    dumbbell_FENE.FENE = true;

    //We set the charactaristic length scale a = 1
    dumbbell_FENE.m = 5.0; //Define constants
    dumbbell_FENE.k = 0.4;

    dumbbell_FENE.delta_t = 1e-2 * sqrt(dumbbell_FENE.m*dumbbell_FENE.k);

    int N_FENE = T/dumbbell_FENE.delta_t;

    dumbbell_FENE.setInitialConditions(initPosition, initMomentum);
    dumbbell_FENE.calculateForce(); //Initial force calculation
    double trajectories_FENE[N_FENE][17]; 
    runSimulation(trajectories_FENE, dumbbell_FENE, N_FENE);
    create_csv(trajectories_FENE, "trajectory_FENE.csv", 17, N_FENE);

    std::cout << "First FENE simulation done. \n";

    //We redefine the time step so we can run the simulation again
    dumbbell_FENE.delta_t = 4 * 1e-2 * sqrt(dumbbell_FENE.m*dumbbell_FENE.k); //Rerun but with the big delta t
    dumbbell_FENE.setInitialConditions(initPosition, initMomentum);
    dumbbell_FENE.calculateForce(); //Initial force calculation
    double trajectories_FENE_big[N_big_step][17]; 
    runSimulation(trajectories_FENE_big, dumbbell_FENE, N_big_step);
    create_csv(trajectories_FENE_big, "trajectory_FENE_big.csv", 17, N_big_step);

    std::cout << "Big step FENE simulation done. \n";

    std::cout << N << " and " << N_big_step << " \n";

    std::cout << "Done" << std::endl;
    return 0;
}
