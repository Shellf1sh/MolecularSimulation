#include<iostream>
#include<cmath>
#include<random>
#include<string>
#include<assert.h>
#include<fstream>
#include<ctime>

#include"functionsLJ.h"


void create_csv(double** data, std::string name, int rows, int columns){
    std::ofstream write_output(name);

    assert(write_output.is_open());

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns-1; j++){
            //write_output << *(data + i * columns + j) << ",";
            write_output << data[i][j] << ",";
        }
        //write_output << *(data + (i+1) * columns - 1); //This line is here to avoid adding a , at the end of the .csv file
        write_output << data[i][columns - 1];
        write_output << "\n";
    }
    write_output.close();

    std::cout << "CSV created: " << name << std::endl;
}

std::seed_seq seed1{std::time(nullptr)};
std::default_random_engine generator(40);

double MaxwellBoltzmann(double kT, double mass){
    double sigma = sqrt(kT/mass);

    std::normal_distribution<double> distribution(0.0, sigma);

    double temp = distribution(generator);
    return temp;
}

double boundary(double distance, double size){
    if(distance > size/2){
        distance -= size;
    }else if(distance < -size/2){
        distance += size;
    }
    return distance;
}

double random_number(double L){
    //std::seed_seq seed{std::time(nullptr)};

    //std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> uniform(-L/2, L/2);
    return uniform(generator);

}

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2, double L) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;

    dx = boundary(dx, L);
    dy = boundary(dy, L);
    dz = boundary(dz, L);

    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    return distance;
}

double LJ_potential(double deltar, double epsilon, double sigma){
    return 4*epsilon*(std::pow(sigma/deltar, 12) - std::pow(sigma/deltar, 6));
}