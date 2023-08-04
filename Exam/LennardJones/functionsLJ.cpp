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
