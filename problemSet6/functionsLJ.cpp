#include<iostream>
#include<cmath>
#include<random>
#include<string>
#include<assert.h>
#include<fstream>
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


std::default_random_engine generator;

double MaxwellBoltzmann(double epsilon, double mass){
    double sigma = sqrt(1.35*epsilon/mass);

    std::normal_distribution<double> distribution(0.0, sigma);

    double temp = distribution(generator);
    return temp;
}

double timeAverage(double* observable, int length, int start){
    double cummalative = 0.0;
    for(int i = start; i < length+start;i++){
        cummalative += *(observable + i);
    }
    return cummalative/length;
}

double fluctuations(double* observable, int length, int start){
    double observable2[length];
    for(int i = 0; i < length; i++){
        observable2[i] = *(observable + i) * *(observable + i);
    }

    double firstMoment = timeAverage(observable, length, start);
    double firstMoment2 = firstMoment*firstMoment2;

    double secondMoment = timeAverage(observable2, length, start);

    return secondMoment - firstMoment2;
}