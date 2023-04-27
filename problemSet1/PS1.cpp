#include <stdio.h> 
#include <fstream>
#include <iostream>


//Initial conditions
int t_0 = 0;
int t_f = 10;
double delta_t_small = 10e-3;
double delta_t_large = 10e-1;


int otherFunc(){
    printf("Test2");
}

int number = 30;

int main(){


    printf("Test\n");

    otherFunc();

    std::cout << number << std::endl;

    return 0;
}
