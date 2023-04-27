#include <stdio.h> 
#include <fstream>
#include <iostream>

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
