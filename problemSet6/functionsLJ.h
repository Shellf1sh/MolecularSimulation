#ifndef MYFILE_H  // Header guard to prevent multiple inclusion
#define MYFILE_H

void create_csv(double** data, std::string name, int rows, int columns);

double MaxwellBoltzmann(double epsilon, double mass);

double timeAverage(double* observable, int length, int start);

double fluctuations(double* observable, int length, int start);

#endif  // End of header guard