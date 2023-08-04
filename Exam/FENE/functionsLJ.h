#ifndef MYFILE_H  // Header guard to prevent multiple inclusion
#define MYFILE_H

void create_csv(double** data, std::string name, int rows, int columns);

double MaxwellBoltzmann(double epsilon, double mass);

double boundary(double distance, double size);


#endif  // End of header guard