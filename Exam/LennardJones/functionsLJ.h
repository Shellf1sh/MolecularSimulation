#ifndef MYFILE_H  // Header guard to prevent multiple inclusion
#define MYFILE_H

void create_csv(double** data, std::string name, int rows, int columns);

double MaxwellBoltzmann(double epsilon, double mass);

double boundary(double distance, double size);

double random_number(double L);

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2, double L);

double LJ_potential(double deltar, double epsilon, double sigma);

#endif  // End of header guard