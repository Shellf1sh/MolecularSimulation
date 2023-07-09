#ifndef MYFILE_H  // Header guard to prevent multiple inclusion
#define MYFILE_H

struct Beads;

void create_csv(double** data, std::string name, int rows, int columns);

double MaxwellBoltzmann(double epsilon, double mass);

double timeAverage(double* observable, int length, int start);

double fluctuations(double* observable, int length, int start);

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);

double lennardJones(double deltar, double epsilon, double sigma);

double LJ_potential(double deltar, double epsilon, double sigma);

double random_pos(double length);
#endif  // End of header guard