#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include "ConfigFile.h"

using namespace std;
class Exercice4{
private:
    double t, dt, tFin, epsilon;
    double m_L,m_T,m_A, g, d,R_T, R_L, r_0;
    double Omega;
    valarray<valarray<double>> y=valarray<valarray<double>>(valarray<double>(0.0, 4),3);
    int sampling;
    unsigned int nsteps;
    int last;
};

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
