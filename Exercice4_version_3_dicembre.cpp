#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include "ConfigFile.h"

using namespace std;
class Exercice4 {
private:
    double t, dt, tFin, epsilon;
    double m_L, m_T, m_A;
    double G=6.67430* pow(10,-11);
    double d=160;
    double R_T=6378.1;
    double D=384748; //distance terre lune
    valarray<double> M=valarray<double>(0.0,3);
    double Omega;
    valarray<valarray<double>> y = valarray<valarray<double>>(valarray<double>(0.0, 4),
                                                              3); //dans l'ordre asteroide terre lune ;
    int sampling;
    unsigned int nsteps;
    int last;
    ofstream *outputFile;

    void printOut(bool force) {
        if ((!force && last >= sampling) || (force && last != 1)) {
            double emec =  0.5*M[0]*(y[0][2]*y[0][2] + y[0][3]*y[0][3]) - G*M[0]*M[2]/(sqrt((y[0][0]-y[2][0])*(y[0][0]-y[2][0]) + (y[0][1]-y[2][1])*(y[0][1]-y[2][1])));
            *outputFile << t << " ";
            for (unsigned int i(0); i < 4; ++i) {
                for (unsigned int j(0); j < 3; ++j) {
                    cout << y[i][j] << " ";
                }
                cout << emec << endl;

                last = 1;
            }
        } else {
            last++;
        }
    }

    valarray<valarray<double>> acceleration(valarray<valarray<double>> const& y_)
  {
	  valarray<valarray<double>> a = valarray<valarray<double>>(valarray<double>(0.0, 4),3);
	  int m_i;
	  int m_j;
	  for(int i = 0; i<3; ++i){
		  if(i==0){
			  m_i = 1;
			  m_j = 2;}
		 if (i == 1){
			 m_i = 0;
			  m_j = 2;}
		 if (i==2){
			 m_i = 1;
			  m_j = 2;}
		 for(int k =0; k<4; ++k){
			 if(k<2){ a[i][k] = y_[i][k+2];		 
		}else{a[i][k] = -G*M[m_i]*(y_[i][k-2]-y_[m_i][k-2])/(pow(sqrt((y_[i][0]-y_[m_i][0])*(y_[i][0]-y_[m_i][0]) + (y_[i][1]-y_[m_i][1])*(y_[i][1]-y_[m_i][1])),3))  - G*M[m_j]*(y_[i][k-2]-y_[m_j][k-2])/(pow(sqrt((y_[i][0]-y_[m_j][0])*(y_[i][0]-y_[m_j][0]) + (y_[i][1]-y_[m_j][1])*(y_[i][1]-y_[m_j][1])),3)); 
	
			}
			 
		 }
		  }
	return a;
  }
  
  
   valarray<valarray<double>> mult( valarray<valarray<double>> const& v, double a){
	   valarray<valarray<double>> r = v;
	   r[0] *= a;
	   r[1] *= a;
	   r[2] *= a;
	   return r;
	   }

    void step(){
		
	valarray<valarray<double>> k1 = mult(acceleration(y),dt);
	valarray<valarray<double>> k2= mult(acceleration(y + mult(k1,1.0/2.0)),2);
	valarray<valarray<double>> k3= mult(acceleration(y+mult(k2,1/2.0)),dt);
	valarray<valarray<double>> k4 = mult(acceleration(y+k3),dt);	
	y += mult((k1 + mult(k2,2) + mult(k3,2) +k4),(1.0/6.0));
		};
public :
    Exercice4(int argc, char* argv[])
    {
        string inputPath("configuration.in"); // Fichier d'input par defaut
        if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice4_2021_student.exe config_perso.in")
            inputPath = argv[1];

        ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

        for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3_2021_student.exe config_perso.in input_scan=[valeur]")
           configFile.process(argv[i]);

        tFin     = configFile.get<double>("tFin");
        nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire la nombre de pas de temps
        Omega    = configFile.get<double>("Omega");
        M[0]     =configFile.get<double>("M_A");
        M[1]     =configFile.get<double>("M_L");
        M[2]     =configFile.get<double>("M_T");
        y[0][0]     =configFile.get<double>("x_A");
        y[0][1]     =configFile.get<double>("y_A");
        y[0][2]     =configFile.get<double>("vx_A");
        y[0][3]     =configFile.get<double>("vy_A");
        y[1][0]     =configFile.get<double>("x_L");
        y[1][1]     =configFile.get<double>("y_L");
        y[1][2]     =configFile.get<double>("vx_L");
        y[1][3]     =configFile.get<double>("vy_L");
        y[2][0]     =configFile.get<double>("x_T");
        y[2][1]     =configFile.get<double>("y_T");
        y[2][2]     =configFile.get<double>("vx_T");
        y[2][3]     =configFile.get<double>("vy_T");
        sampling = configFile.get<int>("sampling");

        dt = tFin / nsteps; // calculer le time step

        // Ouverture du fichier de sortie
        outputFile = new ofstream(configFile.get<string>("output").c_str());
        outputFile->precision(15);
    };
    ~Exercice4()
    {
        outputFile->close();
        delete outputFile;
    };

    void run()
    {
        t = 0.;
        last = 0;
        printOut(true);
        while( t < tFin-0.5*dt )
        {
            step();
            t += dt;
            printOut(false);
        }
        printOut(true);
    };

};

int main(int argc, char* argv[])
{
    Exercice4 engine(argc, argv);
    engine.run();

    cout << "Fin de la simulation." << endl;

    return 0;
}
