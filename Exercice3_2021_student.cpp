#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.h"

using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double m, g, L;
  double r, Omega, kappa;
  double theta, thetadot;
  int sampling;
  unsigned int nsteps;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = L*L*thetadot*thetadot*m*0.5 - m*g*L*cos(theta); // TODO: Evaluer l'energie mecanique
      double pnc =  -kappa*L*L*thetadot*thetadot +m*L*thetadot*Omega*Omega*cos(Omega*t-theta) ; // TODO: Evaluer la puissance des forces non conservatives

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  double acceleration(const double theta_, const double thetadot_, const double t_)
  {
	double a= (1.0/L)*(-g*sin(theta_)-(kappa/m)*L*thetadot_ + r*Omega*Omega*cos(Omega*t-theta_));  // TODO: Evaluer l'accélération
	return a;
  }

  void step()
  {double t_ini = t ;
double theta_ini = theta;
double a = acceleration(theta,thetadot,t);
    theta = theta + thetadot*dt + 0.5*a*dt*dt;     // TODO: Implémentation du schéma de Verlet
    double thetadot_1_2 = thetadot + 0.5*a*dt;
    double a1= acceleration(theta_ini, 0, t_ini);
    double a2=acceleration(theta,0,t+dt);
    double a3=acceleration(0,thetadot_1_2,(1.0/Omega)*M_PI/2.0);
    thetadot = thetadot + 0.5*dt*(a1 + a2) + a3*dt;
    t+=dt;
  }


public:

  Exercice3(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3_2021_student.exe config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3_2021_student.exe config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin     = configFile.get<double>("tFin");
    nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire la nombre de pas de temps
    r        = configFile.get<double>("r");
    Omega    = configFile.get<double>("Omega");
    kappa    = configFile.get<double>("kappa");
    m        = configFile.get<double>("m");
    g        = configFile.get<double>("g");
    L        = configFile.get<double>("L");
    theta    = configFile.get<double>("theta0");
    thetadot = configFile.get<double>("thetadot0");
    sampling = configFile.get<int>("sampling");
    
    dt = tFin / nsteps; // calculer le time step
    
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  };

  ~Exercice3()
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
  Exercice3 engine(argc, argv);
  engine.run();

  cout << "Fin de la simulation." << endl;

  return 0;
}
