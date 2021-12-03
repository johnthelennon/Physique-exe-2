#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include "ConfigFile.h"

using namespace std;
class Exercice4 {
private:
    double t, dt, tFin, eps;
    double m_L, m_T, m_A;
    double G=6.67430* pow(10,-11);
    double d=160;
    double R_T=6378.1;
    double D=384748; //distance terre lune
    bool adaptatif;
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
            double emec = 0.0;
            *outputFile << t << " ";
            for (unsigned int i(0); i < 4; ++i) {
                for (unsigned int j(0); j < 3; ++j) {
                    cout << y[i][j] << " ";
                }
                cout << endl;

                last = 1;
            }
        } else {
            last++;
        }
    }

    valarray<valarray<double>> acceleration(const double theta_, const double thetadot_, const double t_) ;

    void step(valarray<valarray<double>> const& y, double const& dt);

    double norm2(valarray<valarray<double>>y) ;
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
        eps      = configFile.get<double>("eps");
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
        adaptatif   =configFile.get<double>("adaptatif");
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


        void run() //0 pour l'asteroire, 1 pour la lune et 2 pour la terre
        {
            t = 0;
            unsigned int j = 0;
            last = 0; // initialise le parametre d'ecritur

            printOut(true); // ecrire premier pas de temps

            if(adaptatif) {
                valarray<double> x1(0.e0, 12);
                valarray<double> x2(0.e0, 12);
                int n = 4; //ordre du schéma
                double fact = 0.99;
                double d = 0.0;
                while(t<tFin) {

                    j+=1;
                    dt = min(tFin-t, dt); //Pour finir pile au bon temps
                     step(y, dt);
                     valarray<valarray<double >>y_temp=y;
                     step(y, dt/2.);
                     step(y, dt/2.);
                    d = norm2(y-y_temp); //ecrire fonction norme
                    if (d < eps) { //Vérification de la précision : assez précis
                        t += dt;
                        dt = dt*pow(eps/d, 1./(n+1));
                    } else { // Pas assez précis : itérer jusqu'à ce que ça le soit
                        while(d>eps) {
                            dt = fact*dt*pow(eps/d, 1./(n+1)); // Formule du cours
                            dt = min(tFin-t, dt);
                            step(y, dt);
                            valarray<valarray<double >>y_temp=y;
                            step(y, dt/2.);
                            step(y, dt/2.);
                            d = norm2(y-y_temp);
                        }
                        t+= dt;
                    }

                    printOut(false); // ecrire pas de temps actuel

                }
            } else {	// Dans le cas pas de temps fixe.
                dt = tFin/nsteps;
                while(t<tFin) {
                    dt = min(tFin-t, dt);
                    t+=dt;
                    j+=1;
                    step(y, dt);
                    printOut(false);
                }
            }

            printOut(true); // ecrire dernier pas de temps
        }

    };



valarray<valarray<double>> Exercice4::acceleration(const double theta_, const double thetadot_, const double t_) {
    return valarray<valarray<double>>();
}

void Exercice4::step(const valarray<valarray<double>> &y, const double &dt) {

}

double Exercice4::norm2(valarray<valarray<double>> y) {
    return 0;
}

int main(int argc, char* argv[])
{
    Exercice4 engine(argc, argv);
    engine.run();

    cout << "Fin de la simulation." << endl;

    return 0;
}
