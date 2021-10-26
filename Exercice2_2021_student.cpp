#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* definir a fonction template pour calculer le produit interne
   entre deux valarray
   inputs:
     array1: (valarray<T>)(N) vecteur de taille N
     array2: (valarray<T>)(N) vecteur de taille N
   outputs:
     produitInterne: (T) produit entre les deux vecteurs
*/
template<typename T> T produitInterne(valarray<T> const& array1,\
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
} 

/* definir a fonction template pour calculer la norm2 d'un valarray
   inputs:
     array: (valarray<T>)(N) vecteur de taille N
   outputs:
     norm2: (T) norm2 du vecteur
*/
template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

/* definir a fonction template pour calculer le produit vecteur
   entre 2 valarray de dimension 3
   inputs:
     array1, array2: (valarray<T>)(N) vecteurs de taille N
   outputs:
     produitVecteur: (T) produit vectoriel array1 x aray2 
*/
template<typename T> T produitVecteur(valarray<T> const& array1,\
valarray<T> const& array2){
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
} 


/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{

private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  double tfin=0.e0;       // Temps final
  unsigned int nsteps=1;  // Nombre de pas de temps
  double omega0=0.e0; 	  // vitesse initiale de rotation du balon en direction y
  double g=0.e0;	  // acceleration gravitationelle
  double mass=1.e0,R=0.e0;        // Masse et rayon du ballon
  double mu=0.e0,rho=0.e0; // facteur de forme et densite de l'air     
  valarray<double> x0=valarray<double>(0.e0,2); // vecteur contenant la position initiale du ballon en
			 		        // utilisant la sequence index-valeur: 0-x, 1-z
  valarray<double> v0=valarray<double>(0.e0,2); // vecteur contenant la vitesse initiale du ballon en
			  		        // utilisant la sequence index-valeur: 0-vx, 1-vz
  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie


  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
       write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      // TODO calculer l'energie mecanique      
      double energyMecanique=0.0;
      energyMecanique = 0.5*mass*v[1]*v[1]+ 0.5*mass*v[0]*v[0]+ mass*g*x[1];
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << v[0] << " " << v[1] << " " \
      << energyMecanique << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;

protected:

  // donnes internes
  double t,dt;  // Temps courant pas de temps
  double Omega; // vitesse angulaire de la trajectoire (rad/s)
  valarray<double> x=valarray<double>(2); // Position actuelle du ballon 
  valarray<double> v=valarray<double>(2); // Vitesse actuelle du ballon

  
  // TODO
  /* Calcul de l'acceleration
     outputs:
       a: (valarray<double>)(2) acceleration selon (x,z)
  */
  void acceleration(valarray<double>& a) const
  {
    // calcul de l'acceleration
    a[0] =Omega * v[1];
    a[1] =-g-Omega*v[0];
  }

public:

  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */
  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin",tfin);	      // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire la nombre de pas de temps
    mass     = configFile.get<double>("mass",mass);	      // lire la mass du ballon
    R        = configFile.get<double>("R",R);	              // lire le rayon du ballon
    g        = configFile.get<double>("g",g);		      // lire acceleration gravite
    mu       = configFile.get<double>("mu",mu);		      // lire facteur de forme
    rho      = configFile.get<double>("rho",rho);	      // lire densite de l'air
    omega0   = configFile.get<double>("omega0",omega0);	      // lire composante y vitesse angulaire
    x0[0]    = configFile.get<double>("x0",x0[0]);	      // lire composante x position initiale
    x0[1]    = configFile.get<double>("z0",x0[1]);	      // lire composante z position initiale
    v0[0]    = configFile.get<double>("vx0",v0[0]);	      // lire composante x vitesse initiale
    v0[1]    = configFile.get<double>("vz0",v0[1]);	      // lire composante z vitesse initiale
    
    sampling = configFile.get<unsigned int>("sampling",sampling); // lire le parametre de sampling

    dt = tfin / nsteps; // calculer le time step
    // TODO calculer la vitesse angulaire de la trajectoire
    Omega = (mu*rho*R*R*R*omega0)/mass;

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps
    x = x0;   // initialiser la position
    v = v0;   // initialiser la vitesse
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps
    for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
    {
      step();  // faire un pas
      printOut(false); // ecrire le pas de temps actuel
    }
    printOut(true); // ecrire le dernier pas de temps
  };

};

// Extension de la class Engine implementant l'integrateur d'Euler explicite
class EngineEulerExplicite: public Engine
{
public:
  
  // construire la class Engine
  EngineEulerExplicite(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler explicite
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,v.size()); // definir l'acceleration
      acceleration(a);// mis a jour de l'acceleration

    x += v*dt; // mise a jour de la position
    v += a*dt; // mise a jour de la vitesse
    t += dt; // mise a jour du temps
  }
};


// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineRungeKutta2: public Engine
{
public:

  // construire la class Engine
  EngineRungeKutta2(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Runge-kutta 2
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,v.size()); // definir l'acceleration
    // sauvgarder la position et la vitesse
    acceleration(a);
    valarray<double> x_=x ; // position
    valarray<double> v_=v; // vitesse
   valarray<double> kx1= dt * v;
    valarray<double> kv1 =dt* a;
    v += 0.5*kv1;
    x += kx1*0.5;
    acceleration(a);
    valarray<double> kx2 = dt*v;
    valarray<double> kv2 = dt*a;
    x = x_ + kx2;  // mise a jour de la position
    v = v_ + kv2; // mise a jour de la vitesse
    t += dt*0.5; // mise a jour du temps

  }
};


// Extension de la class Engine implementant l'integrateur d'Euler implicite
class EngineEulerImplicite: public Engine
{
private:
  unsigned int maxit=1000; // nombre maximal d'iterations
  double tol=1.e12;        // tolerance methode iterative
public:
  
  // construire la class Engine
  EngineEulerImplicite(ConfigFile configFile): Engine(configFile){

    tol = configFile.get<double>("tol"); // lire la tolerance pour la methode iterative
    maxit = configFile.get<unsigned int>("maxit"); // lire le nombre d'iterations maximal 

  }

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler Implicite, avec la methode du point fixe
     pour resoude le probleme implicite
  */
  void step()
  {
    unsigned int iteration=0;
    double error=999e0;
    valarray<double> a=valarray<double>(0.e0,2);
   
    valarray<double> xi=valarray<double>(x);
    valarray<double> vi=valarray<double>(v);
    valarray<double> xold=valarray<double>(x);
    valarray<double> vold=valarray<double>(v);
  valarray<double> dv=valarray<double>(0.0, error);
  valarray<double> dx=valarray<double>(0.0, error);
    t += dt; // mise a jour du temps
for (int i(0); (i < maxit) and (error > tol); ++i){
   vold= v;
   xold = x;
   acceleration(a);
    x = xi + dt * vold; // mise a jour de la position
    v = vi + dt * a ; // mise a jour de la vitesse
    dv = v - vold;
    dx = x- xold;
    error = sqrt(norm2(dx)*norm2(dx) + norm2(dx));
    ++iteration;
    if(iteration>=maxit){
      cout << "WARNING: maximum number of iterations reached, error: " << error << endl;
    }}
    
  }
};




// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("EulerExplicite"/"EX", "RungeKutta2"/"RK2" ou "EulerImplicite"/"EI")
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser
  if(schema == "EulerExplicite" || schema == "EX")
  {
    // initialiser une simulation avec schema Euler
    engine = new EngineEulerExplicite(configFile);
  }
  else if(schema == "RungeKutta2" || schema == "RK2")
  {
    // initialiser une simulation avec schema runge-kutta 2
    engine = new EngineRungeKutta2(configFile);
  }
  else if(schema == "EulerImplicite" || schema == "EI")
  {
    // initialiser une simulation avec schema Euler Implicite
    engine = new EngineEulerImplicite(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
