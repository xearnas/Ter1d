#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"
#include "/mnt/c/users/simon/desktop/cours/enseirb/2a/c++/eigen-3.3.8/Eigen/Dense"
#include "/mnt/c/users/simon/desktop/cours/enseirb/2a/c++/eigen-3.3.8/Eigen/Sparse"

using namespace std;
using namespace Eigen;
int main(int argc, char const *argv[])
{
  if (argc < 2)
  {
      cout << "Please, enter the name of your data file." << endl;
      exit(0);
  }
  const string data_file_name = argv[1];

  // ----------------------- LECTURE DU FICHIER --------------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->Read_data_file();

  //------------------------- DECLARATION DES POINTEURS ET METHODES -------------------

  Mesh* msh = new Mesh(data_file);
  msh->Read_mesh();
  Function* fct = new Function(data_file);

  FiniteVolume* fin_vol = NULL;

  TimeScheme* time_scheme = NULL;

  if ( data_file->Get_numerical_flux_choice()=="Rusanov")
  {
    fin_vol = new Rusanov(fct,data_file,msh);
  }
  else if ( data_file->Get_numerical_flux_choice()=="LaxF")
  {
    fin_vol = new LaxF(fct,data_file,msh);
  }
  else
  {
    cout << "Please choose an implemented spatial scheme" << endl;
    exit(0);
  }

  if ( data_file->Get_scheme() == "Explicit" )
  {
    time_scheme = new Explicit(data_file,fct,msh,fin_vol);
  }
  else
  {
    cout << "Please enter an implemented time scheme" << endl;
    exit(0);
  }

  //----------------- RESOLUTION ------------------------------------------------

  cout << "Let's time loop" << endl;
  cout << "--------------------------------------" << endl;


  fin_vol->Save_sol(time_scheme->Get_sol(),time_scheme->Get_t());
  while (time_scheme->Get_t()<data_file->Get_tmax())
  {
    time_scheme->Advance();
    fin_vol->Save_sol(time_scheme->Get_sol(),time_scheme->Get_t());
    cout << "Iteration at t = " << time_scheme->Get_t() << " with a dt = " << fin_vol->Get_dt() <<  " ended correctly !" << endl;
  }

  cout << "Seems like it ends correctly" << endl;

  ///system(("atom " + data_file->Get_results() + "h_u_t=" + to_string(time_scheme->Get_t()) + ".txt").c_str());

  //// ---------- Ecrasement des pointeurs- --------------------------------------
  delete data_file;
  delete msh;
  delete fct;
  delete fin_vol;
  delete time_scheme;



  return 0;
}
