#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

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


  ///-----------------MESH-------------------------
  cout << data_file->Get_mesh_name() << endl;
  Mesh1D* msh1d = new Mesh1D(data_file);
  Mesh2D* msh2d = new Mesh2D(data_file);
  if ( data_file->Get_mesh_name()== "Mesh1D")
  {
    msh1d->Initialize1D();
  }
  else if ( data_file->Get_mesh_name()== "Mesh2D")
  {
    msh2d->Initialize2D();
  }
  else
  {
    cout << "Please choose an implemented mesh" << endl;
    exit(0);
  }
  ///------------------------FUNCTION------------------------------
  Function* fct = new Function(data_file);

  ///------------------------FINITE VOLUME------------------------------
  FiniteVolume* fin_vol = NULL;

  if ( data_file->Get_numerical_flux_choice()=="Rusanov")
  {
    fin_vol = new Rusanov(fct,data_file,msh1d,msh2d);
  }
  else if ( data_file->Get_numerical_flux_choice()=="LaxF")
  {
    fin_vol = new LaxF(fct,data_file,msh1d,msh2d);
  }
  else
  {
    cout << "Please choose an implemented spatial scheme" << endl;
    exit(0);
  }

  ///------------------------TIME SCHEME------------------------------

  TimeScheme* time_scheme = NULL;

  if ( data_file->Get_scheme() == "Explicit" )
  {
    time_scheme = new Explicit(data_file,fct,msh1d,msh2d,fin_vol);
  }
  else
  {
    cout << "Please enter an implemented time scheme" << endl;
    exit(0);
  }

  ///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ///-----------------------END OF DECLARATION-----------------------------------------
  ///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //----------------- RESOLUTION ------------------------------------------------

  cout << "Let's time loop" << endl;
  cout << "--------------------------------------" << endl;

  int compteur=1;
  fin_vol->Save_sol(time_scheme->Get_sol(),time_scheme->Get_t());
  while (time_scheme->Get_t()<data_file->Get_tmax())
  {

    time_scheme->Advance();
    fin_vol->Save_sol(time_scheme->Get_sol(),time_scheme->Get_t());
    cout << "Iteration at t = " << time_scheme->Get_t() << " with a dt = " << fin_vol->Get_dt() <<  " ended correctly !" << endl;
    compteur +=1;
    cout << "compteur = " << compteur << endl;
  }

  cout << "Seems like it ends correctly with compteur = " << compteur << endl;

  ///system(("atom " + data_file->Get_results() + "h_u_t=" + to_string(time_scheme->Get_t()) + ".txt").c_str());

  //// ---------- Ecrasement des pointeurs- --------------------------------------
  delete data_file;
  delete msh1d;
  delete msh2d;
  delete fct;
  delete fin_vol;
  delete time_scheme;



  return 0;
}
