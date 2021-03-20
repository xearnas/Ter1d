#include <fstream>
#include <iostream>
#include <cmath>
#include <regex>
#include "DataFile.h"
using namespace std;
using namespace Eigen;

DataFile::DataFile(string file_name)
  : _file_name(file_name),  _if_mesh_name(false), _if_tmin(false), _if_tmax(false) , _if_xmin(false), _if_xmax(false), _if_dx(false) , _if_scheme(false), _if_numerical_flux_choice(false), _if_results(false), _if_g(false)
{
}

void DataFile::Read_data_file()
{
  ifstream data_file(_file_name.data());
  if (!data_file.is_open())
  {
      cout << "Unable to open file " << _file_name << endl;
      exit(0);
  }
  else
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
  }

  string file_line;

  while (!data_file.eof())
  {
      getline(data_file, file_line);

      if (file_line.find("#") !=string::npos)
      {
          // Ignore this line (comment)
      }
      else
      {
          if (file_line.find("mesh") != string::npos)
          {
              data_file >> _mesh_name;
              _if_mesh_name = true;
          }

          if (file_line.find("numerical_flux") != string::npos)
          {
              data_file >> _numerical_flux_choice;
              _if_numerical_flux_choice = true;
              if ((_numerical_flux_choice != "Rusanov") && (_numerical_flux_choice != "LaxF"))
              {
                  cout << "Only Rusanov and Lax Frederichs numerical flux is implemented." << endl;
                  exit(0);
              }
          }

          if (file_line.find("tmin") != string::npos)
          {
              data_file >> _tmin;
              _if_tmin = true;
          }

          if (file_line.find("tmax") != string::npos)
          {
              data_file >> _tmax;
              _if_tmax = true;
          }

          if (file_line.find("xmin") != string::npos)
          {
              data_file >> _xmin;
              _if_xmin = true;
          }

          if (file_line.find("xmax") != string::npos)
          {
              data_file >> _xmax;
              _if_xmax = true;
          }

          if (file_line.find("dx") != string::npos)
          {
              data_file >> _dx;
              _if_dx = true;
          }

          if (file_line.find("which_scenario") != string::npos)
          {
              data_file >> _which_scenario;
          }

          if (file_line.find("scheme") != string::npos)
          {
              data_file >> _scheme;
              _if_scheme = true;
              if ((_scheme != "Explicit"))
              {
                  cout << "Only Explicit scheme is implemented." << endl;
                  exit(0);
              }
          }

          if (file_line.find("g") != string::npos)
          {
              data_file >> _g;
              _if_g = true;
          }

          if (file_line.find("results") != string::npos)
          {
              data_file >> _results;
              _if_results = true;
          }
      }
  }
  if (!_if_tmin)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.) is used for tmin." << endl;
      _tmin = 0.;
  }
  if (!_if_tmax)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.1) is used for tmax." << endl;
      _tmax = 0.1;
  }
  if (!_if_scheme)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default scheme (Explicit scheme) is used." << endl;
      _scheme = "Explicit";
  }
  if (!_if_g)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (0.1) is used for g." << endl;
      _g = 9.81;
  }
  if (!_if_results)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default results folder name (results) is used." << endl;
      _results = "results";
  }
  cout << "-------------------------------------------------" << endl;
  if (!_if_mesh_name)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Do not forget to give the mesh name in the data file." << endl;
      exit(0);
  }
  if (!_if_numerical_flux_choice)
  {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (Rusanov) is used for the numerical flow." << endl;
      _numerical_flux_choice = "Rusanov";
  }
  // Conditions
  if (_scheme == "Explicit")
  {
      cout << "Beware to the CFL condition fixed at every time iteration." << endl;
  }

  if ((_which_scenario == "Barreer"))
  {
      cout << "-------------------------------------------------" << endl;
      cout << "The test case: " << _which_scenario << " has been chosen." << endl;
      cout << "The initial values of h will be 2 for x negative and 1 for other" << endl;
      cout << "-------------------------------------------------" << endl;
  }
  else
  {
      cout << "-------------------------------------------------" << endl;
      cout << "A scenario has to be chosen between the : Barrer , " << endl;
      cout << "-------------------------------------------------" << endl;
      exit(0);
  }

  // Supprimer les anciens résultats
  system(("rm -f ./" +_results + "*.txt").c_str());
  cout << "Last results erased" << endl;
  cout << "-------------------------------------" << endl;

}
