#include <string>
#include <vector>
#include <iostream>
#include "/mnt/c/users/simon/desktop/cours/enseirb/2a/c++/eigen-3.3.8/Eigen/Dense"
#include "/mnt/c/users/simon/desktop/cours/enseirb/2a/c++/eigen-3.3.8/Eigen/Sparse"
using namespace std;
using namespace Eigen;
#pragma once

class DataFile
{
private:
  string _file_name;
  //Temps
  double _tmin, _tmax, _dt;
  //Espace
  double _xmin, _xmax, _dx;
  //Pensanteur
  double _g;

  string _mesh_name;
  string _scheme;
  string _numerical_flux_choice;
  string _results;
  string _which_scenario;


  bool _print_info;
  bool _if_mesh_name;
  bool _if_tmin;
  bool _if_tmax;
  bool _if_xmin;
  bool _if_xmax;
  bool _if_dx;
  bool _if_scheme;
  bool _if_numerical_flux_choice;
  bool _if_results;
  bool _if_g;


public:
  DataFile(string file_name);
  //Lire le fichier
  void Read_data_file();
  //Car non linéaire
  void Adapt_dt(double dt)
  {
    _dt=dt;
  };
  //Les get
  const double Get_tmin() const
  {
    return _tmin;
  };
  const double Get_tmax() const
  {
    return _tmax;
  };
  const double Get_dt() const
  {
    return _dt;
  };
  const double Get_xmin() const
  {
    return _xmin;
  };
  const double Get_xmax() const
  {
    return _xmax;
  };
  const double Get_dx() const
  {
    return _dx;
  };
  const double Get_g() const
  {
    return _g;
  };
  const int Get_Nx() const
  {
    return ((_xmax - _xmin) / _dx);
  };
  const string Get_scenario() const
  {
    return _which_scenario;
  };
  const string Get_mesh_name() const
  {
    return _mesh_name;
  };
  const string Get_scheme() const
  {
    return _scheme;
  };
  const string Get_numerical_flux_choice() const
  {
    return _numerical_flux_choice;
  };
  const string Get_results() const
  {
    return _results;
  };
};
