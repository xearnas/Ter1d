#include "Function.h"
#include "DataFile.h"
#include "Mesh.h"
#include "FiniteVolume.h"

using namespace std;
using namespace Eigen;
#pragma once
class TimeScheme
{
protected:
  Function* _fct;
  DataFile* _df;
  Mesh1D* _msh1d;
  Mesh2D* _msh2d;
  FiniteVolume* _fv;

  MatrixXd _Vn;
  double _t;

public:
  TimeScheme(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume);
  virtual ~TimeScheme();

  virtual void Advance() = 0;
  // Permet de récupérer _sol
  MatrixXd Get_sol()
  {
    return _Vn;
  };
  double Get_t()
  {
    return _t;
  };
};

class Explicit : public TimeScheme
{
public:
  //Constructeur
  Explicit(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume);
  //Une etape en temps
  void Advance();
};

class RK4 : public TimeScheme
{
public:
  //Constructeur
  RK4(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume);
  // Etape en temps
  void Advance();
};
