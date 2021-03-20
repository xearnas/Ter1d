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
  Mesh* _msh;
  FiniteVolume* _fv;

  MatrixXd _Vn;
  double _t;

public:
  TimeScheme(DataFile* data_file, Function* function, Mesh* mesh, FiniteVolume* finite_volume);
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
  Explicit(DataFile* data_file, Function* function, Mesh* mesh, FiniteVolume* finite_volume);
  //Une etape en temps
  void Advance();
};
