#include "Function.h"
#include "DataFile.h"
#include "Mesh.h"

using namespace std;
using namespace Eigen;
#pragma once
class FiniteVolume
{
protected:
  Function* _fct;
  DataFile* _df;
  Mesh* _msh;
  double _dt;
  MatrixXd _flux;
  //Matrice B avec les vap


public:
  FiniteVolume(Function* function, DataFile* data_file, Mesh* mesh);
  virtual ~FiniteVolume();

  MatrixXd Initial_condition();
  VectorXd Source_term(double t);
  void Save_sol(MatrixXd Vn, double t);
  MatrixXd Get_flux()
  {
    return _flux;
  };
  double Get_dt()
  {
    return _dt;
  };

  virtual void MajFlux(MatrixXd Vn) = 0;
  //ATTENTION A LES BOUGER AVANT DE COMPILER
};

class Rusanov : public FiniteVolume
{
private:
  MatrixXd _B;
  double _biundemi;
public:
  Rusanov(Function* function, DataFile* data_file, Mesh* mesh);
  void MajFlux(MatrixXd Vn);
  void MajB(MatrixXd Vn);
  VectorXd F(VectorXd V);
  MatrixXd Ftild(VectorXd Vi, VectorXd Viplus1,double b);
  double Get_bi()
  {
    return _biundemi;
  };

  MatrixXd Get_B()
  {
    return _B;
  };
};

class LaxF : public FiniteVolume
{
private:
  MatrixXd _B;
  double _biundemi;
public:
  LaxF(Function* function, DataFile* data_file, Mesh* mesh);
  void MajFlux(MatrixXd Vn);
  void MajB(MatrixXd Vn);
  VectorXd F(VectorXd V);
  MatrixXd Ftild(VectorXd Vi, VectorXd Viplus1,double b);
  double Get_bi()
  {
    return _biundemi;
  };

  MatrixXd Get_B()
  {
    return _B;
  };
};
