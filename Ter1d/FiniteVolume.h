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
  Mesh1D* _msh1d;
  Mesh2D* _msh2d;
  double _dt;
  MatrixXd _flux;
  //Matrice B avec les vap


public:
  FiniteVolume(Function* function, DataFile* data_file, Mesh1D* mesh1d,Mesh2D* mesh2d);
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
  virtual MatrixXd Get_B()=0;
  virtual void MajFlux(MatrixXd Vn) = 0;
};

class Rusanov : public FiniteVolume
{
private:
  MatrixXd _B;
  double _biundemi;
public:
  Rusanov(Function* function, DataFile* data_file,  Mesh1D* mesh1d,Mesh2D* mesh2d);
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
  LaxF(Function* function, DataFile* data_file, Mesh1D* mesh1d,Mesh2D* mesh2d);
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
