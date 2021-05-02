#include "Function.h"
#include "DataFile.h"
#include "Mesh.h"
#include <vector>

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
  MatrixXd F(VectorXd V);
  MatrixXd Ftild(VectorXd Vi, VectorXd Viplus1,double b, int i);
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
  vector <MatrixXd> _All_Matrices;
  vector <int> _All_Compteurs;
  
public:
  LaxF(Function* function, DataFile* data_file, Mesh1D* mesh1d,Mesh2D* mesh2d);
  void MajFlux(MatrixXd Vn);
  void MajB(MatrixXd Vn);
  MatrixXd F(VectorXd V);
  MatrixXd Ftild(VectorXd Vi, VectorXd Viplus1,double b,int i);
  void InitializeMatrices();
  double Get_bi()
  {
    return _biundemi;
  };

  MatrixXd Get_B()
  {
    return _B;
  };
};
