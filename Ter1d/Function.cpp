#include "Function.h"
#include <cmath>
#include <math.h>

using namespace std;
using namespace Eigen;



Function::Function(DataFile* data_file) :
    _df(data_file), _g(_df->Get_g())
{
}

double Function::Initial_condition_forH(double x, double y)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    if (_df->Get_scenario() == "Barreer")
    {

      if (x<(_df->Get_xmax()+_df->Get_xmin())/2.0)
      {
        return 2.;
      }
      else
      {
        return 1.;
      }
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    if (_df->Get_scenario() == "Barreer")
    {
      if (x<(_df->Get_xmax()+_df->Get_xmin())/2.0)
      {
        return 2.;
      }
      else
      {
        return 1.;
      }
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else
  {
    cout << "Choose an implanted mesh, we can't initialize H" << endl;
    exit(0);
  }

}

VectorXd Function::Initial_condition_forU(double x, double y)
{
  VectorXd speedvector;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    speedvector.resize(1);
    if (_df->Get_scenario() == "Barreer")
    {
      speedvector(0)=0.;
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    speedvector.resize(2);
    if (_df->Get_scenario() == "Barreer")
    {
      speedvector(0)=0.;
      speedvector(1)=0.;
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else
  {
    cout << "Choose an implanted mesh, we can't initialize U" << endl;
    exit(0);
  }
  return speedvector;
}

double Function::Source_term(double x, double y,double t)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    if (_df->Get_scenario() == "Barreer")
    {
      return 0;
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    if (_df->Get_scenario() == "Barreer")
    {
      return 0;
    }
    else
    {
      cout << "Choose an implanted scenario." << endl;
      exit(0);
    }
  }
  else
  {
    cout << "Choose an implanted mesh, we can't calculate S" << endl;
    exit(0);
  }
}

double Function::maxdeB(MatrixXd B)
{
    double max=0;
    for (int i=0; i< B.rows() ; i++)
    {
      for (int j=0;j<B.cols();++j)
      {
        if (B(i,j) > max)
          max = B(i,j);
      }
    }
    return max;
}

VectorXd Function::phi(VectorXd theta)
{
  VectorXd res(2);
  res(0)=(theta(0)+abs(theta(0)))/(1+theta(0));
  res(1)=(theta(1)+abs(theta(1)))/(1+theta(1));
  return res;
}

double Function::Produit_Scalaire(VectorXd U,VectorXd V)
{
  double res=0;
  if (U.rows()!=V.rows())
  {
    cout << "U et V ne sont pas de la même taille abruti !!" << endl;
    exit(0);
  }
  else
  {
    for (int i =0; i<U.rows();++i)
    {
      res+=U(i)*V(i);
    }
  }
  return res;
}
