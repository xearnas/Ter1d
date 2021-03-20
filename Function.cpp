#include "Function.h"
#include <cmath>

using namespace std;
using namespace Eigen;



Function::Function(DataFile* data_file) :
    _df(data_file), _g(_df->Get_g())
{
    cout << "Functions are defined" << endl;
}

double Function::Initial_condition_forH(double x)
{
    if (_df->Get_scenario() == "Barreer")
    {
      if (x<0)
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

double Function::Initial_condition_forU(double x)
{
  if (_df->Get_scenario() == "Barreer")
  {
    return 0.;
  }
  else
  {
    cout << "Choose an implanted scenario." << endl;
    exit(0);
  }
}

double Function::Source_term(double x, double t)
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

double Function::maxdeB(MatrixXd B)
{
    double max=0;
    for (int i=0; i< B.rows() ; i++)
    {
      if (B(i,0) > max)
        max = B(i,0);
      if (B(i,1) > max)
        max = B(i,1);
    }
    return max;
}

VectorXd Function::phi(VectorXd theta)
{
  VectorXd res(2);
  if ((theta(0)>=0) && (theta(1)>=0))
  {
    res(0)=min(1.0,(double)theta(0));
    res(1)=min(1.0,(double)theta(1));
  }
  else
  {
    res(0)=0.0;
    res(1)=0.0;
  }
  //res(0)=(theta(0)+abs(theta(0)))/(1+theta(0));
  //res(1)=(theta(1)+abs(theta(1)))/(1+theta(1));
  return res;
}
