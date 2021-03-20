#include "Function.h"
#include "DataFile.h"

using namespace std;
using namespace Eigen;
#pragma once
class Mesh
{
private:
  DataFile* _df;
  VectorXd _noeuds;
  VectorXd _centres;

public:
  Mesh(DataFile* data_file);
  void Read_mesh();
  void Initialize1D();
  VectorXd Get_noeuds()
  {
    return _noeuds;
  };
  VectorXd Get_centres()
  {
    return _centres;
  };

};
