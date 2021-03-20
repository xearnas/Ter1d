#include "Mesh.h"

using namespace std;
using namespace Eigen;

Mesh::Mesh(DataFile* data_file)
: _df(data_file)
{
  cout << "Mesh is entered" << endl;
}

void Mesh::Read_mesh()
{
  if (_df->Get_mesh_name() == "Mesh1D")
  {
    this->Initialize1D();
  }
  else
  {
    cout << "Mesh inconnue au bataillon" << endl;
    exit(0);
  }
}

void Mesh::Initialize1D()
{
  _noeuds.resize(_df->Get_Nx()+1);
  _centres.resize(_df->Get_Nx());
  for (int i=0; i<_df->Get_Nx() ; i++)
  {
    _noeuds(i)=_df->Get_xmin() + i* _df->Get_dx();
    _centres(i)= _noeuds(i) + _df->Get_dx()/2.;
  }
  _noeuds(_df->Get_Nx()) = _df->Get_xmax();
}
