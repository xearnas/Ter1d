#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

TimeScheme::TimeScheme(DataFile* data_file, Function* function, Mesh* mesh, FiniteVolume* finite_volume)
: _df(data_file), _fct(function), _msh(mesh), _fv(finite_volume), _t(_df->Get_tmin()),_Vn(_fv->Initial_condition())
{
    _Vn.resize(_df->Get_Nx()+1,2);
}

Explicit::Explicit(DataFile* data_file, Function* function, Mesh* mesh, FiniteVolume* finite_volume)
: TimeScheme(data_file,function,mesh,finite_volume)
{
}

TimeScheme::~TimeScheme()
{
}

void Explicit::Advance()
{
  VectorXd TermSource;
  TermSource = _fv->Source_term(_t);

  _fv->MajFlux(_Vn);

  for (int i=1;i<_df->Get_Nx();i++)
  {
    _Vn.row(i) += -_fv->Get_dt()/_df->Get_dx() * (_fv->Get_flux().row(i) - _fv->Get_flux().row(i-1));
  }
  _t += _fv->Get_dt();
}
