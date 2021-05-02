#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

TimeScheme::TimeScheme(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume)
: _df(data_file), _fct(function), _msh1d(mesh1d), _msh2d(mesh2d), _fv(finite_volume), _t(_df->Get_tmin()),_Vn(_fv->Initial_condition())
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    _Vn.resize(_df->Get_Nx()+1,2);
  }
  else if (_df->Get_mesh_name() == "Mesh2D")
  {
    _Vn.resize(_msh2d->Get_triangles().size(),3);
  }
  else
  {
    cout << "Choose an implanted mesh, we can't size correctly the solution vector" << endl;
    exit(0);
  }
}

Explicit::Explicit(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume)
: TimeScheme(data_file,function,mesh1d,mesh2d,finite_volume)
{
}

RK4::RK4(DataFile* data_file, Function* function, Mesh1D* mesh1d, Mesh2D* mesh2d, FiniteVolume* finite_volume)
: TimeScheme(data_file,function,mesh1d,mesh2d,finite_volume)
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
  if (_df->Get_mesh_name() == "Mesh1D")
  {
    for (int i=0;i<_df->Get_Nx()+1;i++)
    {
      _Vn.row(i) += -_fv->Get_dt()/_df->Get_dx() * (_fv->Get_flux().row(i+1) - _fv->Get_flux().row(i));
    }
  }
  else if (_df->Get_mesh_name() == "Mesh2D")
  {
    for (int i=0;i<_msh2d->Get_edges().size();i++)
    {
      int t1 = _msh2d->Get_edges()[i].Get_T1();
      int t2 = _msh2d->Get_edges()[i].Get_T2();
      if (t2==-1)
      {
        _Vn.row(t1) -= _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t1] * _fv->Get_dt() * _fv->Get_flux().row(i);
      }
      else
      {
        _Vn.row(t1) -= _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t1] * _fv->Get_dt() * _fv->Get_flux().row(i);
        _Vn.row(t2) += _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t2] * _fv->Get_dt() * _fv->Get_flux().row(i);
      }

    }

  }
  else
  {
    cout << "Choose an implanted mesh, we can't calculate the solution at the next time step" << endl;
    exit(0);
  }
  _t += _fv->Get_dt();
}

void RK4::Advance()
{
  VectorXd TermSource;
  TermSource = _fv->Source_term(_t);

  _fv->MajFlux(_Vn);
  if (_df->Get_mesh_name() == "Mesh1D")
  {
    MatrixXd k1(1,2),k2(1,2),k3(1,2),k4(1,2);
    for (int i=0;i<_df->Get_Nx()+1;i++)
    {
      k1 = -(_fv->Get_flux().row(i+1) - _fv->Get_flux().row(i));
      k2 = -(_fv->Ftild(_fv->Get_dt()/2 * k1.row(0) + Vn.row(i),_fv->Get_dt()/2 * k1.row(0) + Vn.row(i+1), max(max(_fv->Get_B()(i+1,0),_fv->Get_B()(i+1,1)),max(_fv->Get_B()(i,0),_fv->Get_B()(i,1))) ,i+1) - _fv->Ftild(_fv->Get_dt()/2 * k1.row(0) + Vn.row(i-1),_fv->Get_dt()/2 * k1.row(0) + Vn.row(i), max(max(_fv->Get_B()(i,0),_fv->Get_B()(i,1)),max(_fv->Get_B()(i-1,0),_fv->Get_B()(i-1,1))) ,i));
      k3 = -(_fv->Ftild(_fv->Get_dt()/2 * k2.row(0) + Vn.row(i),_fv->Get_dt()/2 * k2.row(0) + Vn.row(i+1), max(max(_fv->Get_B()(i+1,0),_fv->Get_B()(i+1,1)),max(_fv->Get_B()(i,0),_fv->Get_B()(i,1))) ,i+1) - _fv->Ftild(_fv->Get_dt()/2 * k2.row(0) + Vn.row(i-1),_fv->Get_dt()/2 * k2.row(0) + Vn.row(i), max(max(_fv->Get_B()(i,0),_fv->Get_B()(i,1)),max(_fv->Get_B()(i-1,0),_fv->Get_B()(i-1,1))) ,i));
      k4 = -(_fv->Ftild(_fv->Get_dt() * k3.row(0) + Vn.row(i),_fv->Get_dt() * k3.row(0) + Vn.row(i+1), max(max(_fv->Get_B()(i+1,0),_fv->Get_B()(i+1,1)),max(_fv->Get_B()(i,0),_fv->Get_B()(i,1))) ,i+1) - _fv->Ftild(_fv->Get_dt() * k3.row(0) + Vn.row(i-1),_fv->Get_dt() * k3.row(0) + Vn.row(i), max(max(_fv->Get_B()(i,0),_fv->Get_B()(i,1)),max(_fv->Get_B()(i-1,0),_fv->Get_B()(i-1,1))) ,i));
      _Vn.row(i) += -_fv->Get_dt()/_df->Get_dx() * (_fv->Get_flux().row(i+1) - _fv->Get_flux().row(i));
      //_Vn.row(i) += _fv->Get_dt()/(6*_df->Get_dx()) * (k1.row(0)+2*k2.row(0)+2*k3.row(0)+k4.row(0));
    }
  }
  else if (_df->Get_mesh_name() == "Mesh2D")
  {
    for (int i=0;i<_msh2d->Get_edges().size();i++)
    {
      int t1 = _msh2d->Get_edges()[i].Get_T1();
      int t2 = _msh2d->Get_edges()[i].Get_T2();
      if (t2==-1)
      {
        _Vn.row(t1) -= _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t1] * _fv->Get_dt() * _fv->Get_flux().row(i);
      }
      else
      {
        _Vn.row(t1) -= _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t1] * _fv->Get_dt() * _fv->Get_flux().row(i);
        _Vn.row(t2) += _msh2d->Get_edges_length()[i]*_fv->Get_dt() /_msh2d->Get_triangles_area()[t2] * _fv->Get_dt() * _fv->Get_flux().row(i);
      }

    }
}
}
