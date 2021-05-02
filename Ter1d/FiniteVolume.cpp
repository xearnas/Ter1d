#include "FiniteVolume.h"
#include <fstream>
#include <math.h>

using namespace std;
using namespace Eigen;

FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh1D* mesh1d,Mesh2D* mesh2d)
  : _fct(function) , _df(data_file) , _msh1d(mesh1d) , _msh2d(mesh2d)
{
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//------ATTENTION LE FLUX(0) EST EN FAIT f-1/2 (ON DECALE POUR PRENDRE EN COMPTE LES CL)--------------
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Rusanov::Rusanov(Function* function, DataFile* data_file,  Mesh1D* mesh1d,Mesh2D* mesh2d)
: FiniteVolume(function,data_file,mesh1d,mesh2d)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    _B.resize(_df->Get_Nx()+1,2);
    _flux.resize(_df->Get_Nx()+2,2);
  }
  if (_df->Get_mesh_name()=="Mesh2D")
  {
    _B.resize(_msh2d->Get_triangles().size(),3);
    _flux.resize(_msh2d->Get_edges().size(),3);
  }
  cout << "You choose the Rusanov scheme" << endl;
}

LaxF::LaxF(Function* function, DataFile* data_file,  Mesh1D* mesh1d,Mesh2D* mesh2d)
: FiniteVolume(function,data_file,mesh1d,mesh2d)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    _B.resize(_df->Get_Nx()+1,2);
    _flux.resize(_df->Get_Nx()+2,2);
  }
  if (_df->Get_mesh_name()=="Mesh2D")
  {
    _B.resize(_msh2d->Get_triangles().size(),3);
    _flux.resize(_msh2d->Get_triangles().size(),3);
    this->InitializeMatrices();
  }
  cout << "You choose the Lax Frederichs scheme" << endl;
}

FiniteVolume::~FiniteVolume()
{
}

//:::::::::::::::::::::::::::::::::RUSANOV::::::::::::::::::::::::::::::::::::::::

void Rusanov::MajB(MatrixXd Vn)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    for (int i=0 ; i<this->Get_B().rows() ; i++)
    {
      _B(i,0)=abs(Vn(i,1)/Vn(i,0)+sqrt(_df->Get_g()*Vn(i,0)));
      _B(i,1)=abs(Vn(i,1)/Vn(i,0)-sqrt(_df->Get_g()*Vn(i,0)));
    }
    _biundemi=_fct->maxdeB(_B);
    _dt = 0.9*_df->Get_dx()/(2*_biundemi);
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    for (int i=0 ; i<this->Get_B().rows() ; i++)
    {
      _B(i,0)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0));
      _B(i,1)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0) + sqrt(_df->Get_g()*Vn(i,0)));
      _B(i,2)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0) - sqrt(_df->Get_g()*Vn(i,0)));
    }
    _biundemi=_fct->maxdeB(_B);
    _dt = 0.9*_df->Get_dx()/(2*_biundemi);
  }

}

MatrixXd Rusanov::F(VectorXd V)
{
  MatrixXd res;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    res.resize(2,1);
    res(0)=V(1);
    res(1)=V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    res.resize(3,2);
    res.row(0) << V(1) , V(2);
    res.row(1) << V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2. , V(1)*V(2)/V(0);
    res.row(2) << V(1)*V(2)/V(0) , V(2)*V(2)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  }
  return res;
}

MatrixXd Rusanov::Ftild(VectorXd Vi, VectorXd Viplus1,double b, int i)
{
  MatrixXd res;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    res.resize(1,2);
    res(0,0)= 1.0/2.0 * (F(Vi)(0,0) + F(Viplus1)(0,0)) - b/2.0 * (Viplus1(0)-Vi(0));
    res(0,1)= 1.0/2.0 * (F(Vi)(1,0) + F(Viplus1)(1,0)) - b/2.0 * (Viplus1(1)-Vi(1));

  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    MatrixXd Fvt1= F(Vi);
    MatrixXd Fvt2= F(Viplus1);
    res.resize(1,3);
    res(0,0)= 1.0/2.0 * ( Fvt1(0,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(0,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(0)-Vi(0));
    res(0,1)= 1.0/2.0 * ( Fvt1(1,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(1,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(1)-Vi(1));
    res(0,2)= 1.0/2.0 * ( Fvt1(2,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(2,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(2)-Vi(2));
  }
  return res;
}

void Rusanov::MajFlux(MatrixXd Vn)
{
  this->MajB(Vn);
  if (_df->Get_mesh_name() == "Mesh1D")
  {
    for (int i =0; i<_df->Get_Nx()+2;i++)
    {
      if (i==0)
      {
        _flux(0,0)= Vn(0,1);
        _flux(0,1)= Vn(0,1)*Vn(0,1)/Vn(0,0) + _df->Get_g()/2.0 * Vn(0,0)*Vn(0,0);
      }
      else if (i==_df->Get_Nx()+1)
      {
        _flux(_df->Get_Nx()+1,0)= Vn(_df->Get_Nx(),1);
        _flux(_df->Get_Nx()+1,1)= Vn(_df->Get_Nx(),1)*Vn(_df->Get_Nx(),1)/Vn(_df->Get_Nx(),0) + _df->Get_g()/2.0 * Vn(_df->Get_Nx(),0)*Vn(_df->Get_Nx(),0);
      }
      else
      {
        _flux.row(i)=Ftild(Vn.row(i-1),Vn.row(i), max(max(_B(i,0),_B(i,1)),max(_B(i-1,0),_B(i-1,1))),i);
      }
    }
  }
  else if(_df->Get_mesh_name() == "Mesh2D")
  {
    for (int i =0; i<_msh2d->Get_edges().size();i++)
    {
      int t1 = _msh2d->Get_edges()[i].Get_T1();
      int t2 = _msh2d->Get_edges()[i].Get_T2();
      if (t2==-1)
      {
        ///ON EST SUR UN BORD
        _flux.row(i) =  Ftild(Vn.row(t1),Vn.row(t1), max(max(_B(t1,0),_B(t1,1)),max(_B(t1,0),_B(t1,2))),i) ;
      }
      else
      {
        _flux.row(i) =  Ftild(Vn.row(t1),Vn.row(t2), max(max(max(_B(t1,0),_B(t1,1)),max(_B(t2,0),_B(t2,1))) , max (_B(t1,2),_B(t2,2))),i) ;
      }
    }
  }

}

//:::::::::::::::::::::::::::::::::LAX FREDERICHS::::::::::::::::::::::::::::::::::::::::

void LaxF::InitializeMatrices()
{
  _All_Matrices.resize(_msh2d->Get_triangles().size());
  _All_Compteurs.resize(_msh2d->Get_triangles().size());
  for (int i =0; i<_msh2d->Get_triangles().size();++i)
  {
    _All_Matrices[i].resize(4,3);
    _All_Matrices[i].row(0) <<1.0,0.0,0.0;
    _All_Compteurs[i]=1;
  }
  for (int i =0; i<_msh2d->Get_edges().size();i++)
  {
    int t1 = _msh2d->Get_edges()[i].Get_T1();
    int t2 = _msh2d->Get_edges()[i].Get_T2();

    if (t2==-1)
    {
      _All_Matrices[t1].row(_All_Compteurs[t1]) << 1.0 , 0.0 , 0.0;
      _All_Compteurs[t1]+=1;
    }
    else
    {
      _All_Matrices[t1].row(_All_Compteurs[t1]) << 1.0 , _msh2d->Get_triangles_center()(t2,0) - _msh2d->Get_triangles_center()(t1,0) , _msh2d->Get_triangles_center()(t2,1) - _msh2d->Get_triangles_center()(t1,1) ;
      _All_Matrices[t2].row(_All_Compteurs[t2]) << 1.0 , _msh2d->Get_triangles_center()(t1,0) - _msh2d->Get_triangles_center()(t2,0) , _msh2d->Get_triangles_center()(t1,1) - _msh2d->Get_triangles_center()(t2,1) ;
      _All_Compteurs[t1]+=1;
      _All_Compteurs[t2]+=1;
    }
  }
  for (int i =0; i<_msh2d->Get_triangles().size();++i)
  {
    _All_Matrices[i] = (_All_Matrices[i].transpose()*_All_Matrices[i]).inverse()*_All_Matrices[i].transpose();
    _All_Compteurs[i]=1;
  }
  exit(0);
}

void LaxF::MajB(MatrixXd Vn)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    for (int i=0 ; i<this->Get_B().rows() ; i++)
    {
      _B(i,0)=abs(Vn(i,1)/Vn(i,0)+sqrt(_df->Get_g()*Vn(i,0)));
      _B(i,1)=abs(Vn(i,1)/Vn(i,0)-sqrt(_df->Get_g()*Vn(i,0)));
    }
    _biundemi=_fct->maxdeB(_B);
    _dt = 0.9*_df->Get_dx()/(2*_biundemi);
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    for (int i=0 ; i<this->Get_B().rows() ; i++)
    {
      _B(i,0)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0));
      _B(i,1)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0) + sqrt(_df->Get_g()*Vn(i,0)));
      _B(i,2)=abs(Vn(i,1)/Vn(i,0) + Vn(i,2)/Vn(i,0) - sqrt(_df->Get_g()*Vn(i,0)));
    }
    _biundemi=_fct->maxdeB(_B);
    _dt = 0.9*_df->Get_dx()/(2*_biundemi);
  }

}

MatrixXd LaxF::F(VectorXd V)
{
  MatrixXd res;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    res.resize(2,1);
    res(0)=V(1);
    res(1)=V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    res.resize(3,2);
    res.row(0) << V(1) , V(2);
    res.row(1) << V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2. , V(1)*V(2)/V(0);
    res.row(2) << V(1)*V(2)/V(0) , V(2)*V(2)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  }
  return res;
}

MatrixXd LaxF::Ftild(VectorXd Vi, VectorXd Viplus1,double b, int i)
{
  MatrixXd res;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    res.resize(1,2);
    res(0,0)= 1.0/2.0 * (F(Vi)(0,0) + F(Viplus1)(0,0)) - b/2.0 * (Viplus1(0)-Vi(0));
    res(0,1)= 1.0/2.0 * (F(Vi)(1,0) + F(Viplus1)(1,0)) - b/2.0 * (Viplus1(1)-Vi(1));

  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    MatrixXd Fvt1= F(Vi);
    MatrixXd Fvt2= F(Viplus1);
    res.resize(1,3);
    res(0,0)= 1.0/2.0 * ( Fvt1(0,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(0,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(0)-Vi(0));
    res(0,1)= 1.0/2.0 * ( Fvt1(1,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(1,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(1)-Vi(1));
    res(0,2)= 1.0/2.0 * ( Fvt1(2,0) * _msh2d->Get_edges_normal()(i,0) + Fvt1(2,1) * _msh2d->Get_edges_normal()(i,1)) - b * (Viplus1(2)-Vi(2));
  }
  return res;
}

void LaxF::MajFlux(MatrixXd Vn)
{
  this->MajB(Vn);
  if (_df->Get_mesh_name() == "Mesh1D")
  {
    for (int i =0; i<_df->Get_Nx()+2;i++)
    {
      if (i==0)
      {
        _flux(0,0)= Vn(0,1);
        _flux(0,1)= Vn(0,1)*Vn(0,1)/Vn(0,0) + _df->Get_g()/2.0 * Vn(0,0)*Vn(0,0);
      }
      else if (i==_df->Get_Nx()+1)
      {
        _flux(_df->Get_Nx()+1,0)= Vn(_df->Get_Nx(),1);
        _flux(_df->Get_Nx()+1,1)= Vn(_df->Get_Nx(),1)*Vn(_df->Get_Nx(),1)/Vn(_df->Get_Nx(),0) + _df->Get_g()/2.0 * Vn(_df->Get_Nx(),0)*Vn(_df->Get_Nx(),0);
      }

      //-----------------------------------------------------------------------------------------
      /// NOUVELLE CONDITION POUR EVITER LES OUT OF RANGE POUR LE CALCUL DE THETA
      else if (i==_df->Get_Nx())
      {
        _flux.row(i)=Ftild(Vn.row(i-1),Vn.row(i), max(max(_B(i,0),_B(i,1)),max(_B(i-1,0),_B(i-1,1))) ,i);
      }
      //-----------------------------------------------------------------------------------------
      else
      {
        VectorXd theta(2);

        ///------------VERIFIER QUE LA DIVISION DE THETA EST FAISABLE ET FAIRE DE L'ORDRE 2 SI NON---------------------
        if (Vn(i+1,1)-Vn(i,1)<0.001)
        {
          theta(0)=1.0;
          theta(1)=1.0;
        }
        else if (Vn(i+1,0)-Vn(i,0)<0.001)
        {
          theta(0)=1.0;
          theta(1)=1.0;
        }
        else
        {
          theta(0)=(Vn(i,0)-Vn(i-1,0))/(Vn(i+1,0)-Vn(i,0));
          theta(1)=(Vn(i,1)-Vn(i-1,1))/(Vn(i+1,1)-Vn(i,1));
        }

        ///---------------------CREER LES F1 ET F2---------------------
        MatrixXd Fluxordre1(1,2),Fluxordre2(1,2);
        Fluxordre1=Ftild(Vn.row(i-1),Vn.row(i), max(max(_B(i-1,0),_B(i-1,1)),max(_B(i,0),_B(i,1))),i);
        Fluxordre2=Ftild(-Vn.row(i)/2.0 + 3.0/2.0 * Vn.row(i-1),1.0/2.0*(Vn.row(i-1)+Vn.row(i)), max(max(_B(i-1,0),_B(i-1,1)),max(_B(i,0),_B(i,1))),i);

        ///------------------CALCULER LE FLUX EN FONCTION DE THETA-----------------------
        _flux(i,0)=Fluxordre1(0,0) + _fct->phi(theta)(0) * (Fluxordre2(0,0) - Fluxordre1(0,0));
        _flux(i,1)=Fluxordre1(0,1) + _fct->phi(theta)(1) * (Fluxordre2(0,1) - Fluxordre1(0,1));

      }
    }
  }
  else if (_df->Get_mesh_name() == "Mesh2D")
  {
    for (int i =0; i<_msh2d->Get_edges().size();i++)
    {
      int t1 = _msh2d->Get_edges()[i].Get_T1();
      int t2 = _msh2d->Get_edges()[i].Get_T2();
      if (t2==-1)
      {
        ///ON EST SUR UN BORD
        _flux.row(i) =  Ftild(Vn.row(t1),Vn.row(t1), max(max(_B(t1,0),_B(t1,1)),max(_B(t1,0),_B(t1,2))),i) ;
      }
      else
      {
        double distancex = abs( _msh2d->Get_triangles_center()(t1,0) - _msh2d->Get_triangles_center()(t2,0));
        double distancey = abs( _msh2d->Get_triangles_center()(t1,1) - _msh2d->Get_triangles_center()(t2,1));


        MatrixXd Fluxordre1(1,3) , Fluxordre2(1,3);
        Fluxordre1 =  Ftild(Vn.row(t1),Vn.row(t2), max(max(max(_B(t1,0),_B(t1,1)),max(_B(t2,0),_B(t2,1))) , max (_B(t1,2),_B(t2,2))),i) ;
        Fluxordre2 = Ftild((1  +  _All_Matrices[t1](1,0)*distancex  +  _All_Matrices[t1](2,0)*distancey  )*Vn.row(t1),
         (_All_Matrices[t2](1,1)*distancex + _All_Matrices[t2](2,1)*distancey )*Vn.row(t2), max(max(max(_B(t1,0),_B(t1,1)),max(_B(t2,0),_B(t2,1))) , max (_B(t1,2),_B(t2,2))),i) ;
         _flux(i,0)=Fluxordre1(0,0) + _fct->phi(theta)(0) * (Fluxordre2(0,0) - Fluxordre1(0,0));
         _flux(i,1)=Fluxordre1(0,1) + _fct->phi(theta)(1) * (Fluxordre2(0,1) - Fluxordre1(0,1));
         _flux(i,2)=Fluxordre1(0,2) + _fct->phi(theta)(2) * (Fluxordre2(0,2) - Fluxordre1(0,2));
      }
    }
  }
}




//------------------------------------------------------------------------------------------
//::::::::::::::::::::FONCTIONS COMMUNES AUX FLUX:::::::::::::::::::::::::::::::::::::::::::
//------------------------------------------------------------------------------------------

MatrixXd FiniteVolume::Initial_condition()
{
  MatrixXd Vn0;
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    Vn0.resize(this->Get_B().rows(),2);
    for (int i =0; i<this->Get_B().rows() ; i++)
    {
      Vn0(i,0) = _fct->Initial_condition_forH(_msh1d->Get_noeuds()(i,0),0);
      Vn0(i,1) = Vn0(i,0)*_fct->Initial_condition_forU(_msh1d->Get_noeuds()(i,0),0)(0);
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    Vn0.resize(this->Get_B().rows(),3);
    for (int i =0; i<this->Get_B().rows() ; i++)
    {
      Vn0(i,0) = _fct->Initial_condition_forH(_msh2d->Get_triangles_center()(i,0),_msh2d->Get_triangles_center()(i,1));
      Vn0(i,1) = Vn0(i,0)*_fct->Initial_condition_forU(_msh2d->Get_triangles_center()(i,0),_msh2d->Get_triangles_center()(i,1))(0);
      Vn0(i,2) = Vn0(i,0)*_fct->Initial_condition_forU(_msh2d->Get_triangles_center()(i,0),_msh2d->Get_triangles_center()(i,1))(1);
    }
  }
  else
  {
    cout << "Choose an implanted mesh, we cannot calculate the solution vector t t=0";
    exit(0);
  }
  return Vn0;
}


VectorXd FiniteVolume::Source_term(double t)
{
  VectorXd S(this->Get_B().rows());
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    for (int i =0; i<this->Get_B().rows() ; i++)
    {
      S(i) = _fct->Source_term(_msh1d->Get_noeuds()(i,0),0,t);
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    for (int i =0; i<this->Get_B().rows() ; i++)
    {
      S(i) = _fct->Source_term(_msh2d->Get_triangles_center()(i,0),_msh2d->Get_triangles_center()(i,1),t);
    }
  }
  else
  {
    cout << "Choose an implanted mesh, we cannot calculate the source term vector t t=0";
    exit(0);
  }
  return S;
}

//::::::::::::::::::::::::::::::: SAVE SOL GENERAL POUR AVOIR DES GIFS ::::::::::::::::::::::::::::::::::::::

void FiniteVolume::Save_sol(MatrixXd Vn, double t)
{
  if (_df->Get_mesh_name()=="Mesh1D")
  {
    if (t==0.0)
    {
      ofstream mon_fichier;
      mon_fichier.open(_df->Get_results() + "h.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        mon_fichier << _df->Get_xmin() + j*_df->Get_dx() << "  " << Vn(j,0) <<  endl;//<< _df->Get_xmin() + j*_df->Get_dx()
      }
      mon_fichier.close();

      mon_fichier.open(_df->Get_results() + "u.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        mon_fichier << _df->Get_xmin() + j*_df->Get_dx() << "  " <<  Vn(j,1)/Vn(j,0) << endl;
      }
      mon_fichier.close();

      mon_fichier.open(_df->Get_results() + "hu.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        mon_fichier << _df->Get_xmin() + j*_df->Get_dx() << "  " << Vn(j,1) << endl;
      }
      mon_fichier.close();
    }
    else
    {
      /// ------------CREATION DU FICHIER H -------------------------
      fstream mon_fichier;
      fstream mon_fichier_prime;
      string ligne;
      mon_fichier.open(_df->Get_results() + "h.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "hprime.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,0) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "h.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "hprime.txt",ios::in);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      /// ------------CREATION DU FICHIER U -------------------------


      mon_fichier.open(_df->Get_results() + "u.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "uprime.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,1)/Vn(j,0) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "u.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "uprime.txt",ios::in);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();


      /// ------------CREATION DU FICHIER HU -------------------------


      mon_fichier.open(_df->Get_results() + "hu.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "huprime.txt",ios::out);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,1) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "hu.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "huprime.txt",ios::in);
      for (int j=0;j<_df->Get_Nx()+1;j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();
    }
  }
  else if (_df->Get_mesh_name()=="Mesh2D")
  {
    if (t==0.0)
    {
      ofstream mon_fichier;
      mon_fichier.open(_df->Get_results() + "h.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        mon_fichier << _msh2d->Get_triangles_center()(j,0) << "  " << _msh2d->Get_triangles_center()(j,1) << "  " << Vn(j,0) <<  endl;//<< _df->Get_xmin() + j*_df->Get_dx()
      }
      mon_fichier.close();

      mon_fichier.open(_df->Get_results() + "hux.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        mon_fichier << _msh2d->Get_triangles_center()(j,0) << "  " << _msh2d->Get_triangles_center()(j,1) << "  " << Vn(j,1) << endl;
      }
      mon_fichier.close();

      mon_fichier.open(_df->Get_results() + "huy.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        mon_fichier << _msh2d->Get_triangles_center()(j,0) << "  " << _msh2d->Get_triangles_center()(j,1) << "  " << Vn(j,2) << endl;
      }
      mon_fichier.close();
    }
    else
    {
      /// ------------CREATION DU FICHIER H -------------------------
      fstream mon_fichier;
      fstream mon_fichier_prime;
      string ligne;
      mon_fichier.open(_df->Get_results() + "h.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "hprime.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,0) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "h.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "hprime.txt",ios::in);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();


      /// ------------CREATION DU FICHIER HUX -------------------------


      mon_fichier.open(_df->Get_results() + "hux.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "huxprime.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,1) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "hux.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "huxprime.txt",ios::in);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      /// ------------CREATION DU FICHIER HUY -------------------------


      mon_fichier.open(_df->Get_results() + "huy.txt", ios::in);
      mon_fichier_prime.open(_df->Get_results() + "huyprime.txt",ios::out);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier,ligne);
        mon_fichier_prime << ligne << "  " << Vn(j,2) <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();

      mon_fichier.open(_df->Get_results() + "huy.txt", ios::out);
      mon_fichier_prime.open(_df->Get_results() + "huyprime.txt",ios::in);
      for (int j=0;j<_msh2d->Get_triangles().size();j++)
      {
        getline(mon_fichier_prime,ligne);
        mon_fichier << ligne <<  endl;
      }
      mon_fichier.close();
      mon_fichier_prime.close();
    }
  }
}
