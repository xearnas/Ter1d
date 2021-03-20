#include "FiniteVolume.h"
#include <fstream>
#include <math.h>

using namespace std;
using namespace Eigen;

FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh* mesh)
  : _fct(function) , _df(data_file) , _msh(mesh)
{
}

Rusanov::Rusanov(Function* function, DataFile* data_file, Mesh* mesh)
: FiniteVolume(function,data_file,mesh)
{
  _B.resize(_df->Get_Nx()+1,2);
  _flux.resize(_df->Get_Nx(),2);
  cout << "You choose the Rusanov scheme" << endl;
  cout << "Finite volume functions are ready to use" << endl;
}

LaxF::LaxF(Function* function, DataFile* data_file, Mesh* mesh)
: FiniteVolume(function,data_file,mesh)
{
  _B.resize(_df->Get_Nx()+1,2);
  _flux.resize(_df->Get_Nx(),2);
  cout << "You choose the Lax Frederichs scheme" << endl;
  cout << "Finite volume functions are ready to use" << endl;
}

FiniteVolume::~FiniteVolume()
{
}

MatrixXd FiniteVolume::Initial_condition()
{

  MatrixXd Vn0(_df->Get_Nx()+1,2);
  for (int i =0; i<_df->Get_Nx()+1 ; i++)
  {
    Vn0(i,0) = _fct->Initial_condition_forH(_msh->Get_noeuds()(i));
    Vn0(i,1) = Vn0(i,0)*_fct->Initial_condition_forU(_msh->Get_noeuds()(i));
  }

  return Vn0;
}


VectorXd FiniteVolume::Source_term(double t)
{
  VectorXd S(_df->Get_Nx()+1);
  for (int i =0; i<_df->Get_Nx()+1 ; i++)
  {
    S(i) = _fct->Source_term(_msh->Get_noeuds()(i),t);
  }
  return S;
}

void FiniteVolume::Save_sol(MatrixXd Vn, double t)
{

  if (t==0.0)
  {
    ofstream mon_fichier;
    mon_fichier.open(_df->Get_results() + "h.txt",ios::out);
    for (int j=0;j<_df->Get_Nx()+1;j++)
    {
      mon_fichier << _df->Get_xmin() + j*_df->Get_dx() << "  " << Vn(j,0) <<  endl;
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

void Rusanov::MajB(MatrixXd Vn)
{
  for (int i=0 ; i<_df->Get_Nx()+1 ; i++)
  {
    _B(i,0)=abs(Vn(i,1)/Vn(i,0)+sqrt(_df->Get_g()*Vn(i,0)));
    _B(i,1)=abs(Vn(i,1)/Vn(i,0)-sqrt(_df->Get_g()*Vn(i,0)));
  }
  _biundemi=_fct->maxdeB(_B);
  _dt = 0.9*_df->Get_dx()/(2*_biundemi);
}

VectorXd Rusanov::F(VectorXd V)
{
  VectorXd res(2);
  res(0)=V(1);
  res(1)=V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  return res;
}

MatrixXd Rusanov::Ftild(VectorXd Vi, VectorXd Viplus1,double b)
{
  MatrixXd res(1,2);
  res(0,0)= - b/2.0 * (Viplus1(0)-Vi(0));
  res(0,1)= - b/2.0 * (Viplus1(1)-Vi(1));
  res.row(0)+=1.0/2.0 * (F(Vi) + F(Viplus1));
  return res;
}

void Rusanov::MajFlux(MatrixXd Vn)
{
  this->MajB(Vn);


  for (int i =0; i<_df->Get_Nx();i++)
  {
    // CL continue (même entre 0 et 1 et Nx et Nx-1)
    if (i==0)
    {
      _flux(0,0)= Vn(0,1);
      _flux(0,1)= Vn(0,1)*Vn(0,1)/Vn(0,0) + _df->Get_g()/2.0 * Vn(0,0)*Vn(0,0);
    }
    else if (i==_df->Get_Nx()-1)
    {
      _flux(_df->Get_Nx()-1,0)= Vn(_df->Get_Nx()-1,1);
      _flux(_df->Get_Nx()-1,1)= Vn(_df->Get_Nx()-1,1)*Vn(_df->Get_Nx()-1,1)/Vn(_df->Get_Nx()-1,0) + _df->Get_g()/2.0 * Vn(_df->Get_Nx()-1,0)*Vn(_df->Get_Nx()-1,0);
    }
    else
    {
      _flux.row(i)=Ftild(Vn.row(i),Vn.row(i+1), max(max(_B(i+1,0),_B(i+1,1)),max(_B(i,0),_B(i,1))) );
      //_flux(i,0)=1.0/2.0 * ( Vn(i+1,1)+Vn(i,1) - max(max(_B(i,0),_B(i,1)),max(_B(i+1,0),_B(i+1,1)))*(Vn(i+1,0) - Vn(i,0)) );
      //_flux(i,1)=1.0/2.0 * ( Vn(i+1,1)*Vn(i+1,1)/Vn(i+1,0) + _df->Get_g()/2.0 * Vn(i+1,0)*Vn(i+1,0) + Vn(i,1)*Vn(i,1)/Vn(i,0) + _df->Get_g()/2.0 * Vn(i,0)*Vn(i,0) - max(max(_B(i,0),_B(i,1)),max(_B(i+1,0),_B(i+1,1))) * ( Vn(i+1,1) - Vn(i,1) ) );

    }
  }

}

void LaxF::MajB(MatrixXd Vn)
{
  for (int i=0 ; i<_df->Get_Nx()+1 ; i++)
  {
    _B(i,0)=abs(Vn(i,1)/Vn(i,0)+sqrt(_df->Get_g()*Vn(i,0)));
    _B(i,1)=abs(Vn(i,1)/Vn(i,0)-sqrt(_df->Get_g()*Vn(i,0)));
  }
  _biundemi=_fct->maxdeB(_B);
  _dt = 0.40*_df->Get_dx()/(2*_biundemi);
}

void LaxF::MajFlux(MatrixXd Vn)
{
  this->MajB(Vn);


  for (int i =0; i<_df->Get_Nx();i++)
  {
    // CL continue (même entre 0 et 1 et Nx et Nx-1)
    if (i==0)
    {
      _flux(0,0)= Vn(0,1);
      _flux(0,1)= Vn(0,1)*Vn(0,1)/Vn(0,0) + _df->Get_g()/2.0 * Vn(0,0)*Vn(0,0);
    }
    else if (i==_df->Get_Nx()-1)
    {
      _flux(_df->Get_Nx()-1,0)= Vn(_df->Get_Nx()-1,1);
      _flux(_df->Get_Nx()-1,1)= Vn(_df->Get_Nx()-1,1)*Vn(_df->Get_Nx()-1,1)/Vn(_df->Get_Nx()-1,0) + _df->Get_g()/2.0 * Vn(_df->Get_Nx()-1,0)*Vn(_df->Get_Nx()-1,0);
    }
    else
    {

      VectorXd theta(2);
      if ( (abs(Vn(i+1,0)-Vn(i,0))<0.01) )
      {
        theta(0)=1.0;
        theta(1)=1.0;
      }
      else if ( (abs(Vn(i+1,1)-Vn(i,1))<0.01) )
      {
        theta(0)=1.0;
        theta(1)=1.0;
      }
      else
      {
        theta(0)=(Vn(i,0)-Vn(i-1,0))/(Vn(i+1,0)-Vn(i,0));
        theta(1)=(Vn(i,1)-Vn(i-1,1))/(Vn(i+1,1)-Vn(i,1));
      }

      MatrixXd Fluxordre1(1,2),Fluxordre2(1,2);
      Fluxordre1=Ftild(Vn.row(i),Vn.row(i+1), max(max(_B(i,0),_B(i,1)),max(_B(i+1,0),_B(i+1,1))));
      Fluxordre2=Ftild(-Vn.row(i+1)/2.0 + 3.0/2.0 * Vn.row(i),1.0/2.0*(Vn.row(i)+Vn.row(i-1)), max(max(_B(i,0),_B(i,1)),max(_B(i+1,0),_B(i+1,1))));
      if (std::isnan(theta(0)))
      {
        cout << "c'est theta(0) à i=" << i << endl;
        cout << Vn(i+1,0)-Vn(i,0) << endl;
        cout << Vn(i+1,0) << endl;
        cout << Vn(i,0) << endl;
        cout << " et bah theta est NaN ben joué mec " << endl;
        exit(0);
      }
      if (std::isnan(theta(1)))
      {
        cout << "c'est theta(1 à i=" << i << endl;
        cout << Vn(i+1,1)-Vn(i,1) << endl;
        cout << Vn(i+1,1) << endl;
        cout << Vn(i,1) << endl;
        cout << " et bah theta est NaN ben joué mec " << endl;
        exit(0);
      }
      _flux(i,0)=Fluxordre1(0,0) + _fct->phi(theta)(0) * (Fluxordre2(0,0) - Fluxordre1(0,0));
      _flux(i,1)=Fluxordre1(0,1) + _fct->phi(theta)(1) * (Fluxordre2(0,1) - Fluxordre1(0,1));
    }
  }
}

VectorXd LaxF::F(VectorXd V)
{
  VectorXd res(2);
  res(0)=V(1);
  res(1)=V(1)*V(1)/V(0) + _df->Get_g()* V(0)*V(0)/2.;
  return res;
}

MatrixXd LaxF::Ftild(VectorXd Vi, VectorXd Viplus1,double b)
{
  MatrixXd res(1,2);
  res(0,0)= - b/2.0 * (Viplus1(0)-Vi(0));
  res(0,1)= - b/2.0 * (Viplus1(1)-Vi(1));
  res.row(0)+=1.0/2.0 * (F(Vi) + F(Viplus1));
  return res;
}
