#include "DataFile.h"

using namespace std;
using namespace Eigen;
#pragma once


class Function
{
private:
    // Pointeur de la classe DataFile pour récupérer toutes les
    // valeurs de paramètres
    const DataFile* _df;
    // Diffusion coefficient
    const double _g;

public: // Méthodes et opérateurs de la classe
    Function(DataFile* data_file);
    double Initial_condition_forH(double x, double y);
    VectorXd Initial_condition_forU(double x, double y);
    double Source_term(double x, double y, double t);
    double maxdeB(MatrixXd B);
    double Produit_Scalaire(VectorXd U, VectorXd V);
    VectorXd phi(VectorXd theta);
};
