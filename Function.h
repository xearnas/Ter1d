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
    double Initial_condition_forH(double x);
    double Initial_condition_forU(double x);
    double Source_term(double x, double t);
    double maxdeB(MatrixXd B);
    VectorXd phi(VectorXd theta);
};
