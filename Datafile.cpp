#ifndef _DATAFILE_CPP

#include "Datafile.h"
#include <fstream>
#include <cmath>

double pi(acos(-1));

Datafile::Datafile()
{}

void Datafile::Initialize()
{
    ifstream fichier;
    fichier.open("param.txt");
    fichier >> _xmin;
    fichier >> _xmax;
    fichier >> _ymin;
    fichier >> _ymax;
    fichier >> _dx;
    fichier >> _dy;
    fichier >> _dt;
    fichier >> _rho;
    fichier >> _tmax;
    _nx = ceil(1 + (_xmax-_xmin)/_dx);
    _ny = ceil(1 + (_ymax-_ymin)/_dy);
}


// Definition de Cp(x) 
double Datafile::Cp(double Temp)
{
    return 1.;
}

// Definition de la moyenne harmonique de lambda(T)
double Datafile::Lambda_harm(double Tplus, double Tmoins)
{
    return 2/(1/Lambda(Tplus)+1/Lambda(Tmoins));
}

// Definition de lambda(T)
double Datafile::Lambda(double T)
{
    // return 1. + 2.*T;
    return 1.;
}

// Definition de T ouest(t)
double Datafile::CL_ouest(double t, double y, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else
    {
        return SolExact(t, _xmin-_dx, y);
    }
}

// Definition de T est(t)
double Datafile::CL_est(double t, double y, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else
    {
        return SolExact(t, _xmax+_dx, y);
    }
}

// Definition de T sud(t)
double Datafile::CL_sud(double t, double x, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else
    {
        return SolExact(t, x, _ymin-_dy);
    }
}

// Definition de T est(t)
double Datafile::CL_nord(double t, double x, string CL)
{
    if (CL=="neumann")
    {
        return 1.;
    }
    else
    {
        return SolExact(t, x, _ymax+_dy);
    }
}

// Definition de la donnée initiale
double Datafile::SolInit(double x, double y)
{
    // Solution initial constante
    return 1.;

    // Solution initial non constante
    // return sin(pi*2*x/(_xmax-_xmin));
}

double Datafile::SolExact(double t, double x, double y)
{
    // Validation du schéma
    // return 1.-y;

    // Validation de lambda(T(x)) non constant
    // return -1./2. + sqrt(1+8*y)/2.;

    // Validation code 2D
    return sin(2*pi*x)*sin(2*pi*y)*exp(4*Lambda(0)*pi*pi*t);
}

double Datafile::Source_term(double T)
{
    return 0;
}

#define _DATAFILE_CPP
#endif