#ifndef _MOD_DATAFILE_CPP

#include "mod_datafile.h"
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
    fichier >> _nx;
    fichier >> _ny;
    fichier >> _rho;
    fichier >> _dt;
    fichier >> _tmax;
    _dx = (_xmax-_xmin)/(_nx-1);
    _dy = (_ymax-_ymin)/(_ny-1);
}


// Definition de Cp(x) 
double Datafile::Cp(double Temp)
{
    return 1.;
}

// Definition de lambda(T)
double Datafile::Lambda(double T)
{
    // Cas test 1D, diffu et 2D
    return 1.;

    // Cas test lambda non constant
    // return 1. + 2.*T;
}

// Definition de la moyenne harmonique de lambda(T)
double Datafile::Lambda_harm(double Tplus, double Tmoins)
{
    return 2/(1/Lambda(Tplus)+1/Lambda(Tmoins));
}

// Definition de la donnée initiale
double Datafile::SolInit(double x, double y)
{
    // Validation code 1D, lambda(T(x)) non constant et diffu
    // return 0.;

    // Validation code 2D
    return sin(2*pi*x)*sin(2*pi*y)*exp(0);
}

double Datafile::SolExact(double t, double x, double y)
{
    // Validation code 1D
    // return y;

    // Validation de lambda(T(x)) non constant
    // return -1./2. + sqrt(1+8*y)/2.;

    // Validation code diffu
    // return 1-y;

    // Validation code 2D
    return sin(2*pi*x/(_xmax-_xmin))*sin(2*pi*y/(_ymax-_ymin))*exp(-8*pi*pi*t/(_xmax-_xmin)/(_ymax-_ymin));
}

double Datafile::Source_term(double x, double y, double t)
{
    // Validation cas test
    return 0;
}

// Definition de T ouest(t)
double Datafile::CL_ouest(double t, double y, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else if (CL=="diriclet")
    {
        return SolExact(t, _xmin-_dx, y);
    }
    else
    {
        cout << "No BC found" << endl;
        return 0;
    }
}

// Definition de T est(t)
double Datafile::CL_est(double t, double y, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else if (CL=="diriclet")
    {
        return SolExact(t, _xmax+_dx, y);
    }
    else
    {
        cout << "No BC found" << endl;
        return 0;
    }
}

// Definition de T sud(t)
double Datafile::CL_sud(double t, double x, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else if (CL=="diriclet")
    {
        return SolExact(t, x, _ymin-_dy);
    }
    else
    {
        cout << "No BC found" << endl;
        return 0;
    }
}

// Definition de T est(t)
double Datafile::CL_nord(double t, double x, string CL)
{
    if (CL=="neumann")
    {
        return 0.;
    }
    else if (CL=="diriclet")
    {
        return SolExact(t, x, _ymax+_dy);
    }
    else
    {
        cout << "No BC found" << endl;
        return 0;
    }
}

#define _MOD_DATAFILE_CPP
#endif