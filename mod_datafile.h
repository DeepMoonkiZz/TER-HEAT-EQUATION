#ifndef _MOD_DATAFILE_H

#include <cmath>
#include <iostream>

using namespace std;

class Datafile
{
    private:
        // Nombre de maille en x
        int _nx;
        // Nombre de maille en y
        int _ny;
        // x minimum
        double _xmin;
        // x maximum
        double _xmax;
        // y minimum
        double _ymin;
        // y maximum
        double _ymax;
        // pas d'espace en x
        double _dx;
        // pas d'espace en y
        double _dy;
        // pas de temps
        double _dt;
        // valeur de rho
        double _rho;
        // valeur de tmax
        double _tmax;
        
    public:
        // constructeur par defaut
        Datafile();
        // Initialisation de la classe
        void Initialize();
        // Definition de Cp(x) 
        double Cp(double Temp);
        // Definition de D(x)
        double Lambda_harm(double Tplus, double Tmoins);
        // Definition de D(x)
        double Lambda(double T);
        // Definition de T ouest(t) 
        double CL_ouest(double t, double y, string CL);
        // Definition de T est(t)
        double CL_est(double t, double y, string CL);
        // Definition de T nord(t) 
        double CL_nord(double t, double x, string CL);
        // Definition de T sud(t)
        double CL_sud(double t, double x, string CL);
        // Definition de la donnée initiale
        double SolInit(double x, double y);
        // Definition de la solution exacte a l'instant t et au point x
        double SolExact(double t, double x, double y);
        // Definition du terme source en fonction de T
        double Source_term(double T);
        // recuperer le nombre de maille en x
        int GetNX() {return _nx;};
        // recuperer le nombre de maille en y
        int GetNY() {return _ny;};
        // recuperer x minimum
        double GetXmin() {return _xmin;};
        // recuperer x maximum
        double GetXmax() {return _xmax;};
        // recuperer x minimum
        double GetYmin() {return _ymin;};
        // recuperer x maximum
        double GetYmax() {return _ymax;};
        // recuperer le pas en x
        double GetDeltaX() {return _dx;};
        // recuperer le pas en y
        double GetDeltaY() {return _dy;};
        // recuperer le pas en temps
        double GetDeltaT() {return _dt;};
        // recuperer le rho
        double GetRho() {return _rho;};
        // recuperer le tmax
        double GetTmax() {return _tmax;};
};

#define _MOD_DATAFILE_H
#endif