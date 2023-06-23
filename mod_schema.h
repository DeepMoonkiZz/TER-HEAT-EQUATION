#ifndef _MOD_SCHEMA_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "mod_datafile.h"
#include "mod_algo.h"
#include <vector>

using namespace std;
using namespace Eigen;

class Schema
{
    private:
        // Definition de la matrice A
        SparseMatrix<double> _A;
        // Pointeur vers une classe datafile
        Datafile* _df;
        // Pointeur vers le système de résolution du schéma
        Algo* _algo;
        // Valeur de la temperature a l'instant n
        VectorXd _Tn;
        // Valeur de la temperature temp a l'instant n
        VectorXd _b;
        // Triplets pour construire A
        vector<Triplet<double>> _triplets;
        // Valeur de t, dt, dx, rho
        double _t, _dt, _dx, _dy, _rho;
        // Valeur de n
        int _nx, _ny;

    public:
        // Constructeur de schema
        Schema(Datafile* df, Algo* algo);
        // Initialisation de la classe schema
        void Initialize();
        // Update la temperature un dt plus tard
        void Update();
        // Construction de A
        void Build_A_and_b();
        // Construction des flux F
        void Flux_F(int i, int j, double Cx, string CLO, string CLE);
        // Construction des flux G
        void Flux_G(int i, int j, double Cy, string CLS, string CLN);
        // Coordonnée du plan vers le vecteur 
        Vector2i Coord_vect_to_mat(int N);
        // Coordonnée du vecteur vers le plan 
        int Coord_mat_to_vect(int i, int j);
        // Recuperation de A
        SparseMatrix<double> GetA() {return _A;};
        // Recuperation du temps t
        double GetT() {return _t;};
        // Recuperation de Tn
        VectorXd GetTn();
};

#define _MOD_SCHEMA_H
#endif