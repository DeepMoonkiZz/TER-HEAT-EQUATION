#ifndef _ALGO_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

class Algo
{
    protected:
        // Matrice A
        Eigen::MatrixXd _A;
        // Matrice x et b
        Eigen::VectorXd _x, _b;
        // Dimension de x
        int _n;
        // Epsilon
        double _eps;
        // Nombre d'iteration
        int _k;

    private:

    public:
        // Constructeur
        Algo();
        // Destructeur
        virtual ~Algo();
        // Initialisation de A, b et x 
        void Initialize(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd x, double eps);
        // Permet de recuperer le vecteur x
        const Eigen::VectorXd & GetX() const;
        // Permet de recuperer le nombre d'iteration
        int GetK();
        // Execute l'algorithme voulut  
        virtual void BuildX(){}
        // Algo generique d'une methode iterative
        void Iterative_method(Eigen::MatrixXd M, Eigen::MatrixXd N);
};

class Jacobi : public Algo
{
    public:
        void BuildX();
};

class Gauss_Seidel : public Algo
{
    public:
        void BuildX();
};

class Relaxation : public Algo
{   
    private:
        // Coefficient de relaxation
        double _w;
        
    public:
        // Constructeur relaxation
        Relaxation(double w);
        void BuildX();
};

class Gradient_Optimal : public Algo
{
    private:
        // Tolerance max accepter
        int _kmax;
    
    public:
        Gradient_Optimal(int kmax);
        void BuildX();
};

class Gradient_Conjugue : public Algo
{
    private:
        int _kmax;
    public:
        Gradient_Conjugue(int kmax);
        void BuildX();
};

#define _ALGO_H
#endif