#ifndef _MOD_ALGO_CPP

#include "mod_algo.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Algo::Algo()
{}

Algo::~Algo()
{}

void Algo::Initialize(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd x, double eps)
{
    _A = A;
    _b = b;
    _x = x;
    _n = _A.rows();
    _eps = eps;
    _k = 0;
}

const Eigen::VectorXd & Algo::GetX() const
{
  return _x;
}

int Algo::GetK()
{
    return _k;
}

void Algo::Iterative_method(Eigen::MatrixXd M, Eigen::MatrixXd N)
{
    // Definition de l'erreur 
    double error((_A*_x - _b).array().abs().sum());
    // Boucle sur le systeme jusqu'a avoir une erreur inferieur a epsilon
    while (error > _eps)
    {
        _k += 1;
        _x = (M.inverse()*N)*_x + M.inverse()*_b;
        error = (_A*_x - _b).array().abs().sum();
    }
}

void Jacobi::BuildX()
{
    // Matrice M
    MatrixXd M(_A.diagonal().asDiagonal());
    // Matrice N
    MatrixXd N(M - _A);
    // Boucle de la methode iterative
    this->Iterative_method(M, N);
}

void Gauss_Seidel::BuildX()
{
    // Matrice M
    MatrixXd M(_A.triangularView<Lower>());
    // Matrice N
    MatrixXd N(MatrixXd(_A.diagonal().asDiagonal()) - MatrixXd(_A.triangularView<Upper>()));
    // Boucle de la methode iterative
    this->Iterative_method(M, N);
}

Relaxation::Relaxation(double val_w):
    _w(val_w)
{}

void Relaxation::BuildX()
{
    // Matrice M
    MatrixXd M((1/_w - 1)*MatrixXd(_A.diagonal().asDiagonal()) + MatrixXd(_A.triangularView<Lower>()));
    // Matrice N
    MatrixXd N((1/_w)*MatrixXd(_A.diagonal().asDiagonal()) - MatrixXd(_A.triangularView<Upper>()));
    // Boucle de la methode iterative
    this->Iterative_method(M, N);
}

Gradient_Optimal::Gradient_Optimal(int val_kmax):
    _kmax(val_kmax)
{}

void Gradient_Optimal::BuildX()
{
    VectorXd r(_b - _A*_x);
    VectorXd z;
    double alpha;
    while(r.norm() > _eps && _k <= _kmax)
    {
        z = _A * r;
        alpha = r.dot(r)/z.dot(r);
        _x += alpha*r;
        r += -alpha*z;
        _k += 1;
    }
    if (_k > _kmax)
    {
        cout << "Tolerance non atteinte :" << r.norm() << endl;
    }
}

Gradient_Conjugue::Gradient_Conjugue(int val_kmax):
    _kmax(val_kmax)
{}

void Gradient_Conjugue::BuildX()
{
    VectorXd r(_b - _A*_x);
    VectorXd p(r);
    double beta(r.norm());
    VectorXd z;
    double alpha;
    double gamma;
    while(beta > _eps && _k <= _kmax)
    {
        z = _A * p;
        alpha = r.dot(r)/z.dot(p);
        _x += +alpha*p;
        gamma = (r - alpha*z).dot(r - alpha*z)/r.dot(r);
        p = (r - alpha*z) + gamma*p;
        beta = r.norm();
        r += -alpha*z;
        _k += 1;
        cout << (_b - _A*_x).norm() << endl;
    }
    if (_k > _kmax)
    {
        cout << "Tolerance non atteinte : " << r.norm() << endl;
    }
}

#define _MOD_ALGO_CPP
#endif