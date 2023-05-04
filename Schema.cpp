#ifndef _SCHEMA_CPP

#include "Schema.h"
#include "Algo.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

Schema::Schema(Datafile* df, Algo* algo):
    _df(df), _algo(algo)
{}

void Schema::Initialize()
{
    _df->Initialize();
	_dt = _df->GetDeltaT(), _dx = _df->GetDeltaX(), _dy = _df->GetDeltaY(), _rho = _df->GetRho(), _nx = _df->GetNX(), _ny = _df->GetNY(), _t = 0.;
	_Tn.resize(_nx*_ny);
	_A.resize(_nx * _ny, _nx * _ny);
    for(int i=0; i<_nx; i++)
	{
		for (int j=0; j<_ny; j++)
		{
			_Tn(j + i * _nx) = _df->SolInit(_df->GetXmin() + i*_df->GetDeltaX(), _df->GetYmin() + j*_df->GetDeltaY());
		}
	}
}

void Schema::Update()
{
	_t += _dt;
	_b = _Tn;
    this->Build_A_and_b();
	_algo->Initialize(this->GetA(), _b, _Tn, pow(10, -5));
    _algo->BuildX();
	_Tn = _algo->GetX();
	cout << int(_t/_df->GetTmax()*100) << "%" << endl;
}

void Schema::Build_A_and_b()
{
	double C;
	string neu("neumann"), diri("diriclet");
	
	// Build A and BC of b
	for (int i=0; i < _nx; i++)
	{
		for (int j=0; j< _ny; j++)
		{
			C = _df->GetDeltaT() / (_df->GetRho() * _df->Cp(_Tn(this->Coord_mat_to_vect(i, j))));
			
			// Tij en n+1
			_triplets.push_back({Coord_mat_to_vect(i, j), Coord_mat_to_vect(i, j), 1});
			// Flux F en i
			this->Flux_F(i, j, C / (_df->GetDeltaX() * _df->GetDeltaX()), diri);
			// Flux G en j
			this->Flux_G(i, j, C / (_df->GetDeltaY() * _df->GetDeltaY()), diri);
		}
	}
	_A.setFromTriplets(_triplets.begin(), _triplets.end());
	_triplets.clear();

	// Build b
	for (int i=0; i < _nx; i++) {
    	for (int j=0; j < _ny; j++) {
            _b(Coord_mat_to_vect(i, j)) += _df->Source_term(_Tn(Coord_mat_to_vect(i, j)));
        }
    }
}


void Schema::Flux_F(int i, int j, double Cx, string CL)
{
	double Tij, Timj, Tipj;
	double CLe, CLo;
	int Coord_ij, Coord_imj, Coord_ipj;

	Coord_ij = Coord_mat_to_vect(i, j);
	Coord_imj = Coord_mat_to_vect(i - 1, j);
	Coord_ipj = Coord_mat_to_vect(i + 1, j);


	// Condition limite ouest
	if (i==0) {
		Tij = _Tn(Coord_ij);
		Tipj = _Tn(Coord_ipj);
		CLo = _df->CL_ouest(_t, _df->GetYmin() + j * _df->GetDeltaY(), CL);

		if (CL=="diriclet") {
			// Condition limite
			_b(Coord_ij) += Cx * _df->Lambda_harm(Tij, CLo) * CLo;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, CLo)});
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Tipj)});
			_triplets.push_back({Coord_ij, Coord_ipj, - Cx * _df->Lambda_harm(Tij, Tipj)});
		}
		else if (CL=="neumann") {
			// Condition limite
			_b(Coord_ij) += _df->GetDeltaX() * Cx * _df->Lambda_harm(Tij, CLo) * CLo;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Tipj)});
			_triplets.push_back({Coord_ij, Coord_ipj, - Cx * _df->Lambda_harm(Tij, Tipj)});
		}
	}

	// Condition limite est
	else if (i==_nx-1) {
		Tij = _Tn(Coord_ij);
		Timj = _Tn(Coord_imj);
		CLe = _df->CL_est(_t, _df->GetYmin() + j * _df->GetDeltaY(), CL);

		if (CL=="diriclet") {
			// Condition limite
			_b(Coord_ij) += Cx * _df->Lambda_harm(CLe, Tij) * CLe;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, CLe)});
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Timj)});
			_triplets.push_back({Coord_ij, Coord_imj, - Cx * _df->Lambda_harm(Tij, Timj)});
		}
		else if (CL=="neumann") {
			// Condition limite
			_b(Coord_ij) += _df->GetDeltaX() * Cx * _df->Lambda_harm(CLe, Tij) * CLe;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Timj)});
			_triplets.push_back({Coord_ij, Coord_imj, - Cx * _df->Lambda_harm(Tij, Timj)});
		}
	}

	// Calcul intérieur sans condition limite
	else {
		Tij = _Tn(Coord_ij);
		Timj = _Tn(Coord_imj);
		Tipj = _Tn(Coord_ipj);
		
		_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Timj)});
		_triplets.push_back({Coord_ij, Coord_imj, - Cx * _df->Lambda_harm(Tij, Timj)});
		_triplets.push_back({Coord_ij, Coord_ij, Cx * _df->Lambda_harm(Tij, Tipj)});
		_triplets.push_back({Coord_ij, Coord_ipj, - Cx * _df->Lambda_harm(Tij, Tipj)});
	}
}


void Schema::Flux_G(int i, int j, double Cy, string CL)
{
	double Tij, Tijm, Tijp;
	double CLn, CLs;
	int Coord_ij, Coord_ijm, Coord_ijp;

	Coord_ij = Coord_mat_to_vect(i, j);
	Coord_ijm = Coord_mat_to_vect(i, j - 1);
	Coord_ijp = Coord_mat_to_vect(i, j + 1);

	// Condition limite sud
	if (j==0) {
		Tij = _Tn(Coord_ij);
		Tijp = _Tn(Coord_ijp);
		CLs = _df->CL_sud(_t, _df->GetXmin() + i * _df->GetDeltaX(), CL);

		if (CL=="diriclet") {
			// Condition limite
			_b(Coord_ij) += Cy * _df->Lambda_harm(Tij, CLs) * CLs;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, CLs)});
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, Tijp)});
			_triplets.push_back({Coord_ij, Coord_ijp, - Cy * _df->Lambda_harm(Tij, Tijp)});
		}
		else if (CL=="neumann") {
			// Condition limite
			_b(Coord_ij) += _df->GetDeltaY() * Cy * _df->Lambda_harm(Tij, CLs) * CLs;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, Tijp)});
			_triplets.push_back({Coord_ij, Coord_ijp, - Cy * _df->Lambda_harm(Tij, Tijp)});
		}
	}

	// Condition limite nord
	else if (j==_ny-1) {	
		Tij = _Tn(Coord_ij);
		Tijm = _Tn(Coord_ijm);
		CLn = _df->CL_nord(_t, _df->GetXmin() + i * _df->GetDeltaX(), CL);

		if (CL=="diriclet") {
			// Condition limite
			_b(Coord_ij) += Cy * _df->Lambda_harm(CLn, Tij) * CLn;

			// Calcul intérieur 
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, CLn)});
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, Tijm)});
			_triplets.push_back({Coord_ij, Coord_ijm, - Cy * _df->Lambda_harm(Tij, Tijm)});
		}
		else if (CL=="neumann") {
			// Condition limite
			_b(Coord_ij) += _df->GetDeltaY() * Cy * _df->Lambda_harm(CLn, Tij) * CLn;

			// Calcul intérieur
			_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tij, Tijm)});
			_triplets.push_back({Coord_ij, Coord_ijm, - Cy * _df->Lambda_harm(Tij, Tijm)});
		}
	}

	// Calcul intérieur sans condition limite
	else {
		Tij = _Tn(Coord_ij);
		Tijm = _Tn(Coord_ijm);
		Tijp = _Tn(Coord_ijp);
		
		_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tijm, Tij)});
		_triplets.push_back({Coord_ij, Coord_ijm, - Cy * _df->Lambda_harm(Tijm, Tij)});
		_triplets.push_back({Coord_ij, Coord_ij, Cy * _df->Lambda_harm(Tijp, Tij)});
		_triplets.push_back({Coord_ij, Coord_ijp, - Cy * _df->Lambda_harm(Tijp, Tij)});
	}
}


Vector2i Schema::Coord_vect_to_mat(int N)
{
	Vector2i coord;
	coord(0) = N % _df->GetNX();
	coord(1) = N / _df->GetNX();
	return coord;
}


int Schema::Coord_mat_to_vect(int i, int j)
{
	int N;
	N = j + i * _df->GetNX();
	return N;
}


VectorXd Schema::GetTn()
{
	return _Tn;
}


#define _SCHEMA_CPP
#endif