#include "Algo.h"
#include "Schema.h"
#include "Datafile.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;


int main()
{
    std::srand(time(0));
    VectorXd Tn, Texact;
    double error;

    Datafile* df(0);
    Schema* schema(0);
    Algo* algo(0);

    algo = new Gradient_Optimal(1000);
    df = new Datafile();
    schema = new Schema(df, algo);


    // Solution numérique
    schema->Initialize();
    while(schema->GetT() < df->GetTmax())
    {
        schema->Update();
    }    
    Tn = schema->GetTn();


    // Solution exacte;
    Texact.resize(df->GetNX() * df->GetNY());
    for(int i=0; i < df->GetNX(); i++)
	{
        for(int j=0; j < df->GetNY(); j++)
        {
            Texact(j + i * df->GetNX()) = df->SolExact(0., df->GetXmin() + i*df->GetDeltaX(), df->GetYmin() + j*df->GetDeltaY());
        }
	}


    // Résultat 
    error = (Texact - schema->GetTn()).array().abs().sum()/(df->GetNX()*df->GetNY());
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Erreur du schema entre Tn et la solution exacte : " << error << endl;
    cout << "-------------------------------------------------------------------" << endl;

    cout << "Matrice T:" << endl;
    cout << "-------------------------------------------------------------------" << endl;
    for (int j=0; j<df->GetNY(); j++) {
        for (int i=0; i<df->GetNX(); i++) {
            cout << Tn(j + i * df->GetNX()) << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------------" << endl;
        
    cout << "Matrice T exact:" << endl;
    for (int j=0; j<df->GetNY(); j++) {
        for (int i=0; i<df->GetNX(); i++) {
            cout << Texact(j + i * df->GetNX()) << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------------" << endl;
    //cout << schema->GetA()<< endl;


    return 0;
}