#include "mod_algo.h"
#include "mod_schema.h"
#include "mod_datafile.h"
#include "mod_display.h"

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Eigen;


int main()
{
    cout.precision(2);   
    std::srand(time(0));
    VectorXd Tn, Texact;
    double error;

    Datafile* df(0);
    Schema* schema(0);
    Algo* algo(0);

    algo = new Gradient_Optimal(1000);
    df = new Datafile();
    schema = new Schema(df, algo);

    const char *sol = "num", *exact = "exact";

    // ----------------------------------------------------------- //

    // Solution numérique

    // ----------------------------------------------------------- //

    schema->Initialize();
    Save_Solution(schema->GetT(), schema->GetTn(), df, sol);
    while(schema->GetT() < df->GetTmax())
    {
        schema->Update();
        Save_Solution(schema->GetT(), schema->GetTn(), df, sol);
    }    
    Tn = schema->GetTn();
    Save_Solution(schema->GetT(), schema->GetTn(), df, sol);



    // ----------------------------------------------------------- //

    // Soution exacte

    // ----------------------------------------------------------- //

    double t(0);
    Texact.resize(df->GetNX() * df->GetNY());
    while(t < df->GetTmax())
    {
        for(int i=0; i < df->GetNX(); i++)
        {
            for(int j=0; j < df->GetNY(); j++)
            {
                Texact(j + i * df->GetNX()) = df->SolExact(t, df->GetXmin() + i*df->GetDeltaX(), df->GetYmin() + j*df->GetDeltaY());
            }
        }
        Save_Solution(t, Texact, df, exact);
        t += df->GetDeltaT();
    }    
    Save_Solution(t, Texact, df, exact);



    // ----------------------------------------------------------- //

    // Calcul de l'erreur

    // ----------------------------------------------------------- //

    error = (Texact - schema->GetTn()).array().abs().sum()/(df->GetNX()*df->GetNY());
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Erreur du schema entre Tn et la solution exacte : " << error << endl;
    cout << "-------------------------------------------------------------------" << endl;



    // ----------------------------------------------------------- //

    // Affichage dans le terminal

    // ----------------------------------------------------------- //

    /*

    cout << "Matrice T:" << endl;
    cout << "-------------------------------------------------------------------" << endl;
    for (int j=df->GetNY() - 1; j>=0; j--) {
        for (int i=0; i<df->GetNX(); i++) {
            cout << fixed << Tn(j + i * df->GetNX()) << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------------" << endl;
        
    cout << "Matrice T exact:" << endl;
    for (int j=df->GetNY() - 1; j>=0; j--) {
        for (int i=0; i<df->GetNX(); i++) {
            cout << fixed << Texact(j + i * df->GetNX()) << " ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------------" << endl;
    //cout << schema->GetA()<< endl;

    */

    return 0;
}