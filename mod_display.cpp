#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>

#include "mod_display.h"

using namespace std;
using namespace Eigen;

void Save_Solution(double t, Eigen::VectorXd x, Datafile* df, const char *sol)
{
    char buffer[100];
    if (sol=="num") {
        // Validation 1D
        // std::snprintf(buffer, sizeof(buffer), "Solutions/test_1D/validation_1D/valid_1D_%.3f.dat", t);

        // Validation lambda
        std::snprintf(buffer, sizeof(buffer), "Solutions/test_lambda/validation_lambda/valid_lambda_%.3f.dat", t);

        // Validation diffusion
        // std::snprintf(buffer, sizeof(buffer), "Solutions/test_diffu/validation_diffu/valid_diffu_%.3f.dat", t);

        // Validation 2D
        std::snprintf(buffer, sizeof(buffer), "Solutions/test_2D/validation_2D/valid_2D_%.3f.dat", t);
    }
    else if (sol=="exact") { 
        // Validation 1D
        // std::snprintf(buffer, sizeof(buffer), "Solutions/test_1D/exact_1D/exact_1D_%.3f.dat", t);

        // Validation lambda
        // std::snprintf(buffer, sizeof(buffer), "Solutions/test_lambda/exact_lambda/exact_lambda_%.3f.dat", t);

        // Validation diffusion
        // std::snprintf(buffer, sizeof(buffer), "Solutions/test_diffu/exact_diffu/exact_diffu_%.3f.dat", t);

        // Validation 2D
        std::snprintf(buffer, sizeof(buffer), "Solutions/test_2D/exact_2D/exact_2D_%.3f.dat", t);
    }

    char* fichier = new char[std::strlen(buffer) + 1];
    std::strcpy(fichier, buffer);
    
    // Écriture du fichier
    std::ofstream file_out;
    file_out.open(fichier);

    for(int j = df->GetNY()-5; j >= df->GetNY()-5; j--){
        for(int i = 0; i < df->GetNX(); i++){
            file_out << i*df->GetDeltaX() << " " << j*df->GetDeltaY() << " " << x(j + i * df->GetNX()) << endl;
        }
    }
    file_out.close();
}