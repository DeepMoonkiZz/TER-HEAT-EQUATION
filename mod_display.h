#include "mod_schema.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

void Save_Solution(double t, Eigen::VectorXd x, Datafile* df, const char *sol);