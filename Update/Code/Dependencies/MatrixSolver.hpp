#include "Matrix.hpp"

#include <vector>
#include <algorithm>
#include <Eigen/Sparse>
#include <chrono>
#include <sstream>
#include <fstream>
#include <string>


class MatrixSolver {
public:
    MatrixSolver(Eigen::SparseMatrix<double> sMat, std::vector<double>b,double w, double tolerance = 1e-6, size_t maxIterations = 1e7, double guess = 0, int M = 5);
    bool solveWithGaussSeidel();

    std::vector<double> m_x;
    int M;
private:
    const Eigen::SparseMatrix<double> m_sMat;
    std::vector<double>m_b;
    const double m_tolerance;
    const size_t m_maxIterations;
    double m_guess;
    double m_w;
    std::stringstream statistics;

};

