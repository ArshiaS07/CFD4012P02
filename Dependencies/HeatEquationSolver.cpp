/*
 * HeatEquationSolver.cpp
 *
 *  Created on: Mar 26, 2023
 *      Author: ArshiaSaffari
 */

#include "HeatEquationSolver.hpp"
#include <chrono>
std::ostream& operator<<(std::ostream& os, const std::vector<double>& inVec)
{
	for (size_t i = 0; i < inVec.size(); i++)
		os << inVec.at(i) << "\n";
	return os;
}



namespace ArshiaCFD {
	HeatEquationSolver::HeatEquationSolver(Range2D ComputeDomain,
		Index2D NCount,double m_W) :
		m_computeDomain(ComputeDomain), m_nodeCount(NCount),
		m_eqMatDims((m_nodeCount.xcount - 2)* (m_nodeCount.ycount - 2)),
		m_eqSparseMat(m_eqMatDims, m_eqMatDims),
		m_W(m_W) {
		m_nodeSpacing.x = (m_computeDomain.xrange.second - m_computeDomain.xrange.first) / (m_nodeCount.xcount - 1);
		m_nodeSpacing.y = (m_computeDomain.yrange.second - m_computeDomain.yrange.first) / (m_nodeCount.ycount - 1);
		
	};

	void HeatEquationSolver::run(void) {
		createEquationMatrix();
		MatrixSolver solver(m_eqSparseMat,m_eqKnownVec,m_W,1e-4,1*1e6,25, m_nodeCount.xcount - 2);
		solver.solveWithGaussSeidel();
		try {
			printFileGnuplot(solver.m_x);
		}
		catch (...) {
			std::cerr << "File opening failer";
		}
	}
	HeatEquationSolver::~HeatEquationSolver() {

	}

	void HeatEquationSolver::createEquationMatrix() {
		// 5 elements per row
		for (size_t i = 0; i < m_eqMatDims; i++)
			m_eqKnownVec.push_back(0);
		std::vector<Eigen::Triplet<double> > Vec;
		Vec.reserve(m_eqMatDims * 5);
		size_t N = m_nodeCount.ycount - 2;
		size_t M = m_nodeCount.xcount - 2;
		double L = m_computeDomain.xrange.second - m_computeDomain.xrange.first;
		double H = m_computeDomain.yrange.second - m_computeDomain.yrange.first;
		double gamma = L * (N + 1) / (H * (M + 1));
		//for (size_t i = 0; i < M * N; i++)
		//	for (size_t j = 0; j < M * N; j++)
		//		//m_eqMat.A[i][j] = 0;
		for (size_t i = 0; i < M; i++) {
			for (size_t j = 0; j < N; j++) {
				// for node i,j:
				if (j == 0) {
					//m_eqMat.b[j * M + i] -= 0 * gamma;
					//m_eqMat.A[j * M + i][(j + 1) * M + i] = 1 * gamma;
					m_eqKnownVec[j * M + i] -= 0 * gamma;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, (j + 1) * M + i, 1 * gamma));
				}
				else if (j == N - 1) {
					//m_eqMat.b[j * M + i] -= 1 * gamma;
					//m_eqMat.A[j * M + i][(j - 1) * M + i] = 1 * gamma;
					m_eqKnownVec[j * M + i] -= 1 * gamma;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, (j - 1) * M + i, 1 * gamma));
				}
				else {
					//m_eqMat.A[j * M + i][(j - 1) * M + i] = 1 * gamma;
					//m_eqMat.A[j * M + i][(j + 1) * M + i] = 1 * gamma;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, (j - 1) * M + i, 1 * gamma));
					Vec.push_back(Eigen::Triplet<double>(j * M + i, (j + 1) * M + i, 1 * gamma));

				}
				if (i == 0) {
					//m_eqMat.b[j * M + i] -= 0;
					//m_eqMat.A[j * M + i][j * M + i + 1] = 1;
					m_eqKnownVec[j * M + i] -= 0;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, j * M + i + 1, 1));
				}
				else if (i == M - 1) {
					//m_eqMat.b[j * M + i] -= 0;
					//m_eqMat.A[j * M + i][j * M + i - 1] = 1;
					m_eqKnownVec[j * M + i] -= 0;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, j * M + i - 1, 1));
				}
				else {
					//m_eqMat.A[j * M + i][j * M + i - 1] = 1;
					//m_eqMat.A[j * M + i][j * M + i + 1] = 1;
					Vec.push_back(Eigen::Triplet<double>(j * M + i, j * M + i - 1, 1));
					Vec.push_back(Eigen::Triplet<double>(j * M + i, j * M + i + 1, 1));
				}
				//m_eqMat.A[j * M + i][j * M + i] = -4;
				Vec.push_back(Eigen::Triplet<double>(j * M + i, j * M + i, -4));
			}// for j
		}// for i
		m_eqSparseMat.resize(N * M, N * M);
		m_eqSparseMat.setFromTriplets(Vec.begin(), Vec.end());

	} /* createEquationMatrix */

	void HeatEquationSolver::printFileGnuplot(std::vector<double> xVec) {
		size_t N = m_nodeCount.ycount - 2;
		size_t M = m_nodeCount.xcount - 2;
		std::ofstream outputFile;
		outputFile.exceptions(std::ofstream::badbit | std::ofstream::failbit);
		outputFile.open("data.txt", std::ios_base::out);
		for (size_t i = 0; i < m_nodeCount.xcount; i++) {
			for (size_t j = 0; j < m_nodeCount.ycount; j++) {
				double x = m_computeDomain.xrange.first + i * m_nodeSpacing.x;
				double y = m_computeDomain.xrange.first + j * m_nodeSpacing.y;
				outputFile << x << " " << y << " ";
				if (j == m_nodeCount.ycount - 1) {
					outputFile << 1 << "\n";
				}
				else if (i == 0 || j == 0 || i == m_nodeCount.xcount - 1) {
					outputFile << 0 << "\n";
				}
				else {
					outputFile << xVec[(j - 1) * M + (i - 1)] << "\n";
				}
			}
			outputFile << "\n";
		}
		outputFile.close();
	}

} /* namespace ArshiaCFD */


