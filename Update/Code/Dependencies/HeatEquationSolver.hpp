/*
 * HeatEquationSolver.hpp
 *
 *  Created on: Mar 26, 2023
 *      Author: ArshiaSaffari
 */
#pragma once

#include <utility> //to use std::pair
#include <vector> //to use std::vector
#include <MatrixSolver.hpp> //to solve the matrix
#include "Eigen/SparseCore" //to create sparse matrix. A dense matrix will not work because too much memory is required.

struct Index2D {
	size_t xcount;
	size_t ycount;
};
struct Range2D {
	std::pair<double, double> xrange;
	std::pair<double, double> yrange;
};
struct Vec2d
{
	double x;
	double y;
};
namespace ArshiaCFD {
	class HeatEquationSolver {
	public:
		HeatEquationSolver(Range2D ComputeDomain, Index2D NCount, double m_W);
		void createEquationMatrix();
		void run(void);
		std::vector<double> m_eqKnownVec;
		std::vector<double> m_eqX;
		void printFileGnuplot(std::vector<double>);
		~HeatEquationSolver();
		double m_W;
	protected:
		const Range2D m_computeDomain;
		const Index2D m_nodeCount;
		Vec2d m_nodeSpacing;
		Eigen::SparseMatrix<double> m_eqSparseMat;

		size_t m_eqMatDims;

	};



} /* namespace ArshiaCFD */

