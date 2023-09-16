#include "MatrixSolver.hpp"
//#define DEBUG
MatrixSolver::MatrixSolver(Eigen::SparseMatrix<double> sMat, std::vector<double>b, double w, double tolerance, size_t maxIterations, double guess, int M) :
	m_sMat(sMat),
	m_b(b),
	m_w(w),
	m_tolerance(tolerance),
	m_maxIterations(maxIterations),
	m_guess(guess),
	M(M)
{
	m_x.reserve(sMat.rows());
	for (size_t i = 0; i < sMat.rows(); i++)
		m_x.push_back(m_guess);

}
bool MatrixSolver::solveWithGaussSeidel() {
	const size_t n = m_sMat.rows(); //iteration limit
	size_t k = 0; //iteration counter
	double InfNorm = std::numeric_limits<double>::max();//Just to make sure error > tolerance before the first iteration
	double AbsNorm = std::numeric_limits<double>::max();//Just to make sure error > tolerance before the first iteration
	double EucNorm = std::numeric_limits<double>::max();//Just to make sure error > tolerance before the first iteration

	double Sum = 0;

	std::vector<double> oldX(n,0); // to store values from last iteration.
	std::vector<double> old2X(n,0);
	std::vector<double> old3X(n,0);
	std::cout << "Entering the while loop";
	auto t0 = std::chrono::high_resolution_clock::now();
	while (EucNorm > m_tolerance && k < m_maxIterations)
	{
		auto t1 = std::chrono::high_resolution_clock::now();

		for (size_t i = 0; i < n; i++) {
			Sum = 0;
			if ((int)i + 1 < m_x.size()) {
				Sum += m_sMat.coeff(i, i + 1) * m_x[i + 1];
			}
			if ((int)i - 1 > 0) {
				Sum += m_sMat.coeff(i, i - 1) * m_x[i - 1];
			}
			if ((int)i + M < m_x.size()) {
				Sum += m_sMat.coeff(i, i + M) * m_x[i + M];
			}
			if ((int)i - M > 0) {
				Sum += m_sMat.coeff(i, i - M) * m_x[i - M];
			}

//#ifndef DEBUG
			//for (size_t j = 0; j < i; j++)
			//	Sum += m_sMat.coeff(i, j) * m_x[j];
			//for (size_t j = i + 1; j < n; j++)
			//	Sum += m_sMat.coeff(i, j) * m_x[j];
//#endif // !DEBUG

			
			
			//m_x = xk+1
			old3X[i] = old2X[i];//old3X= xk-2
			old2X[i] = oldX[i]; //old2X= xk-1
			oldX[i] = m_x[i]; //oldX = xk
			m_x[i] = (m_b[i] - Sum) / m_sMat.coeff(i, i);
			m_x[i] = m_w * m_x[i] + (1 - m_w) * oldX[i];

		}

		InfNorm = 0;
		AbsNorm = 0;
		EucNorm = 0;

		for (size_t i = 0; i < n; i++) {
			double deltaX = abs(oldX[i] - m_x[i]);
			AbsNorm += deltaX;
			EucNorm += deltaX*deltaX;
			InfNorm = std::max(InfNorm, deltaX);
		}
		EucNorm = pow(EucNorm, 1.0 / 2.0);


		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> loopTime = t2 - t1;
		std::chrono::duration<double, std::milli> timeTillNow = t2 - t0;
		
		statistics << k << " ";
		statistics << timeTillNow.count() << " ";
		statistics << EucNorm << "\n";

		std::cout << "it: \t" << k << "\n";
		std::cout << "Infinity norm:\t" << InfNorm << "\n";
		std::cout << "Absolute norm:\t" << InfNorm << "\n";
		std::cout << "Euclidean norm:\t" << InfNorm << "\n";
		std::cout << "loop time =\t" << loopTime.count() << "ms\n";
		k++;
		
	}
	double nu = 0;
	double dn = 0;
	for (size_t i = 0; i < n; i++) {
		//Convergnce
		nu += (m_x[i] - oldX[i]) / (oldX[i] - old2X[i]);
		dn += (oldX[i] - old2X[i]) / (old2X[i] - old3X[i]);
	}
	double q = (log10(nu) / log10(dn));
	std::ofstream outputStatistics;
	std::string Filename = "File" + std::to_string(m_w);
	Filename += ".txt";
	outputStatistics.open(Filename, std::ios_base::out);
	outputStatistics << statistics.rdbuf();
	std::cout << "omega=\t" << m_w << "\n";
	std::cout << "q=\t" << q << "\n";
	outputStatistics.close();

	if (k = m_maxIterations)
		return false;
	return true;

}