#include "MetodoIterativo.h"

MetodoIterativo::MetodoIterativo(int maxIter, double epsilon) : 
	maxIter(maxIter), epsilon(epsilon) {}

bool MetodoIterativo::lineCriterion(const vector<vector<double>>& A, int n)
{
	double soma;
	for (int i = 0; i < n; ++i) {
		soma = 0;
		for (int j = 0; j < n; ++j) {
			if (j != i) {
				soma += A[i][j];
			}
		}
		if (soma >= A[i][i]) {
			return false;
		}
	}
	return true;
}

bool MetodoIterativo::sassenfeldCriterion(const vector<vector<double>>& A, int n)
{
	double soma;
	vector<double> multiplicadores(n, 1);
	for (int i = 0; i < n; ++i) {
		soma = 0;
		for (int j = 0; j < n; ++j) {
			if (j != i) {
				soma += abs(A[i][j])*multiplicadores[j];
			}
		}
		soma = soma/abs(A[i][i]);
		if (soma >= 1) {
			return false;
		}
		multiplicadores[i] = soma;
	}
	return true;
}

bool MetodoIterativo::certainlyConverges(const vector<vector<double>>& A, int n)
{
	return (lineCriterion(A, n) || sassenfeldCriterion(A, n));
}

double MetodoIterativo::drk(const vector<double>& x, const vector<double>& x_ant, int n) 
{
	double dk_max = 0;
	double x_max = 0;
	
	// Calculando |d_k| e max{|x_i^(k)|} para encontrar a norma de diferença
	for (int i = 0; i < n; ++i) {
		double dk = abs(x[i] - x_ant[i]); // |d_k|
		dk_max = max(dk_max, dk); // A maior diferença
		x_max = max(x_max, abs(x[i])); // A maior magnitude de |x_i^(k)|
	}

	// Retorna o valor de drk
	return dk_max / x_max;
}

bool MetodoIterativo::canStop(double drk, int iterAtual) {
	if (iterAtual >= maxIter || drk <= epsilon) {
		return true;
	} 
	return false;
}
