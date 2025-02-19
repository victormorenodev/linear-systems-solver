#include "MetodoIterativo.h"

MetodoIterativo::MetodoIterativo(int maxIter, double epsilon) : 
	maxIter(maxIter), epsilon(epsilon) {}

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
