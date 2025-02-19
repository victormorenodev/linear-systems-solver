#ifndef METODO_H
#define METODO_H

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class MetodoIterativo {
protected:
	int maxIter;
	double epsilon;

	// Função para calcular a norma de drk
	double drk(const vector<double>& x, const vector<double>& x_ant, int n);
	
	// Função que retorna se já convergiu
	bool canStop(double drk, int iterAtual);
public:
	
	// Construtor
	MetodoIterativo(int maxIter, double epsilon);

	// Destrutor
	virtual ~MetodoIterativo() = default;

	// Função que deverá ser polimorfizada dependendo da técnica
	virtual vector<double> solve(const vector<vector<double>>& A, const vector<double>& b, int n) = 0;
};

#endif
