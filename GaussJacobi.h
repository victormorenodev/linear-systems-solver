#ifndef GAUSS_JACOBI_H
#define GAUSS_JACOBI_H

#include "MetodoIterativo.h"

class GaussJacobi : public MetodoIterativo {
public:
    // Construtor herdando da classe base
    GaussJacobi(int maxIter, double epsilon);

    // Implementação do método solve()
    vector<double> solve(const vector<vector<double>>& A, const vector<double>& b, int n) override;
};

#endif
