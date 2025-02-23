#include "GaussJacobi.h"

// Construtor chamando a classe base
GaussJacobi::GaussJacobi(int maxIter, double epsilon)
    : MetodoIterativo(maxIter, epsilon) {}

// Implementação do método solve() (Gauss-Jacobi)
vector<double> GaussJacobi::solve(const vector<vector<double>>& A, const vector<double>& b, int n) {

    vector<double> x_ant(n, 0);
    vector<double> x(n, 0);

    // Inicializando o vetor x^(0) com bi/aii
    for (int i = 0; i < n; ++i) {
        x[i] = b[i] / A[i][i];
    }

    int iter = 0;

    do {
        // Salva os valores anteriores antes de começar uma nova iteração
        for (int i = 0; i < n; ++i) {
            x_ant[i] = x[i];
        }

        // Iteração do método Gauss-Jacobi
        for (int i = 0; i < n; ++i) {
            double soma = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    soma += A[i][j] * x_ant[j];  // Usando os valores da iteração anterior
                }
            }
            x[i] = (b[i] - soma) / A[i][i];  // Atualiza x[i] com os valores anteriores
        }

        // Calcular drk (erro entre iterações)
        double dr_k = drk(x, x_ant, n);

        iter++;

        // Exibir valores de dk
        cout << "\nIteração " << iter << ":\n";
        for (int i = 0; i < n; ++i) {
            double dk = abs(x[i] - x_ant[i]);
            cout << "d" << (i + 1) << " = " << dk << endl;
        }
        cout << "Norma da diferença (drk) = " << dr_k << endl;

    } while (!canStop(drk(x, x_ant, n), iter));

    return x;
}
