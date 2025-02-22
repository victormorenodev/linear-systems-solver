#include "GaussJacobi.h"

// Construtor chamando a classe base
GaussJacobi::GaussJacobi(int maxIter, double epsilon)
    : MetodoIterativo(maxIter, epsilon) {}

// Implementação do método solve() (Gauss-Jacobi)
vector<double> GaussJacobi::solve(const vector<vector<double>>& A, const vector<double>& b, int n) {
    vector<double> x_ant(n, 0); // Iteração anterior
    vector<double> x(n, 0);     // Iteração atual

    cout << "Método de Gauss Jacobi" << endl;

    // Inicialização x^(0)
    for (int i = 0; i < n; ++i) {
        x_ant[i] = b[i] / A[i][i];
    }

    int iter = 0;

    do {
        // Copia x_ant para x
        for (int i = 0; i < n; ++i) {
            double soma = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    soma += A[i][j] * x_ant[j];
                }
            }
            x[i] = (b[i] - soma) / A[i][i];
        }

        // Calcular drk
        double dr_k = drk(x, x_ant, n);

        // Exibir valores de dk
        cout << "\nIteração " << iter + 1 << ":\n";
        for (int i = 0; i < n; ++i) {
            double dk = abs(x[i] - x_ant[i]);
            cout << "d" << (i + 1) << " = " << dk << endl;
        }
        cout << "Norma da diferença (drk) = " << dr_k << endl;

        // Atualiza x_ant para a próxima iteração
        x_ant = x;
        iter++;

    } while (!canStop(drk(x, x_ant, n), iter));

    return x;
}
