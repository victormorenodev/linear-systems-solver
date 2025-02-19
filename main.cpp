#include <iostream>
#include "GaussSeidel.h"

int main() {
    int n;
    cout << "Digite o tamanho da matriz A (n): ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n), x(n, 0);

    // Lendo a matriz A
    cout << "\nDigite os termos da matriz A (n x n):" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    // Lendo o vetor b
    cout << "\nDigite os termos do vetor b (n x 1), separados por espaço:" << endl;
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    // Lendo a precisão
    double epsilon;
    cout << "\nDigite a precisão (epsilon): ";
    cin >> epsilon;

    // Inicializando o vetor x^(0) com bi/aii
    for (int i = 0; i < n; ++i) {
        x[i] = b[i] / A[i][i];
    }

    // Exibindo x^(0)
    cout << "\n...Iniciando método de Gauss-Seidel para Ad = b...\n\nValor inicial de x^(0):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << (i + 1) << "^(0) = " << x[i] << endl;
    }

    GaussSeidel gaussSeidel(20, epsilon);
    // Resolvendo usando Gauss-Seidel
    vector<double> d = gaussSeidel.solve(A, b, n);

    // Resultados de gaussJacobi
    cout << "\nResultado do sistema Ad = b:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "d" << (i + 1) << " = " << d[i];

        // Alerta se o valor de d_i for > 0.4
        if (d[i] > 0.4) {
            cout << " (ALERTA: Deslocamento > 0.4)";
        }  else if ((0.4 - d[i]  > 0.99) && ( d[i] > 0 )) {
            cout << " (ATENÇÃO: Deslocamento próximo de 0.4)";
        }

        cout << endl;
    }
    
    return 0;
}
