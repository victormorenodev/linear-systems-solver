#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// função para calcular a inversa de uma matriz A
vector<vector<double>> inverseMatrix(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> invA(n, vector<double>(n, 0.0));
    vector<vector<double>> identity(n, vector<double>(n, 0.0));

    // criar matriz identidade
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }

    // calcular a inversa resolvendo n sistemas lineares (A * coluna_i = identidade_i)
    for (int i = 0; i < n; ++i) {
        // Resolver A * x = identidade_i para cada coluna i
        vector<double> x(n, 0.0);
        vector<double> b = identity[i];

        // método de eliminação de Gauss para resolver o sistema
        for (int k = 0; k < n; ++k) {
            // pivoteamento parcial
            int maxRow = k;
            for (int m = k + 1; m < n; ++m) {
                if (abs(A[m][k]) > abs(A[maxRow][k])) {
                    maxRow = m;
                }
            }
            swap(A[k], A[maxRow]);
            swap(b[k], b[maxRow]);

            // eliminação
            for (int m = k + 1; m < n; ++m) {
                double factor = A[m][k] / A[k][k];
                for (int j = k; j < n; ++j) {
                    A[m][j] -= factor * A[k][j];
                }
                b[m] -= factor * b[k];
            }
        }

        // substituição regressiva
        for (int k = n - 1; k >= 0; --k) {
            x[k] = b[k];
            for (int m = k + 1; m < n; ++m) {
                x[k] -= A[k][m] * x[m];
            }
            x[k] /= A[k][k];
        }

        // Aarmazenar a coluna i da inversa
        for (int j = 0; j < n; ++j) {
            invA[j][i] = x[j];
        }
    }

    return invA;
}

// função para multiplicar matriz por vetor
vector<double> matrixVectorMultiply(const vector<vector<double>>& A, const vector<double>& x) {
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

// função para verificar se algum deslocamento excede 0,4 cm
void checkDisplacements(const vector<double>& d) {
    for (double di : d) {
        if (abs(di) > 0.4) {
            cout << "Atenção: Deslocamento " << di << " excede 0,4 cm. Há risco de danos!" << endl;
            return;
        }
    }
    cout << "Todos os deslocamentos estão dentro do limite seguro." << endl;
}

int main() {
    // matriz A e vetor b fornecidos
    vector<vector<double>> A = {{5, 3, 1}, {5, 6, 1}, {1, 6, 7}};
    vector<double> b = {1, 2, 3};

    // calcular a inversa de A
    vector<vector<double>> invA = inverseMatrix(A);

    cout << "Inversa de A:" << endl;
    for (const auto& row : invA) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    // calcular d usando a inversa de A
    vector<double> d = matrixVectorMultiply(invA, b);

    cout << "Vetor de deslocamentos d:" << endl;
    for (double di : d) {
        cout << di << " ";
    }
    cout << endl;

    // verificar se algum deslocamento excede 0,4 cm
    checkDisplacements(d);

    return 0;
}