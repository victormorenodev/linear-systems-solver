#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Função para calcular a norma de drk
double drk(const vector<double>& x, const vector<double>& x_ant, int n) {
    double dk_max = 0;
    double x_max = 0;

    // Calculando |d_k| e max{|x_i^(k)|} para encontrar a norma da diferença
    for (int i = 0; i < n; ++i) {
        double dk = abs(x[i] - x_ant[i]); // |d_k|
        dk_max = max(dk_max, dk); // A maior diferença
        x_max = max(x_max, abs(x[i])); // A maior magnitude de |x_i^(k)|
    }

    // Retorna o valor de drk
    return dk_max / x_max;
}

// Função para fazer a decomposição LU de uma matriz A
void fatoracao_LU(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int n) {
    // Inicializando as matrizes L e U
    L = vector<vector<double>>(n, vector<double>(n, 0));
    U = vector<vector<double>>(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        // Decomposição para a matriz U
        for (int k = i; k < n; ++k) {
            U[i][k] = A[i][k];
            for (int j = 0; j < i; ++j) {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }

        // Decomposição para a matriz L
        for (int k = i; k < n; ++k) {
            if (i == k)
                L[i][i] = 1;  // A diagonal de L é 1
            else {
                L[k][i] = A[k][i];
                for (int j = 0; j < i; ++j) {
                    L[k][i] -= L[k][j] * U[j][i];
                }
                L[k][i] /= U[i][i];
            }
        }
    }
}

// Função para resolver Ly = b usando substituição para frente
vector<double> subst_L(const vector<vector<double>>& L, const vector<double>& b, int n) {
    vector<double> y(n, 0);
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
    return y;
}

// Função para resolver Ux = y usando substituição para trás
vector<double> subst_U(const vector<vector<double>>& U, const vector<double>& y, int n) {
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}

// Função para calcular a inversa de A usando LU Decomposition
vector<vector<double>> inversa_LU(vector<vector<double>>& A, int n) {
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    fatoracao_LU(A, L, U, n);

    vector<vector<double>> A_inv(n, vector<double>(n, 0));

    // Resolver L * y = I para cada coluna de A^-1
    for (int i = 0; i < n; ++i) {
        vector<double> e(n, 0);  // Vetor coluna da matriz identidade
        e[i] = 1;

        vector<double> y = subst_L(L, e, n);
        vector<double> x = subst_U(U, y, n);

        // A coluna de A^-1 é o vetor x
        for (int j = 0; j < n; ++j) {
            A_inv[j][i] = x[j];
        }
    }

    return A_inv;
}

// Função para resolver o sistema Ax = b usando o método de Gauss-Jacobi
vector<double> gaussJacobi(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, double epsilon, int n) {
    vector<double> x_ant(n, 0);

    int iter = 0;

    do {
        for (int i = 0; i < n; ++i) {
            x_ant[i] = x[i];
        }

        // Iteração do método Gauss-Jacobi
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
        cout << "\nIteração " << ++iter << ":\n";
        for (int i = 0; i < n; ++i) {
            double dk = abs(x[i] - x_ant[i]);
            cout << "d" << (i + 1) << " = " << dk << endl;
        }
        cout << "Norma da diferença (drk) = " << dr_k << endl;

    } while (drk(x, x_ant, n) > epsilon);  // Convergência com base em drk

    return x;
}

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
    cout << "\n...Iniciando método de Gauss-Jacobi para Ad = b...\n\nValor inicial de x^(0):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << (i + 1) << "^(0) = " << x[i] << endl;
    }

    // Resolvendo usando gaussJacobi
    vector<double> d = gaussJacobi(A, b, x, epsilon, n);

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
