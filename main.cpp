#include <iostream>
#include <vector>
#include "GaussSeidel.h"
#include "GaussJacobi.h"
#include "MetodoIterativo.h"
#include "utils.h"
#define MAX_ITER 20

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

    utils::printMatriz(A, "A");
    utils::printMatrizColuna(b, "B");

    if (!MetodoIterativo::certainlyConverges(A, n))
    {
      cout << "A matriz A não passou nos critérios de convergência, ela PODE não convergir" << endl;  
    }

    //criando a matriz identidade
    vector<vector<double>> identidade = utils::identidade(n);
    
    //criando os objetos dos métodos
    GaussJacobi gaussJacobi(MAX_ITER, epsilon);
    GaussSeidel gaussSeidel(MAX_ITER, epsilon);

    //criando as matrizes que serão as inversas
    vector<vector<double>> inversaJacobi(n, vector<double>(n));
    vector<vector<double>> inversaSeidel(n, vector<double>(n));

    cout << "CALCULANDO A INVERSA POR JACOBI E SEIDEL" << endl;

    for (int i = 0; i < n; i++)
    {
        //calculando a matriz inversa coluna por coluna
        cout << "\nCOLUNA " << i << " MÉTODO JACOBI";
        vector<double> colunaInversaJacobi = gaussJacobi.solve(A, identidade[i], n);
        cout << "\nCOLUNA " << i << " MÉTODO SEIDEL";
        vector<double> colunaInversaSeidel = gaussSeidel.solve(A, identidade[i], n);

        for (int j = 0; j < n; j++)
        {
            //preenchendo as matrizes inversas por coluna
            inversaJacobi[j][i] = colunaInversaJacobi[j];
            inversaSeidel[j][i] = colunaInversaSeidel[j];
        }
    }

    //por fim, multiplicando as inversas por b, para obter d
    vector<double> dJacobi = utils::multMatriz(inversaJacobi, b);
    vector<double> dSeidel = utils::multMatriz(inversaSeidel, b);

    utils::printMatriz(inversaJacobi, "Inversa Jacobi");
    utils::printMatriz(inversaSeidel, "Inversa Seidel");
    utils::printMatrizColuna(dJacobi, "d Jacobi");
    utils::printMatrizColuna(dSeidel, "d Seidel");

    utils::checkDeslocamentos(dJacobi, n , "Jacobi");
    utils::checkDeslocamentos(dSeidel, n , "Seidel");
}
