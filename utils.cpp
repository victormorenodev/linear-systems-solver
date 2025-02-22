#include "utils.h"
#include <iostream>

using namespace std;
vector<double> utils::multMatriz(
    vector<vector<double>>& A, 
    vector<double>& B)
{
    int N = A.size(); //eu sei que B tem o mesmo número de colunas que A
    vector<double> resultado(N, 0); // Vetor resultado inicializado com zeros

    // Multiplicação otimizada
    for (int i = 0; i < N; ++i) {
        double soma = 0;
        for (int j = 0; j < N; ++j) {
            soma += A[i][j] * B[j];
        }
        resultado[i] = soma;
    }

    return resultado;
}

vector<vector<double>> utils::identidade(int n){
    vector<vector<double>> resultado(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++)
    {
        resultado[i][i] = 1;
    }
    
    return resultado;
}

void utils::printMatriz(const vector<vector<double>>& matriz, const string& name){
    cout << "\nMatriz " << name << ":\n";
    for (const auto& row : matriz) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

void utils::printMatrizColuna(const vector<double>& matriz, const string& name){
    cout << "\nMatriz " << name << ":\n";
    for (const auto& val : matriz) {
        cout << val << "\t";
    }
    cout << endl;
}