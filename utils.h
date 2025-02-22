#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
using namespace std;

class utils
{
private:
public:
    static vector<double> multMatriz(vector<vector<double>>& A, vector<double>& B);
    static vector<vector<double>> identidade(int n);
    static void printMatriz(const vector<vector<double>>& matriz, const string& name);
    static void printMatrizColuna(const vector<double>& matriz, const string& name);
};
#endif