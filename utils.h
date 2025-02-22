#ifndef UTILS_H
#define UTILS_H

#include <vector>
using namespace std;

class utils
{
private:
public:
    static vector<double> multMatriz(vector<vector<double>>& A, vector<double>& B);
    static vector<vector<double>> identidade(int n);
};
#endif