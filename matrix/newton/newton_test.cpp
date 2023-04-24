#pragma once

#include <iostream>
#include <chrono>
#include "newton.cpp"

using namespace std;
using namespace chrono;

void displayMatrix(const realValued &matrix)
{
    cout << "\n";
    for (size_t i = 0; i < matrix.Matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix.Matrix.at(i).size(); j++)
        {
            cout << " " << left << matrix.Matrix[i][j];
        }
        cout << "\n\n";
    }
}

int main()
{
    double epsilon = 4;
    NewtonINversion inv;
    cout
        << " *** Matrix tests ***\n\n";
    cout << " Test #1\n Epsilon = " << epsilon << "\n\n";
    vector<vector<double>> A{
        vector<double>{1, 4, 0},
        vector<double>{7, 1, 7},
        vector<double>{4, 1, 1}};
    cout << " Matrix A:\n";
    realValued mtrx(A);
    auto start = high_resolution_clock::now();
    realValued matrix = inv.inverse(mtrx, epsilon);
    auto stop = high_resolution_clock::now();
    cout << " Inversed A:\n";
    displayMatrix(matrix);
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Execution time: " << duration.count() << endl;
    cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n";
}
