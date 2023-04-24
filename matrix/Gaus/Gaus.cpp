#include "Gaus.h"
#include <cmath>
#include <chrono>

/**
 * @brief Calculates the inverse of a given real-valued square matrix using Gaussian elimination.
 *
 * @param A: The input matrix to be inverted
 *
 * @return: The inverse of the input matrix
 *
 * @throws: std::invalid_argument if the input matrix is not square or if it is singular (has a determinant of 0)
 */

realValued GausInv(realValued A)
{
    if (A.Matrix.size() != A.Matrix[0].size())
        throw std::invalid_argument("this matrix is not square");

    realValued inversed(A.Matrix.size(), A.Matrix.size());

    realValued src(A);

    for (int i = 0; i < A.Matrix.size(); i++)
        inversed.Matrix[i][i] = 1;

    for (int k = 0; k < A.Matrix.size(); k++)
    {
        double diag = src.Matrix[k][k];
        if (!abs(diag))
            throw std::invalid_argument("Inverse matrix doesn't exist");

        for (int i = 0; i < A.Matrix.size(); i++)
        {
            inversed.Matrix[k][i] /= diag;
            src.Matrix[k][i] /= diag;
        }

        src.Matrix[k][k] = 1;
        for (int i = 0; i < A.Matrix.size(); i++)
        {
            if (k == i)
                continue;
            else if (abs(src.Matrix[i][k]) == 0)
                continue;
            else
            {
                double div = src.Matrix[i][k] / src.Matrix[k][k];
                for (int j = 0; j < A.Matrix.size(); j++)
                {
                    src.Matrix[i][j] -= div * src.Matrix[k][j];
                    inversed.Matrix[i][j] -= div * inversed.Matrix[k][j];
                }
            }
        }
    }

    return inversed;
}

void gausTest(int beg, int end, int step)
{
    using namespace std::chrono;

    for (int i = beg; i < end; i += step)
    {

        realValued init(i, i, true);

        auto start = high_resolution_clock::now();
        realValued inv = GausInv(init);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        realValued res = init.multiplication(inv);

        double error = 0;
        for (int y = 0; y < i; ++y)
        {
            for (int x = 0; x < i; ++x)
            {
                // cout << y << ' ' << x << endl;
                if (x != y)
                    error += fabsf64(res.Matrix[y][x]);
                else
                    error += fabsf64(1 - res.Matrix[y][x]);
            }
        }

        error /= (i * i);

        cout << "Size: " << i << 'x' << i << "\n  time: " << ((double)duration.count()) / 1000 << " sec \n  error: " << error << '\n'
             << endl;
    }
}