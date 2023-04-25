

#include <iostream>
#include <cfloat>
#include "newton.hpp"

using namespace std;

/**
 * @brief Determines if a given real-valued matrix is square (i.e. has the same number of rows and columns).
 *
 * @param m: The matrix to be checked
 *
 * @return: True if the matrix is square, false otherwise
 */

bool NewtonINversion::isSquare(const realValued &m)
{
    size_t n = m.Matrix.size();
    for (size_t i = 0; i < n; i++)
    {
        if (m.Matrix[i].size() != n)
            return false;
    }
    return true;
}

/**
 * @brief Calculates the inverse of a given real-valued matrix using the iterative method of finding the inverse.
 *
 * The method involves multiplying the input matrix by a scaling factor, then repeatedly using the Strassen
 * algorithm to update the scaling factor until a certain threshold of accuracy is reached.
 *
 * @param A: The input matrix to be inverted
 * @param epsilon: The desired threshold of accuracy for the inverse calculation
 *
 * @return: The inverse of the input matrix, or an empty matrix if the input matrix is invalid
 */

realValued NewtonINversion::inverse(const realValued A, double epsilon)
{
    if (!isValid(A))
    {
        return {};
    }
    double t = 1.0 / (getMaxRowSum(A) * getMaxColumnSum(A));
    realValued B = multiplyMatrixByNumber(A, t);
    realValued I = getUnitMatrix(A.Matrix.size());
    realValued E;
    do
    {
        E = I.subtraction(B.Strassen_multiplication(A));
        B = I.addition(E).Strassen_multiplication(B);
    } while (getAverageSum(E) > epsilon);
    return B;
}

/**
 * @brief Determines if a given real-valued matrix is valid for matrix operations.
 *
 * A matrix is considered valid if it is non-empty and square (i.e. has the same number of rows and columns).
 *
 * @param A: The matrix to be checked
 *
 * @return: True if the matrix is valid, false otherwise
 */

bool NewtonINversion::isValid(const realValued &A)
{
    if (A.Matrix.empty())
    {
        cout << "realValued A has no elements.\n";
        return false;
    }
    else if (!isSquare(A))
    {
        cout << "realValued A is not a square realValued.\n";
        return false;
    }
    return true;
}

/**
 * @brief Multiplies a given real-valued matrix by a scalar value.
 *
 * Each element in the matrix is multiplied by the given scalar value.
 *
 * @param A: The matrix to be multiplied
 * @param val: The scalar value to multiply by
 *
 * @return: The result of the matrix multiplication as a new real-valued matrix
 */

realValued NewtonINversion::multiplyMatrixByNumber(const realValued &A, double val)
{
    realValued result = realValued(A);
    for (size_t i = 0; i < A.Matrix.size(); i++)
        for (size_t j = 0; j < A.Matrix.at(0).size(); j++)
        {
            result.Matrix[i][j] = A.Matrix[i][j] * val;
        }
    return result;
}
/**
 * @brief Computes the average of the sum of all elements in a given matrix.
 *
 * This function takes a real-valued matrix as input and returns the average value of the sum of all its elements.
 *
 * @param E: The matrix to compute the average sum of
 *
 * @return: The average sum of all elements in the matrix
 */

double NewtonINversion::getAverageSum(const realValued &E)
{
    int numberOfElements = E.Matrix.size() * E.Matrix[0].size();
    double number = 0;
    for (int i = 0; i < E.Matrix.size(); ++i)
    {
        for (int j = 0; j < E.Matrix[0].size(); ++j)
        {
            number += E.Matrix.at(i).at(j);
        }
    }
    return number / numberOfElements;
}
/**
 *
 *   @brief Generates a square unit matrix of the specified size.
 *
 *   This function takes an integer n as input and returns a square unit matrix of size n x n.
 *   A unit matrix is a matrix with diagonal elements as 1 and all other elements as 0.
 *
 *   @param n: The size of the square unit matrix to be generated.
 *   @return: A square unit matrix of size n x n.
 *
 */

realValued NewtonINversion::getUnitMatrix(int n)
{
    realValued I = realValued(vector<vector<double>>(n, vector<double>(n, double())));
    for (size_t i = 0; i < n; i++)
    {
        I.Matrix[i][i] = 1;
    }
    return I;
}
/**
 *
 *    @brief Computes the maximum row sum of a given real-valued matrix.
 *    This function takes a real-valued matrix as input and returns the maximum sum of its rows.
 *    @param A: The matrix to compute the maximum row sum of
 *
 *    @return: The maximum row sum of the matrix
 */
double NewtonINversion::getMaxRowSum(const realValued &A)
{
    auto max = DBL_MIN;
    for (const auto &i : A.Matrix)
    {
        double sum = 0;
        for (size_t j = 0; j < i.size(); j++)
        {
            sum += i[j];
        }
        if (max < sum)
        {
            max = sum;
        }
    }
    return max;
}

/**

 *    @brief Computes the maximum sum of all columns in a given real-valued matrix.
 *   This function takes a real-valued matrix as input and returns the maximum sum of all columns in the matrix.
 *   @param A: The matrix to compute the maximum sum of all columns
 *   @return: The maximum sum of all columns in the matrix

*/

double NewtonINversion::getMaxColumnSum(const realValued &A)
{
    auto max = DBL_MIN;
    for (size_t j = 0; j < A.Matrix.at(0).size(); j++)
    {
        double sum = 0;
        for (const auto &i : A.Matrix)
        {
            sum += i[j];
        }
        if (max < sum)
        {
            max = sum;
        }
    }
    return max;
}
