#include "basic_functions.h"

realValued realValued::Strassens_algorithm(realValued matrix, size_t n)
{
    if (n <= 64)
    {
        return (this->multiplication(matrix));
    }
    n = n / 2;
    realValued a11(n), a12(n), a21(n), a22(n);
    realValued b11(n), b12(n), b21(n), b22(n);
    this->split(a11, a12, a21, a22);
    matrix.split(b11, b12, b21, b22);
    realValued p1 = (a11.addition(a22)).Strassens_algorithm(b11.addition(b22), n);
    realValued p2 = (a21.addition(a22)).Strassens_algorithm(b11, n);
    realValued p3 = a11.Strassens_algorithm(b22.subtraction(b12), n);
    realValued p4 = a22.Strassens_algorithm((b11.subtraction(b21)), n);
    realValued p5 = (a11.addition(a12)).Strassens_algorithm(b22, n);
    realValued p6 = (a11.subtraction(a21)).Strassens_algorithm(b11.addition(b12), n);
    realValued p7 = (a22.subtraction(a12)).Strassens_algorithm(b21.addition(b22), n);

    realValued c11 = (p1.addition(p4)).addition(p5.subtraction(p7));
    realValued c12 = p3.addition(p5);
    realValued c21 = p2.addition(p4);
    realValued c22 = (p2.subtraction(p1)).addition(p3.addition(p6));
    realValued result(2 * n);
    result.collect(c11, c12, c21, c22);
    return result;
}

/**
 *
 *	@brief Splits a real-valued matrix into four submatrices.
 *
 *	This method takes a real-valued matrix and splits it into four submatrices, each of size n/2 x n/2, where n is the size of the original matrix. The submatrices are stored in the passed parameters a1, a2, a3, and a4, respectively.
 *	@param a1: The first submatrix to store
 *
 *	@param a2: The second submatrix to store
 *	@param a3: The third submatrix to store
 *	@param a4: The fourth submatrix to store
 *	@return: None
 *
 */
void realValued::split(realValued &a1, realValued &a2, realValued &a3, realValued &a4)
{
    size_t n = this->Matrix.size() / 2;
    for (size_t i = 0; i < n; ++i)
    {
        a1.Matrix[i].resize(n);
        a2.Matrix[i].resize(n);
        a3.Matrix[i].resize(n);
        a4.Matrix[i].resize(n);

        copy(begin(this->Matrix[i]), begin(this->Matrix[i]) + n, begin(a1.Matrix[i]));
        copy(begin(this->Matrix[i]) + n, end(this->Matrix[i]), begin(a2.Matrix[i]));
        copy(begin(this->Matrix[i + n]), begin(this->Matrix[i + n]) + n, begin(a3.Matrix[i]));
        copy(begin(this->Matrix[i + n]) + n, end(this->Matrix[i + n]), begin(a4.Matrix[i]));
    }
}
/**
 *
 *	@brief Collects 4 submatrices to form a larger matrix.
 *	This method takes four realValued matrices (c11, c12, c21, c22) as input and combines them to form a larger matrix,
 *	with c11 in the top-left quadrant, c12 in the top-right quadrant, c21 in the bottom-left quadrant, and c22 in the bottom-right quadrant.
 *
 *	@param c11: The top-left quadrant matrix to be combined.
 *	@param c12: The top-right quadrant matrix to be combined.
 *	@param c21: The bottom-left quadrant matrix to be combined.
 *	@param c22: The bottom-right quadrant matrix to be combined.
 *	@return: void, as the method modifies the calling realValued object.
 *
 * */

void realValued::collect(realValued &c11, realValued &c12, realValued &c21, realValued &c22)
{
    size_t n = this->Matrix.size() / 2;
    for (size_t i = 0; i < n; ++i)
    {
        copy(begin(c11.Matrix[i]), end(c11.Matrix[i]), begin(this->Matrix[i]));
        copy(begin(c12.Matrix[i]), end(c12.Matrix[i]), begin(this->Matrix[i]) + n);
        copy(begin(c21.Matrix[i]), end(c21.Matrix[i]), begin(this->Matrix[i + n]));
        copy(begin(c22.Matrix[i]), end(c22.Matrix[i]), begin(this->Matrix[i + n]) + n);
    }
}
/**
 *
 *	@brief Computes the next power of two greater than or equal to a given number.
 *	This function takes a positive integer as input and returns the next power of two greater than or equal to it.
 *	@param number: The number to compute the next power of two of
 *	@return: The next power of two greater than or equal to the given number
 */
size_t degrees_two(const size_t &number)
{
    size_t d_two = 2;
    while (number > d_two)
    {
        d_two = d_two << 1;
    }
    return d_two;
}
/**
 *
 *	@brief Implements Strassen's algorithm for multiplying two real-valued matrices.
 *	This function uses Strassen's algorithm to compute the product of two real-valued matrices of size n x n. If n is smaller than or equal to 64, the function uses standard matrix multiplication.
 *	@param matrix: The real-valued matrix to be multiplied with the current object
 *	@param n: The size of the matrices
 *	@return: The product of the two matrices
 */

/**
 *
 *	@brief Performs matrix multiplication using the Strassen algorithm.
 *	This function performs matrix multiplication using the Strassen algorithm, which is an algorithm for computing the product of two matrices. If the size of matrices is less than or equal to 64, it performs the regular multiplication. If the size is greater than 64, it expands the matrices to the nearest power of 2, and then applies the Strassen algorithm.
 *
 *	@param matrix: The matrix to be multiplied with the calling object
 *
 *	@return: The resulting matrix obtained after multiplication using Strassen algorithm
 *	@throws runtime_error: If the matrices are not of correct size for multiplication
 *	*/

realValued realValued::Strassen_multiplication(realValued matrix)
{
    if (!(this->correct_size_for_multiplication(matrix)))
    {
        throw runtime_error("");
    }
    if (this->rows <= 64)
    {
        return (this->multiplication(matrix));
    }
    size_t first_rows = this->rows, first_columns = this->columns, second_rows = matrix.rows, second_columns = matrix.columns;
    size_t new_size = max(max(first_rows, first_columns), second_columns); // first_columns=second_rows
    new_size = degrees_two(new_size);

    if (this->rows != new_size || this->columns != new_size)
    {
        this->expand(new_size);
        matrix.expand(new_size);
    }
    realValued result = this->Strassens_algorithm(matrix, new_size);
    if (first_rows != new_size)
    {
        result.Matrix.resize(first_rows);
    }
    if (second_columns != new_size)
    {
        for (auto &r : result.Matrix)
        {
            r.resize(second_columns);
        }
    }
    return result;
}

auto realValued::test_time(const realValued &matrix)
{
    auto start = high_resolution_clock::now();
    realValued result = this->multiplication(matrix);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    long int result_time = duration.count();
    return result_time;
}

auto realValued::test_time_Strassen(const realValued &matrix)
{
    auto start = high_resolution_clock::now();
    realValued result = this->Strassen_multiplication(matrix);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    long int result_time = duration.count();
    return result_time;
}

void realValued::test(const realValued &matrix)
{

    realValued result1 = this->multiplication(matrix);
    long int result_time = this->test_time(matrix);
    cout << "Execution time of the usual algorithm: ";
    cout << result_time / 1000000 << " seconds, "
         << (result_time % 1000000) / 1000 << " milliseconds, "
         << result_time % 1000 << " microseconds";
    cout << endl;

    realValued result2 = this->Strassen_multiplication(matrix);
    long int result_time_str = this->test_time_Strassen(matrix);
    cout << "Execution time of the Strassen algorithm: ";
    cout << result_time_str / 1000000 << " seconds, "
         << (result_time_str % 1000000) / 1000 << " milliseconds, "
         << result_time_str % 1000 << " microseconds";
    cout << endl;
    cout << endl;
    assert(result1.compare(result2));
    /*result1.print_matrix();
    cout << endl;
    result2.print_matrix();*/
}

void examples_tests()
{
    realValued a1(15, 44);
    realValued b1(44, 15);
    a1.test(b1);

    realValued a2({{2, 1, 1, 1, 1, 1, 1},
                   {2, 3, 4.5, 4.5, 1, 1, 1},
                   {2, 3, 1, 12.4, 1, 1, 1},
                   {2, 2.3, 5, 3.4, 1, 1, 1},
                   {1, 1, 1, 1, 1, 1, 1}});
    realValued b2({{1, 1},
                   {2, 4},
                   {2, 2},
                   {3.5, 3.4},
                   {1, 1},
                   {1, 1},
                   {1, 1}});
    a2.test(b2);

    realValued a3(1, 1);
    realValued b3(1, 1);
    a3.test(b3);

    realValued a4({{0, 43434, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {2323, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {2323, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});
    realValued b4({{0, 43434, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 2, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 3, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 5, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 2},
                   {0, 0, 0, 0, 0, 0},
                   {0, 0, 0, 0, 0, 0}});
    a4.test(b4);
}