#pragma once

#include <vector>
#include <iostream>

using namespace std;
class realValued
{
public:
	vector<vector<double>> Matrix;
	size_t rows;	//= matrix.size()
	size_t columns; // matrix[0].size()

	bool compare(realValued matrix)
	{
		if ((this->Matrix.size() != matrix.Matrix.size()) || (this->Matrix[0].size() != matrix.Matrix[0].size()))
		{
			return false;
		}
		for (size_t i = 0; i < this->Matrix.size(); i += 1)
		{
			for (size_t j = 0; j < this->Matrix[0].size(); j += 1)
			{
				if (size_t(this->Matrix[i][j] - matrix.Matrix[i][j]) != 0)
				{
					return false;
				}
			}
		}
		return true;
	}

	void print_matrix()
	{

		for (size_t i = 0; i < this->Matrix.size(); i++)
		{
			for (size_t j = 0; j < this->Matrix[0].size(); j++)
			{
				cout << this->Matrix[i][j] << " ";
			}
			cout << endl;
		}
	}

	bool correct_size_for_multiplication(const realValued matrix)
	{
		if (this->columns != matrix.rows)
		{
			cout << this->columns << " " << matrix.rows << endl;
			return false;
		}
		return true;
	}

	bool correct_sizes(const realValued matrix)
	{
		if (this->rows != matrix.rows || this->columns != matrix.columns)
		{
			false;
		}
		return true;
	}

	realValued multiplication(const realValued matrix)
	{
		this->rows = this->Matrix.size();
		this->columns = this->Matrix[0].size();
		if (!this->correct_size_for_multiplication(matrix))
		{
			throw runtime_error("");
		}
		realValued result(this->rows);
		for (size_t i = 0; i < this->rows; ++i)
		{
			result.Matrix[i].resize(matrix.columns);
			for (size_t j = 0; j < matrix.columns; ++j)
			{
				for (size_t k = 0; k < this->columns; ++k)
				{
					result.Matrix[i][j] += this->Matrix[i][k] * matrix.Matrix[k][j];
				}
			}
		}
		return result;
	}

	realValued addition(realValued matrix)
	{
		if (!this->correct_sizes(matrix))
		{
			throw runtime_error("");
		}
		for (size_t i = 0; i < matrix.rows; ++i)
		{
			for (size_t j = 0; j < matrix.columns; ++j)
			{
				matrix.Matrix[i][j] += this->Matrix[i][j];
			}
		}
		return matrix;
	}

	realValued subtraction(realValued matrix)
	{
		if (!correct_sizes(matrix))
		{
			throw runtime_error("");
		}
		for (size_t i = 0; i < matrix.rows; ++i)
		{
			for (size_t j = 0; j < matrix.columns; ++j)
			{
				matrix.Matrix[i][j] -= this->Matrix[i][j];
			}
		}
		return matrix;
	}

	void expand(const size_t &size)
	{
		this->Matrix.resize(size);
		for (auto &row : this->Matrix)
		{
			row.resize(size, 0);
		}
	}

	// Constructors

	realValued() = default;

	realValued(const int &n)
	{
		vector<vector<double>> matrix;
		matrix.resize(n);
		for (int i = 0; i < n; i++)
		{
			matrix[i].resize(n);
		}
		Matrix = matrix;
		rows = matrix.size();
		columns = matrix[0].size();
	}

	realValued(vector<vector<double>> matrix)
	{
		Matrix = matrix;
		rows = matrix.size();
		columns = matrix[0].size();
	}

	realValued(const size_t &N, const size_t &M, bool randInit = false)
	{
		vector<vector<double>> matrix;
		rows = N;
		columns = M;
		srand(time(NULL));

		for (int i = 0; i < this->rows; ++i)
		{

			vector<double> line;
			for (int j = 0; j < this->columns; ++j)
			{
				double number = rand() % 1000 + 1;
				line.push_back((randInit) ? number : 0);
			}
			matrix.push_back(line);
		}
		Matrix = matrix;
	}

	// for Strassen multiplication
	void split(realValued &a1, realValued &a2, realValued &a3, realValued &a4);
	void collect(realValued &c11, realValued &c12, realValued &c21, realValued &c22);
	realValued Strassens_algorithm(realValued matrix, size_t n);
	realValued Strassen_multiplication(realValued matrix);
	auto test_time(const realValued &matrix);
	auto test_time_Strassen(const realValued &matrix);
	void test(const realValued &matrix);
};
