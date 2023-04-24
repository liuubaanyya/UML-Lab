#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include "../basic_functions.h"
#include "Strassen.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

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
