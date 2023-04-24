#include <vector>
#include "../basic_functions.h"

class NewtonINversion
{
public:
    realValued inverse(realValued A, double epsilon);

private:
    bool isValid(const realValued &A);

    realValued multiplyMatrixByNumber(const realValued &A, double val);

    double getAverageSum(const realValued &E);

    realValued getUnitMatrix(int n);

    double getMaxRowSum(const realValued &A);

    double getMaxColumnSum(const realValued &A);

    bool isSquare(const realValued &m);
};