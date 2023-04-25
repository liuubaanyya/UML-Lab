
#include "Gaus/Gaus.h"
#include "newton/newton.hpp"

#include <iostream>
using namespace std;

inline void print(std::string s, bool skip = true)
{
    cout << s;
    if (skip)
        cout << endl;
}
template <typename T>
inline void read(T &sth)
{
    cin >> sth;
}

void senMenu()

{
    print("What you wanna do?");

    print("Gausian inversion: 1");
    print("Strassen multiplication: 2");
    print("Newton inversion: 3");
    print("Quit: 4");
}

void sendMatrixOptions()
{
    print("Choose option to enter the matrix:");
    print("Enter from console: 1");
    print("Read from file: 2");
    print("Randomly generated: 3");
}
realValued readConsoleMatrix()
{
    size_t w, h;
    print("width: ", false);
    cin >> w;

    print("height: ", false);
    cin >> h;

    vector<vector<double>> v(h, vector<double>(w));
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
            read(v[i][j]);
    }

    return realValued(v);
}

realValued readFileMatrix()
{
    string filename;
    print("filename: ", false);
    read(filename);

    freopen(filename.c_str(), "r", stdin);

    size_t w, h;

    cin >> w;

    cin >> h;

    vector<vector<double>> v(h, vector<double>(w));
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
            read(v[i][j]);
    }

    return realValued(v);
}

realValued generateRandom()
{
    size_t w, h;
    print("width: ", false);
    cin >> w;

    print("height: ", false);
    cin >> h;

    return realValued(w, h);
}

realValued giveUserMatrix()
{
    sendMatrixOptions();
    int option;
    read(option);

    switch (option)
    {
    case 1:
        return readConsoleMatrix();
        break;
    case 2:
        return readFileMatrix();

    case 3:
        return generateRandom();
    default:
        print("Not a correct option");
        return generateRandom();
        break;
    }
    return generateRandom();
}

int main()
{
    senMenu();

    int option;
    read(option);
    realValued a, b, res;
    NewtonINversion invertor;

    switch (option)
    {
    case 1:
        a = giveUserMatrix();
        a.print_matrix();
        res = GausInv(a);
        res.print_matrix();
        break;

    case 2:
        a = giveUserMatrix();
        a.print_matrix();
        res = invertor.inverse(a, 0.1);
        res.print_matrix();
        break;
    case 3:
        a = giveUserMatrix();
        b = giveUserMatrix();
        a.print_matrix();
        b.print_matrix();

        res = a.Strassen_multiplication(b);
        res.print_matrix();

        break;
    case 4:
        return 0;
        break;

    default:
        senMenu();
    }
}