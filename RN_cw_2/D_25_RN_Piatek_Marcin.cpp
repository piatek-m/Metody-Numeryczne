// Visual Studio Code, wersja 1.105.0

#include <iostream>
#include <math.h>

const double eps = 1.0e-7;

// funkcja pomocnicza do metody iteracji prostej, wylicza wart. funkcji g
double g(double x)
{
    return x = 1. / 3 * cos(x);
}

// funkcja pomocnicza do metody Newtona, wylicza wart. funkcji f
double f(double x)
{
    return x = 3 * x + (sin(x) - exp(x));
}

// funkcja pomocnicza do metody Newtona, wylicza wart. pochodnej
double df(double x)
{
    return x = 3 + (cos(x) - exp(x));
}

void iteracjaProsta()
{
    int nrIter = 0; // nr iteracji
    double x0 = 0, x1 = 0;
    do
    {
        x0 = x1;
        x1 = g(x0);

        nrIter++;
    } while ((fabs(x1 - x0)) >= eps); // fabs ~ floating absolute, funkcja zwraca wartosc bezwzgledna z wartości typu zmiennoprzecinkowego

    std::cout << "Wynik metody iteracji prostej po " << nrIter << " iteracjach to: " << x1 << ".\n\n\n\n\n";
}

void metodaNewtona()
{
    int nrIter{0};
    double x0{}, x1{0};

    do
    {
        x0 = x1;
        x1 = x0 - (f(x0) / df(x0));

        nrIter++;
    } while (fabs(x1 - x0) >= eps);

    std::cout << "Wynik metody Newtona po " << nrIter << " iteracjach to: " << x1 << ".\n\n\n\n\n";
}

int main()
{
    std::cout.precision(10);
    iteracjaProsta();

    metodaNewtona();

    return 0;
}
