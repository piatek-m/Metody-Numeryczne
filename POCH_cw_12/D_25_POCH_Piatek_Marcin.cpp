// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <utility>
#include <iomanip>
#include <numbers>

struct struct_RSP
{
    double first_degree;
    double second_degree;
    double third_degree;
};

double f(double x, int i, double h)
{
    x = x + (i * h);
    return exp(-0.5 * x) * sin(5 * std::numbers::pi * x);
}

// Funkcja obliczająca pochodne metodą różnic skończonych w przód
struct_RSP RSP(double x = 1, double h = 1e-4, double eps = 1e-4)
{
    int i{0};

    double deriv_1_deg{}; // Pochodna 1. stopnia
    double deriv_2_deg{}; // Pochodna 2. stopnia
    double deriv_3_deg{}; // Pochodna 3. stopnia

    deriv_1_deg = (-1.0 * f(x, i + 2, h) + 4.0 * f(x, i + 1, h) - 3.0 * f(x, i, h)) / (2.0 * h);

    deriv_2_deg = (-1.0 * f(x, i + 3, h) + 4.0 * f(x, i + 2, h) - 5.0 * f(x, i + 1, h) + 2.0 * f(x, i, h)) / (h * h);

    deriv_3_deg = (-3.0 * f(x, i + 4, h) +
                   14.0 * f(x, i + 3, h) -
                   24.0 * f(x, i + 2, h) +
                   18.0 * f(x, i + 1, h) -
                   5.0 * f(x, i, h)) /
                  (2.0 * h * h * h);

    struct_RSP wynik = struct_RSP(deriv_1_deg, deriv_2_deg, deriv_3_deg);

    return wynik;
}

int main()
{
    constexpr double epsilon{1e-4};

    std::cout.precision(10);

    // Przedział całki
    int a{0}, b{5};

    std::cout
        << "\n\t Dokładność eps = ";
    std::cout << epsilon << "\n";

    double x = 1;

    // Wypisanie wyników metody RSP
    struct_RSP RSP_wynik = RSP();
    std::cout << "\nf'(" << x << ") = " << RSP_wynik.first_degree;
    std::cout << "\nf''(" << x << ") = " << RSP_wynik.second_degree;
    std::cout << "\nf'''(" << x << ") = " << RSP_wynik.third_degree;

    std::cout << "\n\n\n";
    return 0;
}