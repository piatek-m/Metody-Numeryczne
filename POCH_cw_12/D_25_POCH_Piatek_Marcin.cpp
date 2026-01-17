// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <numbers>

// Wartości pochodnych
struct deriv_values
{
    double first_degree;  // Wartość pochodnej 1. stopnia
    double second_degree; // Wartość pochodnej 2. stopnia
    double third_degree;  // Wartość pochodnej 3. stopnia
};

double f(double x, int i, double h)
{
    x = x + (i * h);
    return exp(-0.5 * x) * sin(5 * std::numbers::pi * x);
}

// Funkcja obliczająca pochodne metodą różnic skończonych w przód
deriv_values RSP(double x = 1, double h = 1e-4, int i = 0)
{
    // Funkcja lambda jako wrapper funkcji f(x)
    auto y{
        [x, h](int i)
        // [x, h] -> przechwyt zmiennych x, h do funkcji
        // (int i) -> zdefiniowanie parametru funkcji lambda
        {
            return f(x, i, h);
        }};

    double deriv_1_deg{}; // Pochodna 1. stopnia
    double deriv_2_deg{}; // Pochodna 2. stopnia
    double deriv_3_deg{}; // Pochodna 3. stopnia

    deriv_1_deg = (-1.0 * y(i + 2) +
                   4.0 * y(i + 1) -
                   3.0 * y(i)) /
                  (2.0 * h);

    deriv_2_deg = (-1.0 * y(i + 3) +
                   4.0 * y(i + 2) -
                   5.0 * y(i + 1) +
                   2.0 * y(i)) /
                  (h * h);

    deriv_3_deg = (-3.0 * y(i + 4) +
                   14.0 * y(i + 3) -
                   24.0 * y(i + 2) +
                   18.0 * y(i + 1) -
                   5.0 * y(i)) /
                  (2.0 * h * h * h);

    deriv_values wynik = deriv_values(deriv_1_deg, deriv_2_deg, deriv_3_deg);

    return wynik;
}

// Funkcja obliczająca pochodne metodą różnic skończonych w tył
deriv_values RST(double x = 1, double h = 1e-4, int i = 0)
{
    // Funkcja lambda jako wrapper funkcji f(x)
    auto y{
        [x, h](int i)
        // [x, h] -> przechwyt zmiennych x, h do funkcji
        // (int i) -> zdefiniowanie parametru funkcji lambda
        {
            return f(x, i, h);
        }};

    double deriv_1_deg{}; // Pochodna 1. stopnia
    double deriv_2_deg{}; // Pochodna 2. stopnia
    double deriv_3_deg{}; // Pochodna 3. stopnia

    deriv_1_deg = (3.0 * y(i) -
                   4.0 * y(i - 1) +
                   1.0 * y(i - 2)) /
                  (2.0 * h);

    deriv_2_deg = (2.0 * y(i) -
                   5.0 * y(i - 1) +
                   4.0 * y(i - 2) -
                   1.0 * y(i - 3)) /
                  (h * h);

    deriv_3_deg = (5.0 * y(i) -
                   18.0 * y(i - 1) +
                   24.0 * y(i - 2) -
                   14.0 * y(i - 3) +
                   3.0 * y(i - 4)) /
                  (2.0 * h * h * h);

    deriv_values wynik = deriv_values(deriv_1_deg, deriv_2_deg, deriv_3_deg);

    return wynik;
}

// Funkcja obliczająca pochodne metodą różnic skończonych centralnych
deriv_values RSC(double x = 1, double h = 1e-4, int i = 0)
{
    // Funkcja lambda jako wrapper funkcji f(x)
    auto y{
        [x, h](int i)
        // [x, h] -> przechwyt zmiennych x, h do funkcji
        // (int i) -> zdefiniowanie parametru funkcji lambda
        {
            return f(x, i, h);
        }};

    double deriv_1_deg{}; // Pochodna 1. stopnia
    double deriv_2_deg{}; // Pochodna 2. stopnia
    double deriv_3_deg{}; // Pochodna 3. stopnia

    deriv_1_deg = (1.0 * y(i + 1) -
                   1.0 * y(i - 1)) /
                  (2.0 * h);

    deriv_2_deg = (1.0 * y(i + 1) -
                   2.0 * y(i) +
                   1.0 * y(i - 1)) /
                  (h * h);

    deriv_3_deg = (1.0 * y(i + 2) -
                   2.0 * y(i + 1) +
                   2.0 * y(i - 1) -
                   1.0 * y(i - 2)) /
                  (2.0 * h * h * h);

    deriv_values wynik = deriv_values(deriv_1_deg, deriv_2_deg, deriv_3_deg);

    return wynik;
}

int main()
{
    constexpr double x{1.};
    constexpr int wezel{0};
    constexpr double krok_h{1e-4};

    std::cout.precision(15);

    std::cout << "\n\n\t Krok h = ";
    std::cout << krok_h << "\n";

    // Wypisanie wyników metody RSP
    std::cout << "\n Metoda różnic skończonych w przód (RSP): \n";

    deriv_values d_vals = RSP(x, krok_h, wezel);

    // Wypisywanie wyników różniczkowania
    auto printResult{
        [x, &d_vals]()
        {
            std::cout << "\n\t\tf'(" << x << ") = " << d_vals.first_degree;
            std::cout << "\n\t\tf''(" << x << ") = " << d_vals.second_degree;
            std::cout << "\n\t\tf'''(" << x << ") = " << d_vals.third_degree;
            std::cout << "\n\n";
        }};

    printResult();

    // Wypisanie wyników metody RST
    std::cout << "\n Metoda różnic skończonych w tył (RST): \n";
    d_vals = RST();
    printResult();

    // Wypisanie wyników metody RSC
    std::cout << "\n Metoda różnic skończonych centralnych (RSC): \n";
    d_vals = RSC();
    printResult();

    std::cout << "\n\n";
    return 0;
}