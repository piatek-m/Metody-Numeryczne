// Środowisko: Visual Studio Code wersja 1.105.1

#include <iostream>
#include <vector>
#include <math.h>
#include <utility>
#include <iomanip>

// Alias dla typu wektora dwuwymiarowego
using Matrix = std::vector<std::vector<double>>;
// Alias dla typu wektora typu double
using Vector = std::vector<double>;

// Funkcja do wypisywania macierzy
void printMatrix(const Matrix &mtr)
{
    std::cout << std::setprecision(2);
    std::cout << "\n";
    for (int i = 0; i < mtr.size(); i++)
    {
        std::cout << " | ";
        for (int j = 0; j < mtr.size() + 1; j++)
        {
            if (j == mtr.size())
                std::cout << std::setw(6) << " | ";
            std::cout << std::setw(6) << mtr[i][j];
        }
        std::cout << std::setw(6) << " | \n";
    }
}

// Funkcja pomocnicza do wyświetlania wyniku
void printVector(const Vector &vector)
{
    std::cout << "\n";
    for (int i = 0; i < vector.size(); i++)
    {
        std::cout << " | " << vector[i] << " | \n";
    }
}

// Funkcja rozwiązująca układ równań metodą eliminacji Gaussa-Jordana
Vector gaussJordan(Matrix matrix, Vector vector)
{
    int dimension(matrix.size()); // Wymiar macierzy i l. zmiennych wektora

    double max{};   // Największa znaleziona wartość w macierzy
    int rowOfMax{}; // Wiersz, w którym jest największa wartość

    constexpr double eps{1.0e-10}; // Dokładność epsilon

    // Sprawdzenie wymiarów

    // Czy macierz jest kwadratowa
    if (matrix[0].size() != dimension)
    {
        throw std::invalid_argument("Macierz musi być kwadratowa!");
    }

    // Czy ilość elementów wektora jest równa wymiarom macierzy (elementy wektora = wiersze macierzy)
    if (vector.size() != dimension)
    {
        throw std::invalid_argument("Wektor ma nieprawidłowy wymiar!");
    }

    // Utworzenie macierzy rozszerzonej (dodanie wektora stałych jako ostatniej kolumny)
    for (int i = 0; i < dimension; i++)
    {
        // Do każdego wiersza dopisujemy element z wektora
        matrix[i].push_back(vector[i]);
    }

    for (int k = 0; k < dimension; k++)
    {
        // Przesuwanie się po przekątnej, ustawienie max na obecny element
        max = matrix[k][k];
        rowOfMax = k;

        // Przechodzenie po elementach w dół od przekątnej (ruch w jednej kolumnie, iterowanie wierszy)
        for (int i = k; i < dimension; i++)
        {
            // Znalezienie bezwględnie największej wartości w kolumnie (od przekątnej [włącznie] w dół)
            if (abs(matrix[i][k]) > abs(max))
            {
                max = matrix[i][k];
                rowOfMax = i; // Zapamiętanie, w którym wierszu jest największa wartość
            }
        }

        // Jeżeli bezwględnie największa wartość jest zerem to mamy macierz osobliwą - przerwanie algorytmu
        if (abs(max) < eps)
        {
            throw std::runtime_error("Macierz układu osobliwa - brak jednoznacznego rozwiązania!");
        }

        // Zamiana wierszy, jeśli pivot nie jest na przekątnej
        if (rowOfMax != k)
        {
            // Zamiana całych wierszy
            std::swap(matrix[k], matrix[rowOfMax]);
        }

        // Normalizacja wiersza pivotowego
        for (int j = k; j <= dimension; j++)
        {
            // Dzielenie przez pivot
            matrix[k][j] *= (1.0 / max);
        }

        // Przechodzenie po wierszach
        for (int i = 0; i < dimension; i++)
        {

            // Jeśli nie jesteśmy na przekątnej
            if (i != k)
            {
                // Element używany do zerowania
                double p{matrix[i][k]};

                // Przechodzenie po kolumnach
                for (int j = k; j <= dimension; j++)
                {

                    // Eliminacja elementów
                    matrix[i][j] = matrix[i][j] - p * (matrix[k][j]);
                }
            }
        }
    }

    // Wyciągamy rozwiązanie (ostatnia kolumna)
    Vector solution(dimension);

    // Przechodzenie po wierszach
    for (int i = 0; i < dimension; i++)
    {
        // W każdym wierszu wyciągnięcie ostatniej wartości (ostatniej kolumny) do wektora wynikowego
        solution[i] = matrix[i][dimension];
    }

    return solution;
}

double f(double x, double y)
{
    return (pow(x, 3) * pow(y, 2)) + // x^3*y^2
           (x * y) +                 // xy
           (x * pow(y, 3)) -         // xy^3
           3;                        // -3
}
double g(double x, double y)
{
    return pow(x, 2) +               // x^2
           (pow(x, 2) * pow(y, 2)) - // x^2*y^2
           2 * x * y;                // -2xy
}
double fx(double x, double y)
{
    return (3 * pow(x, 2) * pow(y, 2)) + // 3x^2*y^2
           y +                           // y
           pow(y, 3);                    // y^3
}
double fy(double x, double y)
{
    return (2 * pow(x, 3) * y) + // 2x^3*y
           x +                   // x
           3 * x * pow(y, 2);    // 3x*y^2
}
double gx(double x, double y)
{
    return 2 * x +             // 2x
           2 * x * pow(y, 2) - // 2x*y^2
           2 * y;              // -2y
}
double gy(double x, double y)
{
    return 2 * pow(x, 2) * y - // 2x^2*y
           2 * x;              // 2x
}

double max(double a, double b)
{
    if (a > b)
        return a;

    return b;
}

Vector solveLinearSystem(double x, double y, double epsX = 1e-8, double epsF = 1e-8)
{
    double hx{};
    double hy{};
    int i{0};
    do
    {
        double w{fx(x, y) * gy(x, y) - gx(x, y) * fy(x, y)};
        double wx{g(x, y) * fy(x, y) - f(x, y) * gy(x, y)};
        double wy{f(x, y) * gx(x, y) - g(x, y) * fx(x, y)};
        hx = wx / w;
        hy = wy / w;
        double xp{x};
        double yp{y};
        x = xp + hx;
        y = yp + hy;
        i++;
    } while (max(fabs(hx), fabs(hy)) > epsX && max(fabs(f(x, y)), fabs(g(x, y))) > epsF);

    Vector solution{};
    solution.push_back(x);
    solution.push_back(y);

    std::cout << "\n\nRozwiązanie po " << i << " iteracjach: \n|x|\n|y|";
    return solution;
}

int main()
{
    // const double epsX{1e-8};
    // const double epsF{1e-8};

    double x0{-1.0};
    double y0{-1.5};

    printVector(solveLinearSystem(x0, y0));

    double x1{4.0};
    double y1{-2.5};
    printVector(solveLinearSystem(x1, y1));

    return 0;
}