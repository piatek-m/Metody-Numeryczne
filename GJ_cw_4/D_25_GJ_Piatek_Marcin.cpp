// Visual Studio Code, wersja 1.105.0

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
            /* równoważne:
                for(int j = k; j <= dimension; j++)
                {
                    std::swap(matrix[k][j], matrix[rowOfMax][j])
                }
            */
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

            /* Sprawdzanie

            printMatrix(matrix);
            std::cout << "i=" << i << " k=" << k << "\n";

            */
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

// Funkcja pomocnicza do wyświetlania wyniku
void printSolution(const Vector &vector)
{
    std::cout << "\nRozwiązanie:\n";
    for (int i = 0; i < vector.size(); i++)
    {
        std::cout << "x[" << i << "] = " << vector[i] << "\n";
    }
}

int main()
{
    // Układ równań
    // 0x - y - 2z = -2
    // 2x + 2y + 1z = 3
    // 4x + y + 3z = -1

    Matrix A = {
        {0, -1, -2},
        {2, 2, 1},
        {4, 1, 3}};

    Vector b = {-2,
                3,
                -1};

    try
    {
        Vector solution = gaussJordan(A, b);
        printSolution(solution);
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nBłąd: " << e.what() << "\n";
    }

    Matrix A2 = {
        {3, 6},
        {2, 4},
    };

    Vector b2 = {1,
                 2};

    try
    {
        Vector solution = gaussJordan(A2, b2);
        printSolution(solution);
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nBłąd: " << e.what() << "\n";
    }

    Matrix A3 = {
        {0, 5},
        {-2, 3},
    };

    Vector b3 = {15,
                 7};

    try
    {
        Vector solution = gaussJordan(A3, b3);
        printSolution(solution);
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nBłąd: " << e.what() << "\n";
    }

    std::cout << "\n\n\n\n";
    return 0;
}
