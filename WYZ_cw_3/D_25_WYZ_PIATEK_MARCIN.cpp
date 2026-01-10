// Visual Studio Code, wersja 1.105.0
#include <iostream>
#include <vector>
#include <math.h>
#include <utility>

// Alias dla typu wektora dwuwymiarowego
using Matrix = std::vector<std::vector<double>>;

double wyznacznik(Matrix &matrix)
{
    // Jeśli macierz nie jest kwadratowa to przerwij program
    if (matrix.size() != matrix[0].size())
        throw std::invalid_argument("Dla macierzy niekwadratowej nie można wyliczyć wyznacznika!");

    int mtrDimension = matrix.size(); // Wymiar(y) macierzy

    int sign{1};   // Znak wyznacznika
    double max{};  // Największa wartość w kolumnie
    int r{};       // Indeks wiersza z największym elementem
    double det{1}; // Wyznacznik
    double pom{};

    // Wybór elem. głównego - p
    for (int k = 0; k < mtrDimension; k++)
    {
        max = matrix[k][k];
        r = k;
        for (int i = k + 1; i < mtrDimension; i++)
        {
            if (abs(matrix[i][k]) > abs(max))
            {
                max = matrix[i][k];
                r = i;
            }
        }

        // Czy macierz jest osobliwa
        if (max == 0)
        {
            return det = 0;
        }

        // Jeśli największy elem. nie jest w wierszu k, to zamieniamy wiersze
        if (r != k)
        {
            sign = -sign; // Zmiana znaku przy zmianie wierszy
            for (int j = k; j < mtrDimension; j++)
            {
                std::swap(matrix[k][j], matrix[r][j]);
            }
        }

        // Zerowanie poniżej pivota
        for (int i = k + 1; i < mtrDimension; i++)
        {
            pom = matrix[i][k];

            for (int j = k; j < mtrDimension; j++)
                matrix[i][j] = matrix[i][j] - (pom * matrix[k][j] / max);
        }
    }

    // Wyliczanie wyznacznika po przekątnej
    for (int k = 0; k < mtrDimension; k++)
        det = det * matrix[k][k];

    det *= sign;

    return det;
}

int main()
{
    Matrix A1 = {{3, 4},
                 {2, 5}};
    std::cout << "det A1 = " << wyznacznik(A1) << std::endl;

    Matrix A2 = {{2, 0, 0},
                 {0, 5, 0},
                 {7, 8, -3}};
    std::cout << "det A2 = " << wyznacznik(A2) << std::endl;

    Matrix A3 = {{6, 4},
                 {3, 2}};
    std::cout << "det A3 = " << wyznacznik(A3) << std::endl;

    Matrix A4 = {{1, -1, 2, 2, 0},
                 {0, 0.5, 1, 3, 0},
                 {-3, -2, 0, 5, 0.5},
                 {0, 2, 1, 3, 0},
                 {2, 1, 0, -4, 2}};
    std::cout << "det A4 = " << wyznacznik(A4) << std::endl;

    return 0;
}