// Visual Studio Code 1.105.1

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Funkcja wypisująca wielomian
void printPolynomial(const Vector &a);

// Funkcja do wypisywania macierzy
void printMatrix(const Matrix &mtr);

// Funkcja implementująca Gaussa-Jordana
Vector gaussJordan(Matrix matrix, Vector vector);

// Funkcja budująca macierz układu równań dla metody sklejania
Matrix buildSplineMatrix(int n, double h, double alfa, double beta);

// Funkcja budująca wektor prawej strony
Vector buildRightHandSide(const Vector &x, const Vector &y, double h, double alfa, double beta);

// Funkcja obliczająca wartość funkcji sklejanej S(x)
double evaluateSpline(double x, int n, double h, const Vector &coefficients, double alfa, double beta, double a);

// Funkcja implementująca interpolację metodą sklejania
Vector splineInterpolation(const Vector &x, const Vector &y, double alfa, double beta);

int main()
{
    /*
    std::cout << std::fixed << std::setprecision(6);

    // Dane wejściowe - zad. 1
    Vector x = {1.5, 2.0, 2.5, 3.5, 3.8, 4.1};
    Vector y = {2.0, 5.0, -1.0, 0.5, 3.0, 7.0};

    // Dane wejściowe - zad. 2
    /*
    Vector x = {-2.0, -0.5, 1.2, 3.0, 3.5, 5.0, 5.5};
    Vector y = {7.0, 5.0, 1.0, -0.5, 2.0, 1.0, -1.0};
    */

    std::cout << std::fixed << std::setprecision(6);

    // Parametry zadania
    const int n = 20;         // liczba przedziałów
    const double a = -1.0;    // początek przedziału
    const double b = 1.0;     // koniec przedziału
    const double alfa = 0.5;  // warunek brzegowy: S'(a) = alfa
    const double beta = -0.5; // warunek brzegowy: S'(b) = beta

    // Obliczenie kroku
    double h = (b - a) / n;

    // Tworzenie węzłów interpolacji
    Vector x(n + 1);
    Vector y(n + 1);

    std::cout << "--- INTERPOLACJA METODĄ SKLEJANIA WIELOMIANÓW ---\n\n";
    std::cout << "Parametry:\n";
    std::cout << "  Przedział: [" << a << ", " << b << "]\n";
    std::cout << "  Liczba przedziałów: " << n << "\n";
    std::cout << "  Krok h = " << h << "\n";
    std::cout << "  Warunki brzegowe: S'(" << a << ") = " << alfa
              << ", S'(" << b << ") = " << beta << "\n\n";

    std::cout << "Funkcja interpolowana: f(x) = 1 / (1 + x^2)\n\n";

    // Wypełnienie węzłów
    std::cout << "Węzły interpolacji:\n";
    for (int i = 0; i <= n; ++i)
    {
        x[i] = a + i * h;
        y[i] = 1.0 / (1.0 + x[i] * x[i]);

        if (i < 5 || i > n - 3)
        {
            std::cout << "  x[" << std::setw(2) << i << "] = " << std::setw(8) << x[i]
                      << ",  y[" << std::setw(2) << i << "] = " << y[i] << "\n";
        }
        else if (i == 5)
        {
            std::cout << "  ...\n";
        }
    }

    try
    {
        // Obliczanie współczynników metodą sklejania
        Vector coefficients = splineInterpolation(x, y, alfa, beta);

        std::cout << "\n--- Współczynniki (C_0, C_1, ..., C_n) ---\n";
        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            if (i < 5 || i > coefficients.size() - 3)
            {
                std::cout << "C[" << std::setw(2) << i << "] = "
                          << std::setw(15) << coefficients[i] << "\n";
            }
            else if (i == 5)
            {
                std::cout << "...\n";
            }
        }

        // Weryfikacja w węzłach
        std::cout << "\n--- Weryfikacja w wybranych węzłach ---\n";
        std::cout << std::setw(8) << "x_i"
                  << std::setw(15) << "S(x_i)"
                  << std::setw(15) << "y_i"
                  << std::setw(18) << "błąd [S-y]\n";
        std::cout << std::string(56, '-') << "\n";

        for (int i = 0; i <= n; i += (n / 10))
        {
            double sxi = evaluateSpline(x[i], n, h, coefficients, alfa, beta, a);
            double error = sxi - y[i];
            std::cout << std::setw(8) << x[i]
                      << std::setw(15) << sxi
                      << std::setw(15) << y[i]
                      << std::setw(18) << error << "\n";
        }

        // Interakcja z użytkownikiem
        std::cout << "\n"
                  << std::string(60, '=') << "\n";
        std::cout << "----- Obliczanie wartości S(x) -----\n";
        std::cout << std::string(60, '=') << "\n";
        std::cout << "\nWspółczynniki są zapisane.\n";
        std::cout << "Można obliczać wart. S(x) dla wybranego x należącego [" << a << ", " << b << "]\n\n";

        char choice;
        do
        {
            std::cout << "Podaj wartość x (lub 'q' aby zakończyć): ";
            double userX;

            if (!(std::cin >> userX))
            {
                // Jeśli użytkownik podał coś innego niż liczbę
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                break;
            }

            // Sprawdzenie czy x jest w przedziale
            if (userX < a || userX > b)
            {
                std::cout << "\t UWAGA: x poza przedziałem interpolacji!\n";
            }

            // Oblicz wartość S(x)
            double result = evaluateSpline(userX, n, h, coefficients, alfa, beta, a);
            double actual = 1.0 / (1.0 + userX * userX);
            double error = result - actual;

            std::cout << "\tS(" << userX << ") = " << result << "\n";
            std::cout << "\tf(" << userX << ") = " << actual << "\n";
            std::cout << "\tbłąd = " << error << "\n\n";

            std::cout << "Czy chcesz obliczyć kolejną wartość? (t/n): ";
            std::cin >> choice;
            std::cout << "\n";

        } while (choice == 't' || choice == 'T');

        std::cout << "\nZakończono program.\n";
    }
    catch (const std::exception &ex)
    {
        std::cerr << "\n!!! BŁĄD: " << ex.what() << "\n";
    }

    return 0;
}

// Funkcja do wypisywania macierzy
void printMatrix(const Matrix &mtr)
{
    std::cout << std::setprecision(4);
    for (size_t i = 0; i < mtr.size(); i++)
    {
        std::cout << " | ";
        for (size_t j = 0; j < mtr[i].size(); j++)
        {
            std::cout << std::setw(10) << mtr[i][j] << " ";
        }
        std::cout << "|\n";
    }
    std::cout << std::setprecision(6);
}

// Funkcja rozwiązująca układ równań metodą eliminacji Gaussa-Jordana
Vector gaussJordan(Matrix matrix, Vector vector)
{
    int dimension = matrix.size();
    constexpr double eps = 1.0e-10;

    // Sprawdzenie wymiarów
    if (matrix[0].size() != dimension)
        throw std::invalid_argument("Macierz musi być kwadratowa!");

    if (vector.size() != dimension)
        throw std::invalid_argument("Wektor ma nieprawidłowy wymiar!");

    // Utworzenie macierzy rozszerzonej [A|b]
    // Dodajemy wektor y jako dodatkową kolumnę
    for (int i = 0; i < dimension; i++)
    {
        matrix[i].push_back(vector[i]);
    }

    for (int k = 0; k < dimension; k++)
    {
        // Znajdź wiersz z największym elementem
        double max = fabs(matrix[k][k]);
        int rowOfMax = k;

        for (int i = k + 1; i < dimension; i++)
        {
            if (fabs(matrix[i][k]) > max)
            {
                max = fabs(matrix[i][k]);
                rowOfMax = i;
            }
        }

        // Sprawdź czy macierz nie jest osobliwa
        if (max < eps)
            throw std::runtime_error("Macierz układu osobliwa - brak jednoznacznego rozwiązania!");

        // Zamień wiersze jeśli potrzeba
        if (rowOfMax != k)
        {
            std::swap(matrix[k], matrix[rowOfMax]);
        }

        // Normalizacja wiersza pivotowego
        double pivot = matrix[k][k];
        for (int j = k; j <= dimension; j++)
        {
            matrix[k][j] /= pivot;
        }

        // Eliminacja w pozostałych wierszach (Gauss-Jordan)
        for (int i = 0; i < dimension; i++)
        {
            if (i != k)
            {
                double factor = matrix[i][k];
                for (int j = k; j <= dimension; j++)
                {
                    matrix[i][j] -= factor * matrix[k][j];
                }
            }
        }
    }

    // Wyciągnij rozwiązanie (ostatnia kolumna)
    Vector solution(dimension);
    for (int i = 0; i < dimension; i++)
    {
        solution[i] = matrix[i][dimension];
    }

    return solution;
}

// Funkcja wypisująca wielomian
void printPolynomial(const Vector &a)
{
    int n = static_cast<int>(a.size()) - 1;

    std::cout << "P(x) = ";
    for (int i = 0; i <= n; ++i)
    {
        if (i > 0 && a[i] >= 0)
            std::cout << " + ";
        else if (i > 0 && a[i] < 0)
            std::cout << " ";

        std::cout << a[i] << "*u_" << i << "(x)";
    }
    std::cout << "\n";
}

// Funkcja budująca macierz układu równań dla metody sklejania
Matrix buildSplineMatrix(int n, double h, double alfa, double beta)
{
    // Tworzymy macierz trójdiagonalną (n+1) x (n+1)
    Matrix matrix(n + 1, Vector(n + 1, 0.0));

    for (int i = 0; i <= n; ++i)
    {
        if (i == 0)
        {
            // Pierwszy wiersz: warunek brzegowy na lewym końcu
            matrix[0][0] = 4.0;
            matrix[0][1] = 2.0;
        }
        else if (i == n)
        {
            // Ostatni wiersz: warunek brzegowy na prawym końcu
            matrix[n][n - 1] = 2.0;
            matrix[n][n] = 4.0;
        }
        else
        {
            // Wiersze środkowe: warunki ciągłości
            matrix[i][i - 1] = 1.0;
            matrix[i][i] = 4.0;
            matrix[i][i + 1] = 1.0;
        }
    }

    return matrix;
}

// Funkcja budująca wektor prawej strony
Vector buildRightHandSide(const Vector &x, const Vector &y, double h, double alfa, double beta)
{
    int n = static_cast<int>(x.size()) - 1;
    Vector rhs(n + 1);

    for (int i = 0; i <= n; ++i)
    {
        // Wartość funkcji w węźle
        rhs[i] = y[i];

        // Modyfikacja dla warunków brzegowych
        if (i == 0)
        {
            rhs[i] += (h / 3.0) * alfa;
        }

        if (i == n)
        {
            rhs[i] -= (h / 3.0) * beta;
        }
    }

    return rhs;
}

// Funkcja implementująca interpolację metodą sklejania
Vector splineInterpolation(const Vector &x, const Vector &y, double alfa, double beta)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("x oraz y powinny mieć równe wymiary");

    int n = static_cast<int>(x.size()) - 1;

    // Obliczenie kroku (zakładamy równomierne węzły)
    double h = (x[n] - x[0]) / n;

    // Zbuduj macierz układu równań
    Matrix matrix = buildSplineMatrix(n, h, alfa, beta);

    std::cout << "\n--- Macierz układu równań (fragment) ---\n";

    // Wypisz tylko fragment dla dużych macierzy
    if (n + 1 <= 10)
    {
        printMatrix(matrix);
    }
    else
    {
        std::cout << "Macierz " << (n + 1) << "x" << (n + 1) << " (zbyt duża do wyświetlenia)\n";
    }

    // Zbuduj wektor prawej strony
    Vector rhs = buildRightHandSide(x, y, h, alfa, beta);

    // Rozwiąż układ równań metodą Gaussa-Jordana
    Vector coefficients = gaussJordan(matrix, rhs);

    return coefficients;
}

// Funkcja obliczająca wartość funkcji sklejanej S(x)
double evaluateSpline(double x, int n, double h, const Vector &coefficients,
                      double alfa, double beta, double a)
{
    // Znajdź przedział, w którym leży x
    int m = static_cast<int>(std::floor((x - a) / h));

    // Obsługa przypadku brzegowego (x = b)
    if (m == n)
        m = n - 1;

    // Obliczenie X_m i parametru t
    double Xm = a + m * h;
    double t = (x - Xm) / h;

    // Obliczenie wartości A i B zgodnie z warunkami brzegowymi
    double A, B;

    if (m == 0)
    {
        // Lewy warunek brzegowy
        A = coefficients[1] - (h / 3.0) * alfa;
    }
    else
    {
        A = coefficients[m - 1];
    }

    if (m == n - 1)
    {
        // Prawy warunek brzegowy
        B = coefficients[n] + (h / 3.0) * beta;
    }
    else
    {
        B = coefficients[m + 2];
    }

    // Obliczenie wartości funkcji sklejanej w punkcie x
    // S(x) = A*(1-t)^3 + C_m*((2-t)^3 - 4*(1-t)^3) + C_{m+1}*((1+t)^3 - 4*t^3) + B*t^3
    double result = A * std::pow(1.0 - t, 3) +
                    coefficients[m] * (std::pow(2.0 - t, 3) - 4.0 * std::pow(1.0 - t, 3)) +
                    coefficients[m + 1] * (std::pow(1.0 + t, 3) - 4.0 * std::pow(t, 3)) +
                    B * std::pow(t, 3);

    return result;
}
