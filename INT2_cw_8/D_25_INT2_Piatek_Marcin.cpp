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

// Funkcja obliczająca wartość u_k(x) według definicji rekurencyjnej
double u(double x, unsigned k)
{
    // Definicja rekurencyjna:
    // u_0(x) = 1
    // u_1(x) = 2x
    // u_k(x) = 2x * u_[k-1](x) - 2(k-1) * u_[k-2](x) dla k>=2

    if (k == 0)
        return 1.0;

    if (k == 1)
        return 2.0 * x;

    // Dla k≥2 używamy rekurencji iteracyjnej
    double u_k[3];

    u_k[2] = 1.0;     // u_[k-2] = u_0
    u_k[1] = 2.0 * x; // u_[k-1] = u_1
    u_k[0] = 0.0;

    for (unsigned j = 2; j <= k; ++j)
    {
        u_k[0] = 2.0 * x * u_k[1] - 2.0 * (j - 1) * u_k[2];
        u_k[2] = u_k[1];
        u_k[1] = u_k[0];
    }

    return u_k[0];
}

// Tworzenie macierzy układu korzystając z definicji rekurencyjnej funkcji bazowych
Matrix buildAugmentedMtrx(const Vector &x)
{
    int n = static_cast<int>(x.size()) - 1; // stopień wielomianu

    // Tworzenie macierzy kwadratowej (n+1) x (n+1)
    // Ostatnia kolumna dopisywana jest w funkcji realizującej metodę Gaussa-Jordana
    Matrix augMtrx(n + 1, Vector(n + 1));

    // Wypełnienie macierzy: A[i][j] = u_j(x_i)
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            augMtrx[i][j] = u(x[i], j);
        }
    }

    return augMtrx;
}

// Obliczanie wartości P(x) = suma a_k * u_k(x)
double evaluatePolynomial(double x, unsigned n, const Vector &a)
{

    double result = 0.0;

    // P(x) = a_0 * u_0(x) + a_1 * u_1(x) + ... + a_n * u_n(x)
    for (unsigned k = 0; k <= n; ++k)
    {
        result += a[k] * u(x, k);
    }

    return result;
}

// Funkcja implementująca interpolację wielomianu
Vector interpolatePolynomial(const Vector &x, const Vector &y)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("x oraz y powinny mieć równe wymiary");

    // Zbuduj macierz układu równań
    Matrix mtrx = buildAugmentedMtrx(x);

    std::cout << "\n--- Macierz układu równań ---\n";
    std::cout << "Struktura: A[i][j] = u_j(x_i)\n";
    printMatrix(mtrx);

    // Rozwiąż układ równań metodą Gaussa-Jordana
    Vector wspolczynniki = gaussJordan(mtrx, y);

    return wspolczynniki;
}

int main()
{
    std::cout << std::fixed << std::setprecision(6);

    // Dane wejściowe - zad. 1
    Vector x = {1.5, 2.0, 2.5, 3.5, 3.8, 4.1};
    Vector y = {2.0, 5.0, -1.0, 0.5, 3.0, 7.0};

    // Dane wejściowe - zad. 2
    /*
    Vector x = {-2.0, -0.5, 1.2, 3.0, 3.5, 5.0, 5.5};
    Vector y = {7.0, 5.0, 1.0, -0.5, 2.0, 1.0, -1.0};
    */

    std::cout << "--- INTERPOLACJA FUNKCJI ---\n\n";
    std::cout << "Liczba punktów: " << x.size() << "\n";
    std::cout << "Stopień wielomianu: " << (x.size() - 1) << "\n\n";
    std::cout << "Punkty (węzły):\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "  i=" << i << ": (" << x[i] << ", " << y[i] << ")\n";
    }

    try
    {
        // Obliczanie współczynników i zapisanie ich w wektorze a
        Vector a = interpolatePolynomial(x, y);

        std::cout << "\n--- Współczynniki (a_0, a_1, ..., a_n) ---\n";
        for (size_t i = 0; i < a.size(); ++i)
        {
            std::cout << "a[" << i << "] = " << std::setw(15) << a[i] << "\n";
        }

        std::cout << "\n";
        printPolynomial(a);

        // Sprawdzenie w węzłach
        std::cout << "\n--- Weryfikacja w węzłach ---\n";
        std::cout << std::setw(8) << "x_i"
                  << std::setw(15) << "P(x_i)"
                  << std::setw(15) << "y_i"
                  << std::setw(18) << "błąd [P-y]\n";
        std::cout << std::string(56, '-') << "\n";

        int n = static_cast<int>(a.size()) - 1;
        for (size_t i = 0; i < x.size(); ++i)
        {
            double pxi = evaluatePolynomial(x[i], n, a);
            double error = pxi - y[i];
            std::cout << std::setw(8) << x[i]
                      << std::setw(15) << pxi
                      << std::setw(15) << y[i]
                      << std::setw(18) << error << "\n";
        }

        // Interakcja z użytkownikiem
        std::cout << "\n"
                  << std::string(60, '=') << "\n";
        std::cout << "----- Obliczanie wartości P(x) -----\n";
        std::cout << std::string(60, '=') << "\n";
        std::cout << "\nWspółczynniki są zapisane.\n";
        std::cout << "Można obliczać wart. P(x) dla wybranego x\n\n";

        char choice;
        do
        {
            std::cout << "Podaj wartość x: ";
            double userX;

            if (!(std::cin >> userX))
            {
                // Jeśli użytkownik podał coś innego niż liczbę
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                break;
            }

            // Oblicz wartość P(x) używając zapisanych współczynników a
            double result = evaluatePolynomial(userX, n, a);
            std::cout << "\tP(" << userX << ") = " << result << "\n\n";

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