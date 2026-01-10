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

// Funkcja obliczająca wartość fi_k(x) według definicji rekurencyjnej
double fi(double x, unsigned k);

// Funkcja obliczająca wartość S(p_0, ..., p_m) = suma [P(x_j) - y_j]^2
double S(Vector x, Vector y, unsigned n, unsigned m, const Vector &a);

// Tworzenie macierzy układu korzystając z definicji rekurencyjnej funkcji bazowych
Matrix buildAugmentedMtrx(const Vector &x, unsigned m, unsigned n);

// Obliczanie wartości P(x) = suma p_i * fi_i(x)
double evaluatePolynomial(double x, unsigned m, const Vector &p);

// Funkcja implementująca aproksymacje wielomianu
Vector aprroximatePolynomial(const Vector &x, const Vector &y, const int m, const int n);

int main()
{
    std::cout << std::fixed << std::setprecision(20);

    // Dane wejściowe - zad. 1
    Vector x = {1.0, 1.5, 2.0, 4.0, 6.0, 7.0};
    Vector y = {-0.7, -0.35, -0.05, 1.45, 3.95, 5.35};

    const int n{static_cast<int>(x.size())};

    std::cout << "--- APROKSYMACJA FUNKCJI ---\n\n";
    std::cout << "Liczba punktów: " << x.size() << "\n";
    std::cout << "Punkty (węzły):\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "  i=" << i << ": (" << x[i] << ", " << y[i] << ")\n";
    }

    std::cout << "\nWprowadź m: ";
    int m{};
    std::cin >> m;
    std::cout << "\n";

    try
    {
        // Obliczanie współczynników i zapisanie ich w wektorze p
        Vector p = aprroximatePolynomial(x, y, m, n);

        std::cout << "\n--- Współczynniki (a_0, a_1, ..., a_n) ---\n";
        for (size_t i = 0; i < p.size(); ++i)
        {
            std::cout << "a[" << i << "] = " << std::setw(15) << p[i] << "\n";
        }

        // Sprawdzenie w węzłach
        std::cout << "\n--- Weryfikacja w węzłach ---\n";
        std::cout << std::setw(8) << "x_i"
                  << std::setw(15) << "P(x_i)"
                  << std::setw(15) << "y_i"
                  << std::setw(18) << "błąd S\n";
        std::cout << std::string(56, '-') << "\n";

        int n = static_cast<int>(p.size()) - 1;
        for (size_t i = 0; i < x.size(); ++i)
        {
            double pxi = evaluatePolynomial(x[i], n, p);
            double error = S(x, y, n, m, p);
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

            // Oblicz wartość P(x) używając zapisanych współczynników p
            double result = evaluatePolynomial(userX, n, p);
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

double S(Vector x, Vector y, unsigned n, unsigned m, const Vector &a)
{
    double sum{};

    for (int j = 0; j < n; ++j)
    {
        sum += (evaluatePolynomial(x[j], m, a) - y[j]) * (evaluatePolynomial(x[j], m, a) - y[j]);
    }

    return sum;
};

Vector aprroximatePolynomial(const Vector &x, const Vector &y, const int m, const int n)
{
    // Zbuduj macierz układu równań
    Matrix mtrx = buildAugmentedMtrx(x, m, n);

    Vector b(m + 1);
    for (int i = 0; i <= m; i++)
    {
        double sum{0};
        for (int t = 0; t < n; t++)
        {
            b[i] += 1.0 * y[t] * fi(x[t], i);
        }
    }

    for (int i = 0; i < m + 1; i++)
        mtrx[i][m + 1] = b[i];

    std::cout << "\n--- Macierz układu równań ---\n";
    printMatrix(mtrx);

    // Rozwiąż układ równań metodą Gaussa-Jordana
    Vector wspolczynniki = gaussJordan(mtrx, y);

    return wspolczynniki;
}

double evaluatePolynomial(double x, unsigned m, const Vector &p)
{

    double result = 0.0;

    // P(x) = p_0 * fi_0(x) + p1_1 * fi_1(x) + ... + p_m * fi_m(x)
    for (unsigned k = 0; k <= m; ++k)
    {
        result += p[k] * fi(x, k);
    }

    return result;
}

double fi(double x, unsigned k)
{

    if (k == 0)
        return 2.0;

    if (k == 1)
        return x + 1;

    // Dla k≥2 używamy rekurencji iteracyjnej
    double fi_k{};

    fi_k = (k - 1.0 / k + 1.0) * x * fi(x, k - 1) -
           1.0 / k * fi(x, k - 2);

    return fi_k;
}

Matrix buildAugmentedMtrx(const Vector &x, unsigned m, unsigned n)
{
    // Tworzenie macierzy kwadratowej (n+1) x (n+1)
    // Ostatnia kolumna dopisywana jest w funkcji realizującej metodę Gaussa-Jordana
    Matrix augMtrx(m + 1, Vector(m + 1));

    // Wypełnienie macierzy
    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= m; ++j)
        {
            double sum = 0;
            for (int t = 0; t < n; ++t)
            {
                // 1.0 = w(x)
                sum += fi(x[t], i) * fi(x[t], j);
            }
            augMtrx[i][j] = fi(x[i], j);
        }
    }

    return augMtrx;
}