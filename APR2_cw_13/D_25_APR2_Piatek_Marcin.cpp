// Visual Studio Code 1.105.1

#include <iostream>
#include <math.h>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Wyliczanie funkcji g
double g(double x, u_int j)
{
    // Dla j=0: g(x,j) = f(x)
    if (j == 0)
        return f(x);

    // Dla j>0: g(x) = f(x) * x^i
    double f_x = f(x);
    for (int i = 0; i < j; ++i)
    {
        f_x *= x;
    }
    return f_x;
}

// Wyliczanie funkcji h
double h(double x, u_int k, Vector &p)
{
    double temp = f(x) - P(x, k, p); // Wartość wyrażenia wewnętrznęgo

    // Wartość funkcji to kwadrat wyrażenia wewnętrznego
    return temp * temp;
};

double P(double x, u_int k, Vector &p)
{
    double sum{0};

    for (int j = 0; j < k; ++j)
    {
        // Wyliczanie x^j
        double x_ToPow_j{x}; // x do potęgi j-tej
        for (int i = 0; i < j; ++i)
        {
            x_ToPow_j *= x;
        }

        sum += p[j] * x_ToPow_j;
    }

    return sum;
}

// Funkcja obliczająca całki oznaczone metodą trapezów oraz zliczająca liczbę kroków
std::pair<double, int> calka_h(double a, double b, double x, u_int k, Vector &p, double eps = 1.0e-5)
{

    int krok{0};            // Numer kroku
    int m{1};               // Ilość przedziałów
    double interval{b - a}; // Długość przedziału

    double f_a{h(a, k, p)}; // Wartość funkcji f_a
    double f_b{h(b, k, p)}; // Wartość funkcji f_b

    double T_nowe{(interval / 2.) * (f_a + f_b)};
    double T_stare{};

    do
    {
        m = 2 * m;
        krok++;
        interval = interval / 2.;

        // Wyliczanie sumy
        double suma{};
        for (int i = 1; i < m; ++i)
        {
            suma += f(a + i * interval); // x = a+(i*h); y_i = f(x_i)
        }

        T_stare = T_nowe;
        T_nowe = interval / 2. * (f_a + f_b + 2.0 * suma);

    } while (fabs(T_nowe - T_stare) >= eps);

    return std::make_pair(T_nowe, krok);
}

std::pair<double, int> calka_g(double a, double b, double x, u_int j, double eps = 1.0e-5)
{

    int krok{0};            // Numer kroku
    int m{1};               // Ilość przedziałów
    double interval{b - a}; // Długość przedziału

    double f_a{g(a, j)}; // Wartość funkcji f_a
    double f_b{g(b, j)}; // Wartość funkcji f_b

    double T_nowe{(interval / 2.) * (f_a + f_b)};
    double T_stare{};

    do
    {
        m = 2 * m;
        krok++;
        interval = interval / 2.;

        // Wyliczanie sumy
        double suma{};
        for (int i = 1; i < m; ++i)
        {
            suma += f(a + i * interval); // x = a+(i*h); y_i = f(x_i)
        }

        T_stare = T_nowe;
        T_nowe = interval / 2. * (f_a + f_b + 2.0 * suma);

    } while (fabs(T_nowe - T_stare) >= eps);

    return std::make_pair(T_nowe, krok);
}

// Wyliczanie funkcji f
double f(double x)
{
    return exp(-x);

    /* funkcje testowe
    return -1+2*x;
    return 3-2*x+x*x;
    */
}

// Wartość funkcji błędu aproksymacji
double S(double a, double b, double x, u_int k, Vector &p)
{
    return calka_h(a, b, x, k, p).first;
}

double I_j(double a, double b, double x, u_int j)
{
    return calka_g(a, b, x, j).first;
};

Matrix buildAugmentedMtrx(const Vector &x, unsigned m, unsigned n, double a, double b)
{
    // Tworzenie macierzy kwadratowej (n+1) x (m+1)
    // Ostatnia kolumna dopisywana jest w funkcji realizującej metodę Gaussa-Jordana
    Matrix augMtrx(n + 1, Vector(m + 1));

    // Wypełnienie macierzy
    for (int i = 0; i <= m; ++i)
    {
        double a_i{0}, b_i{0};
        for (int z = 0; z < i; ++z)
        {
            a_i *= a;
            b_i *= b;
        }

        for (int j = 0; j <= m; ++j)
        {
            double delta_i{b_i - a_i};
            augMtrx[i][j] = delta_i / n + 1;
        }
    }

    return augMtrx;
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

int main()
{
    constexpr double eps{1.0e-5};

    double a{}, b{};
    std::cout << "\n\nPodaj koniec przedziału a: ";
    std::cin >> a;
    std::cout << "\tPodaj b: ";
    std::cin >> b;
    u_int k{};
    std::cout << "\nPodaj stopień k wielomianu aproksymującego P_k: ";
    std::cin >> k;

    double x{};
    char wybor = 't';
    while (wybor == 't' || wybor == 'T')
    {
        std::cout << "\nPodaj x = ";
        std::cin >> x;
        std::cout << "P(" << x << ") = " << P(x, k, p);
        << "\n";
        std::cout << "\nNastepne x? (t/n): ";
        std::cin >> wybor;
    }

    return 0;
}