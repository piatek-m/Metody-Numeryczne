// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Wyliczanie wartości funkcji g
double g(double x, unsigned j);

// Wyliczanie wartości funkcji h
double h(double x, unsigned k, const Vector &p);

// Wyliczanie wielomianu
double P(double x, unsigned k, const Vector &p);

// Funkcja obliczająca całki oznaczone metodą trapezów oraz zliczająca liczbę kroków
std::pair<double, int> calkaTrapezow(std::function<double(double)> func, double a, double b, double eps = 1.0e-5); // Przyjmuje jako argument funkcję, aby możliwe było wygodne definiowanie funkcji całkowanej

// Wyliczanie funkcji f
double f(double x);

// Wartość funkcji błędu aproksymacji
double S(unsigned k, double a, double b, const Vector &p);

// Całkowanie funkcji g
double I_j(unsigned j, double a, double b);

// Budowa macierzy układu równań
std::pair<Matrix, Vector> buildMatrix(unsigned k, double a, double b);

// Rozwiązanie układu równań metodą Gaussa-Jordana
Vector gaussJordan(Matrix matrix, Vector vector);

int main()
{
    std::cout << std::setprecision(20);

    double a{}, b{};
    std::cout << "\n\nPodaj początek przedziału a: ";
    std::cin >> a;
    std::cout << "\nPodaj koniec przedziału b: ";
    std::cin >> b;

    unsigned k{};
    std::cout << "\nPodaj stopień wielomianu k: ";
    std::cin >> k;

    // Budowa układu równań
    auto [macierz, wektor_I] = buildMatrix(k, a, b);

    // Rozwiązanie układu
    Vector p = gaussJordan(macierz, wektor_I);

    // Wypisanie współczynników
    std::cout << "\n\nWspółczynniki wielomianu P_" << k << "(x):\n";
    for (unsigned j = 0; j <= k; ++j)
    {
        std::cout << "\t\tp" << j << " = " << p[j] << "\n";
    }

    // Obliczenie błędu aproksymacji
    double error = S(k, a, b, p);
    std::cout << "\nBłąd aproksymacji S(" << k << ") = " << error << "\n";

    // Testowanie wartości wielomianu
    char wybor = 't';
    while (wybor == 't' || wybor == 'T')
    {
        double x;
        std::cout << "\nPodaj x = ";
        std::cin >> x;

        std::cout << "P_" << k << "(" << x << ") = " << P(x, k, p) << "\n";

        std::cout << "\nNastępne x? (t/n): ";
        std::cin >> wybor;
    }

    return 0;
}

double g(double x, unsigned j)
{
    if (j == 0)
        return f(x);
    return f(x) * pow(x, j);
}

double P(double x, unsigned k, const Vector &p)
{
    double sum{0};

    for (int j = 0; j <= k; ++j)
        sum += p[j] * pow(x, j);

    return sum;
}

double h(double x, unsigned k, const Vector &p)
{
    double blad = f(x) - P(x, k, p);

    return blad * blad;
}

std::pair<double, int> calkaTrapezow(std::function<double(double)> func, double a, double b, double eps)
{

    int krok{0};     // Numer kroku
    int m{1};        // Ilość przedziałów
    double h{b - a}; // Długość przedziału

    double f_a{func(a)}; // Wartość funkcji f_a
    double f_b{func(b)}; // Wartość funkcji f_b

    double T_nowe{(h / 2.) * (f_a + f_b)};
    double T_stare{};

    do
    {
        m = 2 * m;
        krok++;
        h = h / 2.;

        // Wyliczanie sumy
        double suma{};
        for (int i = 1; i < m; ++i)
        {
            suma += func(a + i * h); // x = a+(i*h); y_i = f(x_i)
        }

        T_stare = T_nowe;
        T_nowe = h / 2. * (f_a + f_b + 2.0 * suma);

    } while (fabs(T_nowe - T_stare) >= eps);

    return std::make_pair(T_nowe, krok);
}

double f(double x)
{
    // return exp(-x);

    // funkcje testowe
    // return -1 + 2 * x;
    return 3 - 2 * x + x * x;
}

double S(unsigned k, double a, double b, const Vector &p)
{
    // Lambda funkcji h, aby działało "podkładanie" funkcji pod całkę
    auto func_h = [k, &p](double x)
    {
        return h(x, k, p);
    };

    return calkaTrapezow(func_h, a, b).first;
}

double I_j(unsigned j, double a, double b)
{
    // Lambda funkcji g, aby działało "podkładanie" funkcji pod całkę
    auto func_g = [j](double x)
    {
        return g(x, j);
    };

    return calkaTrapezow(func_g, a, b).first;
};

std::pair<Matrix, Vector> buildMatrix(unsigned k, double a, double b)
{

    // Tworzenie macierzy kwadratowej (k+1) x (k+1)
    Matrix Mtrx(k + 1, Vector(k + 1));

    // Wypełnienie macierzy
    for (int row = 0; row <= k; ++row)
    {
        for (int col = 0; col <= k; ++col)
        {
            // Wykładnik i = j+1, ..., j+(k+1)
            int i{row + col};
            if (i == 0)
                Mtrx[row][col] = b - a;
            else
                Mtrx[row][col] = (pow(b, i + 1) - pow(a, i + 1)) / (i + 1);
        }
    }

    // Wektor b macierzy rozszerzonej [A|b]
    Vector Vec(k + 1);

    // b = I_j = całka z [x^j * f(x)]dx
    for (unsigned j = 0; j <= k; ++j)
    {
        Vec[j] = I_j(j, a, b);
    }

    return std::make_pair(Mtrx, Vec);
}

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