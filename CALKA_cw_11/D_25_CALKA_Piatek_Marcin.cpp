// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <utility>
#include <iomanip>

double f(double x)
{
    return sin(sqrt(x));
}

// Funkcja obliczająca całki oznaczone metodą trapezów oraz zliczająca liczbę kroków
std::pair<double, int> calkaTrapezow(double a, double b, double eps = 1e-4)
{

    int k{0};        // Numer kroku
    int m{1};        // Ilość przedziałów
    double h{b - a}; // Długość przedziału

    double f_a{f(a)}; // Wartość funkcji f_a
    double f_b{f(b)}; // Wartość funkcji f_b

    double T_nowe{(h / 2.) * (f_a + f_b)};
    double T_stare{};

    do
    {
        m = 2 * m;
        k++;
        h = h / 2.;

        // Wyliczanie sumy
        double suma{};
        for (int i = 1; i < m; ++i)
        {
            suma += f(a + i * h); // x = a+(i*h); y_i = f(x_i)
        }

        T_stare = T_nowe;
        T_nowe = h / 2. * (f_a + f_b + 2.0 * suma);

    } while (fabs(T_nowe - T_stare) >= eps);

    return std::make_pair(T_nowe, k);
}

// Funkcja obliczająca całki oznaczone metodą Simpsona oraz zliczająca liczbę kroków
std::pair<double, int> calkaSimpsona(double a, double b, double eps = 1e-4)
{

    int k{1};                // Numer kroku
    int m{2};                // Ilość przedziałów
    double h{(b - a) / 2.0}; // Długość przedziału

    double f_a{f(a)}; // Wartość funkcji f_a
    double f_b{f(b)}; // Wartość funkcji f_b

    double S_nowe{(h / 3.0) * (f_a + 4. * f((a + b) / 2.0) + f_b)};
    double S_stare{};

    do
    {
        m = 2 * m;
        k++;
        h = h / 2.;

        // Wyliczanie sum
        double sumaParzyste{}, sumaNieparzyste{};
        for (int i = 1; i < m; ++i)
        {
            i % 2 == 0 ? sumaParzyste += f(a + i * h) : sumaNieparzyste += f(a + i * h);
        }

        S_stare = S_nowe;

        S_nowe = h / 3. * (f_a + f_b + 4. * sumaNieparzyste + 2. * sumaParzyste);
    } while (fabs(S_nowe - S_stare) >= eps);

    return std::make_pair(S_nowe, k);
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

    // Wypisanie wyników metody trapezów

    std::pair<double, int> MetodaTrapezow = {calkaTrapezow(a, b)}; // Wartości całki (first) oraz liczba kroków (second)

    std::cout << "\n\t Metoda trapezów";
    std::cout << "\n\t\t wartość całki: " << MetodaTrapezow.first;
    std::cout << "\n\t\t liczba kroków: " << MetodaTrapezow.second;

    // Wypisanie wyników metody Simpsona

    // Równoważne do std::pair<double, int> z metody trapezów powyżej
    auto [integral, steps] = calkaSimpsona(a, b); // Wartości całki oraz liczba kroków

    std::cout << "\n\t Metoda Simpsona";
    std::cout << "\n\t\t wartość całki: " << integral;
    std::cout << "\n\t\t liczba kroków: " << steps;

    std::cout << "\n\n\n";
    return 0;
}