// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <utility>

double f(double x, double y = 0)
{
    return x + y;
};

// Wzory Rungego-Kutty - formuła czwartego rzędu, obliczanie szukanej wartości funkcji
std::pair<double, double> RungeKutty(double x, double y, double h)
{
    double k1 = h * f(x, y);
    double k2 = h * f(x + (h * 0.5), y + (k1 * 0.5));
    double k3 = h * f(x + (h * 0.5), y + (k2 * 0.5));
    double k4 = h * f(x + h, y + k3);

    y = y + (k1 + (2 * k2) + (2 * k3) + k4) / 6.;
    x = x + h;
    return std::make_pair(y, x);
}

double RozwiazRownanie(double x0, double y0, double h, double x)
{
    int n = ceil(((x - x0) / h));
    for (int i = 0; i < n; ++i)
    {
        y0 = RungeKutty(x0, y0, h).first;
        x0 = RungeKutty(x0, y0, h).second;
    }
    return y0;
};

int main()
{
    std::cout.precision(15);
    // Warunek początkowy
    double y0{1.};
    // Początek przedziału
    double x0{0.};
    // Punkt, dla którego szukamy wartość funkcji szukanej y
    double x{};
    double h{1e-2};

    std::cout << "\n\nWpisz x: ";
    std::cin >> x;
    std::cout << "\nWpisz h: ";
    std::cin >> h;

    double hdiff{3};
    for (int i = 0; i < hdiff; ++i)
    {
        h = h * 0.1;
        double y{RozwiazRownanie(x0, y0, h, x)};
        std::cout << "\n\n\t\th = " << h;
        std::cout << "\n\t\ty(" << x << ") = " << y;
    }
    return 0;
}