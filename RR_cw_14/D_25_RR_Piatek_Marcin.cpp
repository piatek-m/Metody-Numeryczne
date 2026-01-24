// Visual Studio Code 1.105.1

#include <iostream>
#include <cmath>
#include <utility>
#include <iomanip>

// Podana funkcja: y'(x) = f(x,y) = x+y
// Z warunkiem początkowym: y(x0) = y0 = 1
double f(double x, double y = 0)
{
    return x + y;
};

// Metoda Rungego-Kutty czwartego rzędu, obliczanie szukanej wartości funkcji.
// Przyjmujemy referencję do x, aby móc zmienić x_n na x_n+1, bez jawnego zwracania x_n
double RungeKutta(double &x, double y, double h)
{
    double k1 = h * f(x, y);
    double k2 = h * f(x + (h * 0.5), y + (k1 * 0.5));
    double k3 = h * f(x + (h * 0.5), y + (k2 * 0.5));
    double k4 = h * f(x + h, y + k3);

    // Wykorzystanie referencji: x_n+1 = x_n + h
    x += h;

    // Wyliczenie y_n+1
    y += (k1 + (2 * k2) + (2 * k3) + k4) / 6.;

    // Zwracamy y_n+1
    return y;
}

double RozwiazRownanie(double x0, double y0, double h, double x)
{
    // Wyliczenie liczby kroków
    int n = ceil(((x - x0) / h));
    // Poprawka kroku
    h = (x - x0) / n;

    // Wyliczanie wartości y
    for (int i = 0; i < n; ++i)
        y0 = RungeKutta(x0, y0, h);

    // Zwrócenie wartości funkcji szukanej y(x_n) = y_n
    return y0;
};

int main()
{
    std::cout.precision(16);

    // Warunek początkowy
    double y0{1.};
    // Początek przedziału
    double x0{0.};
    // Punkt, dla którego wyliczamy wartość funkcji szukanej y
    double x{1.};
    // Pierwszy krok, dla którego wyliczać wartość y
    double h_poczatkowe{1e-2};
    // Ostatni krok, dla którego wyliczać wartość y
    double h_koncowe{1e-5};

    // Pobieranie wartości od użytkownika (wartości domyślne kiedy nie wpiszemy nic)
    auto getInput = [](const std::string &tekst, double &wartosc)
    {
        std::string input;
        std::cout << "\n\n"
                  << tekst << " (domyslnie " << wartosc << "):\t";
        std::getline(std::cin, input);

        if (!input.empty())
            std::stringstream(input) >> wartosc;

        std::cout << "\tPrzyjęto: " << wartosc;
    };

    getInput("Wpisz x", x);
    getInput("Wpisz h poczatkowe", h_poczatkowe);
    getInput("Wpisz h koncowe", h_koncowe);

    // Zamiana końców przedziału kroku, aby zmniejszana była wyłącznie wartość większa
    if (h_poczatkowe < h_koncowe)
        std::swap(h_poczatkowe, h_koncowe);

    // Wyliczenie różnicy wykładników, aby uzyskać odpowiednią ilość iteracji
    unsigned int exponent_diff{static_cast<unsigned int>(log10(h_poczatkowe) - log10(h_koncowe))};

    //  Tablice wartości do porównania zależności wyników y od różnych kroków h
    double y_values[2];

    // Referencja w celu poprawy czytelności pętli poniżej
    double &h = h_poczatkowe;

    // Wyliczenie i wypisanie wartości y dla zakresu kroków
    for (unsigned int i = 0; i <= exponent_diff; ++i)
    {
        // Zapamiętanie obecnego y
        y_values[1] = RozwiazRownanie(x0, y0, h, x);

        // h wypisujemy w notacji naukowej
        std::cout << std::setprecision(0) << std::scientific;
        std::cout << "\n\n\n\th = " << h;

        // f wypisujemy w postaci dziesiętnej
        std::cout << std::setprecision(16) << std::defaultfloat;
        std::cout << "\n\t\ty(" << x << ") = " << y_values[1];

        if (i > 0)
        {
            // Różnicę wypisujemy w notacji naukowej
            std::cout << std::setprecision(2) << std::scientific;
            std::cout << "\n\t\tzmiana dy: = " << fabs(y_values[0] - y_values[1]);
        }
        // Zapamiętanie poprzedniego y
        y_values[0] = y_values[1];
        // Zmniejszenie wykladnika kroku
        h *= 0.1;
    }

    std::cout << "\n\n\n";
    return 0;
}