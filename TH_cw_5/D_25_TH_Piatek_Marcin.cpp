// Visual Studio Code, wersja 1.105.0

#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

/*
    Metoda Thomasa do rozwiązywania trójdiagonalnego układu równań liniowych.

    Układ:
    a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i], dla i = 1, 2, ..., n
    gdzie a[1] = 0 i c[n] = 0
     Parametry:
    a - dolna przekątna (a[0] powinno być 0)
    b - główna przekątna
    c - górna przekątna (c[n-1] powinno być 0)
    d - wektor prawych stron
    x - wektor rozwiązań (wynik)
*/
void thomasMethod(const vector<double> &a, const vector<double> &b,
                  const vector<double> &c, const vector<double> &d,
                  vector<double> &x)
{
    int n = d.size();

    // Wektory pomocnicze dla współczynników beta i gamma
    vector<double> beta(n);
    vector<double> gamma(n);

    // Obliczanie współczynników beta i gamma (przejście "w przód")

    // Warunki początkowe dla i=0
    beta[0] = -c[0] / b[0];
    gamma[0] = d[0] / b[0];

    // Dla i = 1, 2, ..., n-1 (odpowiada i=2,3,...,n)
    for (int i = 1; i < n; i++)
    {
        double denominator = a[i] * beta[i - 1] + b[i];
        beta[i] = -c[i] / denominator;
        gamma[i] = (d[i] - a[i] * gamma[i - 1]) / denominator;
    }

    // Obliczanie rozwiązań x (przejście "od tyłu")

    // Warunek końcowy: x[n-1] = gamma[n-1]
    x[n - 1] = gamma[n - 1];

    // Dla i = n-2, n-3, ..., 0
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = beta[i] * x[i + 1] + gamma[i];
    }
}

// Funkcja pomocnicza do wyświetlania układu równań
void printSystem(const vector<double> &a, const vector<double> &b,
                 const vector<double> &c, const vector<double> &d)
{
    int n = d.size();
    cout << "\nUkład równań:" << endl;
    for (int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            cout << setw(10) << b[i] << " * x[" << i << "] + "
                 << setw(10) << c[i] << " * x[" << i + 1 << "] = "
                 << setw(10) << d[i] << endl;
        }
        else if (i == n - 1)
        {
            cout << setw(10) << a[i] << " * x[" << i - 1 << "] + "
                 << setw(10) << b[i] << " * x[" << i << "] = "
                 << setw(10) << d[i] << endl;
        }
        else
        {
            cout << setw(10) << a[i] << " * x[" << i - 1 << "] + "
                 << setw(10) << b[i] << " * x[" << i << "] + "
                 << setw(10) << c[i] << " * x[" << i + 1 << "] = "
                 << setw(10) << d[i] << endl;
        }
    }
}

// Funkcja pomocnicza do wyświetlania rozwiązania
void printVector(const vector<double> &x)
{
    cout << "\nRozwiązanie:" << endl;
    for (int i = 0; i < x.size(); i++)
    {
        cout << "x[" << i << "] = " << setw(12) << setprecision(6)
             << fixed << x[i] << endl;
    }
}

int main()
{
    // Przykład 1: Układ 4x4
    cout << "=== PRZYKŁAD 1: Układ 4x4 ===" << endl;

    int n = 4;
    vector<double> a = {0.0, 1.0, 1.0, 1.0}; // dolna przekątna (a[0]=0)
    vector<double> b = {2.0, 2.0, 2.0, 2.0}; // główna przekątna
    vector<double> c = {1.0, 1.0, 1.0, 0.0}; // górna przekątna (c[n-1]=0)
    vector<double> d = {5.0, 6.0, 6.0, 5.0}; // prawa strona
    vector<double> x(n);                     // rozwiązanie

    printSystem(a, b, c, d);
    thomasMethod(a, b, c, d, x);
    printVector(x);

    // Przykład 2: Inny 5x5
    cout << "\n\n=== PRZYKŁAD 2: Układ 5x5 ===" << endl;

    n = 5;
    a = {0.0, 2.0, 2.0, 2.0, 2.0};
    b = {4.0, 4.0, 4.0, 4.0, 4.0};
    c = {1.0, 1.0, 1.0, 1.0, 0.0};
    d = {2.0, 6.0, 4.0, 9.0, 8.0};
    x.resize(n); // Funkcja obiektu vector, zmienia rozmiar do n

    printSystem(a, b, c, d);
    thomasMethod(a, b, c, d, x);
    printVector(x);

    // Przykład z Teams
    cout << "=== PRZYKŁAD Z PDF ===" << endl;

    n = 5;

    // Dane z macierzy (interpretacja jako trójdiagonalna)
    a = {0.0, 3.0, 1.0, 1.0, 2.0}; // dolna przekątna
    b = {2.0, 1.0, 2.0, 1.0, 2.0}; // główna przekątna
    c = {2.0, 1.0, 4.0, 1.0, 0.0}; // górna przekątna
    d = {2.0, 6.0, 4.0, 1.0, 4.0}; // prawa strona
    x.resize(n);

    printSystem(a, b, c, d);
    thomasMethod(a, b, c, d, x);
    printVector(x);

    return 0;
}
