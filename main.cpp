/*
    Zaawansowane Metody Obliczeniowe
    Demon Creutza na regularnej sieci dwuwymiarowej (model Isinga 2D)

    Katarzyna Stańczyk
    06.11.2025
*/

#include <iostream>
#include <filesystem> //c++ 17
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <string>
#include <limits>

using namespace std;

// stałe
const int J = 1;                       //sprzężenie między spinami
const double kB = 1.0;                 //stała Boltzmanna

// Parametry histogramu i analizy
const int binWidth = 4;                //szerokość przedziału histogramu
const int minCount = 10;

// Parametry symulacji Monte Carlo
const double MCfraction = 0.20;        //część iteracji przeznaczona na równoważenie
const int E_demon_max = 1000;
const int histogramStep = 100;         //co ile jednostek energii zapisywać szczegółowe dane

//struktura do wyników
struct ResultPoint {
    double T;
    double m;     // <m>
    double slope;
    int E_demon;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    namespace fs = std::filesystem;

    try {
        fs::create_directories("histogram");
        fs::create_directories("magnetization");
    } catch (const fs::filesystem_error& e) {
        cerr << "Błąd przy tworzeniu folderów: " << e.what() << "\n";
    }


    int X, Y, numSweeps, numEnergies;
    if (!(cin >> X >> Y >> numSweeps >> numEnergies)) {
        cerr << "Błąd: nie udało się wczytać parametrów z wejścia.\n";
        return 1;
    }

    vector<int> initialDemonEnergies(numEnergies);
    for (int i = 0; i < numEnergies; i++) {
        if (!(cin >> initialDemonEnergies[i])) {
            cerr << "Błąd: nie udało się wczytać energii demona.\n";
            return 1;
        }
    }

    // Generator losowy
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> randX(0, X - 1);
    uniform_int_distribution<> randY(0, Y - 1);

    ofstream mT_file("mT.txt");
    mT_file << "E_demon;T;<m>;slope\n";

    for (int Ed0 : initialDemonEnergies) {
        int L = X * Y;
        int demonEnergy = max(0, min(E_demon_max, Ed0));

        // Ferromagnetyk wszystkie spiny +1
        vector<vector<int>> spin(X, vector<int>(Y, +1));
        vector<double> magnetization(numSweeps, 0.0);
        map<int, long long> hist;

        long long accept_count = 0;

        // Parametry równoważenia
        int numEquil = int(numSweeps * MCfraction);
        if (numEquil >= numSweeps) numEquil = numSweeps / 5;
        if (numEquil < 1) numEquil = 1;

        // petla Monte Carlo
        for (int sweep = 0; sweep < numSweeps; ++sweep) {
            for (int attempt = 0; attempt < L; ++attempt) {
                int i = randX(gen);
                int j = randY(gen);

                int up = spin[(i + X - 1) % X][j];
                int down = spin[(i + 1) % X][j];
                int left = spin[i][(j + Y - 1) % Y];
                int right = spin[i][(j + 1) % Y];
                int sumNeighbours = up + down + left + right;

                int deltaE = 2 * spin[i][j] * sumNeighbours*J;

                bool flipped = false;
                if (deltaE <= 0) {
                    spin[i][j] = -spin[i][j];
                    demonEnergy -= deltaE;
                    if (demonEnergy > E_demon_max) demonEnergy = E_demon_max;
                    flipped = true;
                } else if (demonEnergy >= deltaE) {
                    spin[i][j] = -spin[i][j];
                    demonEnergy -= deltaE;
                    if (demonEnergy < 0) demonEnergy = 0;
                    flipped = true;
                }

                if (flipped) accept_count++;
            }


            if (sweep >= numEquil) {
                int binE = (demonEnergy / binWidth) * binWidth;
                hist[binE]++;
            }


            long long sumSpin = 0;
            for (int x = 0; x < X; ++x)
                for (int y = 0; y < Y; ++y)
                    sumSpin += spin[x][y];
            magnetization[sweep] = double(sumSpin) / double(L);
        }


        vector<double> xs, ys, ws;
        for (auto &p : hist) {
            if (p.second >= minCount) {
                xs.push_back((double)p.first);
                ys.push_back(log((double)p.second));
                ws.push_back((double)p.second);
            }
        }
        double a = 0.0, b = 0.0, T = 0.0;
        size_t n = xs.size();
        if (n >= 2) {
            double Sw = 0, Swx = 0, Swy = 0, Swxx = 0, Swxy = 0;
            for (size_t k = 0; k < n; ++k) {
                double w = ws[k];
                Sw += w;
                Swx += w * xs[k];
                Swy += w * ys[k];
                Swxx += w * xs[k] * xs[k];
                Swxy += w * xs[k] * ys[k];
            }
            double denom = Sw * Swxx - Swx * Swx;
            if (fabs(denom) > 1e-12) {
                a = (Sw * Swxy - Swx * Swy) / denom;
                b = (Swxx * Swy - Swx * Swxy) / denom;
                T = -1.0 / a;
            } else {
                T = numeric_limits<double>::quiet_NaN();
            }
        }

        //Średnia magnetyzacja
        double sumMag = 0.0;
        for (int i = numEquil; i < numSweeps; ++i)
            sumMag += magnetization[i];
        double avgMag = sumMag / (numSweeps - numEquil);

        if (Ed0 % histogramStep == 0) {

            ofstream d_en_hist_file("histogram/histogram_E=" + to_string(Ed0) + ".txt");
            d_en_hist_file << "E_demon;N(E);ln(N)\n";
            for (auto &p : hist) {
                if (p.second > 0)
                    d_en_hist_file << p.first << ";" << p.second << ";" << log((double)p.second) << "\n";
            }
            d_en_hist_file.close();


            ofstream sweep_mag_file("magnetization/magnetization_E=" + to_string(Ed0) + ".txt");
            sweep_mag_file << "sweep;<m>;\n";
            for (int i = 0; i < numSweeps; ++i)
                sweep_mag_file << (i + 1) << ";" << magnetization[i] << "\n";
            sweep_mag_file.close();
        }

        // Wyniki zbiorcze mT.txt
        mT_file << Ed0 << ";" << T << ";" << avgMag << ";" << a << "\n";


        cout << "E_demon=" << Ed0 << "  <m>=" << avgMag << "  slope=" << a << "  T=" << T;
        if (Ed0 % histogramStep == 0)
            cout << "  -> zapisano histogram i magnetyzacje";
        cout << endl;
    }

    mT_file.close();
    cout << "Zapisano pliki w katalogach ./magnetization/ i ./histogram/ oraz mT.txt\n";
    return 0;
}

