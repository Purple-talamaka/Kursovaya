#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> 
#include <iterator> 
#include <fstream>
#include <complex>

using namespace std;
const double pi = 3.1415;
const double e = 2.71828;


int drawGraph(vector <double> t, vector <double> M, vector <double> S);
vector <double> create_meander(double amplitude, double period, vector <double> t);
vector <double> create_sawtooth_signal(double amplitude, double period, vector <double> t);
vector <double> create_triangle_signal(double amplitude, double period, vector <double> t);
vector <double> Fourier_series_meander(double amplitude, double period, int number_of_members, vector <double> t);
vector <double> Fourier_series_sawtooth_signal(double amplitude, double period, int number_of_members, vector <double> t);
vector <double> Fourier_series_triangle_signal(double amplitude, double period, int number_of_members, vector <double> t);

vector <double> symmetric_triangle_signal(double amplitude, double period, vector <double> t);
vector <double> symmetric_triangle_signal_spectrum(double amplitude, double period, vector <double> t);
vector <double> asymmetrical_triangle_signal(double amplitude, double period, vector <double> t);
vector <double> asymmetrical_triangle_signal_spectrum(double amplitude, double period, vector <double> t);
int drawGraph2(vector<double> t, vector<double> M, vector<double> S, double amplitude, double period);

int main() {
    int N = 20; // Количество членов ряда Фурье, участвующих в расчете
    double A = 20; // Амплитуда функций
    double T = 4; // Период функций
    double a = -4; // Нижняя граница графиков по x
    double b = 4; // Верхняя граница графиков по x
    double d = 0.001; // Шаг разбиения на точки

    vector <double> t;
    for (double i = a; i <= b; i += d) {
        t.push_back(i);
    }

    cout << "Введите цифру" << endl;
    cout << "1 - меандр" << endl;
    cout << "2 - пилообразный сигнал" << endl;
    cout << "3 - треугольный сигнал" << endl;
    cout << "4 - симметричный треугольный сигнал и его спектр" << endl;
    cout << "5 - несимметричный треугольный сигнал и его спектр" << endl;
    int selector;
    cin >> selector;

    vector <double> M;
    vector <double> S;
    switch (selector)
    {
    case 1:
        M = create_meander(A, T, t);
        S = Fourier_series_meander(A, T, N, t);
        drawGraph(t, M, S);
        break;
    case 2:
        M = create_sawtooth_signal(A, T, t);
        S = Fourier_series_sawtooth_signal(A, T, N, t);
        drawGraph(t, M, S);
        break;
    case 3:
        M = create_triangle_signal(A, T, t);
        S = Fourier_series_triangle_signal(A, T, N, t);
        drawGraph(t, M, S);
        break;
    case 4:
        M = symmetric_triangle_signal(A, T, t);
        S = symmetric_triangle_signal_spectrum(A, T, t);
        drawGraph2(t, M, S, A, T);
        break;
    case 5:
        M = asymmetrical_triangle_signal(A, T, t);
        S = asymmetrical_triangle_signal_spectrum(A, T, t);
        drawGraph2(t, M, S, A, T);
        break;
    default:
        cout << "Вы неправильно ввели цифру(((" << endl;
        break;
    }

    return 0;
}

// Возвращает массив значений графика меандра
vector <double> create_meander(double amplitude, double period, vector <double> t) {
    vector <double> M;

    for (int i = 0; i < t.size(); i++) {
        if (t[i] + period / 4 >= 0 && (int)((t[i] + period / 4) * 2 / period) % 2 == 0)
            M.push_back(amplitude);
        else if (t[i] + period / 4 < 0 && (int)((t[i] + period / 4) * 2 / period) % 2 != 0)
            M.push_back(amplitude);
        else
            M.push_back(0);
    }

    return M;
}

// Возвращает массив значений графика ряда Фурье для меандра
vector <double> Fourier_series_meander(double amplitude, double period, int number_of_members, vector <double> t) {
    vector <double> S;
    for (int i = 0; i < t.size(); i++) {
        S.push_back(amplitude / 2);
    }
    for (double n = 0; n < number_of_members; n++) {
        for (int i = 0; i < t.size(); i++) {
            S[i] += 2 * amplitude / (2 * (double)n + 1) / pi * pow(-1, (double)n) * cos(2 * pi * (2 * (double)n + 1) * t[i] / period);
        }
    }
    return S;
}

// Возвращает массив значений графика пилообразного сигнала
vector <double> create_sawtooth_signal(double amplitude, double period, vector <double> t) {
    vector <double> M;
    for (int i = 0; i < t.size(); i++) {
        int n = (t[i] + period / 2) / period;
        if ((t[i] + period / 2) >= 0)
            M.push_back(((t[i] + period / 2) - (double)n * period) * amplitude * 2 / period - amplitude);
        else
            M.push_back(((t[i] + period / 2) - (double)n * period) * amplitude * 2 / period + amplitude);
    }
    return M;
}

// Возвращает массив значений графика ряда Фурье для пилообразного сигнала
vector <double> Fourier_series_sawtooth_signal(double amplitude, double period, int number_of_members, vector <double> t) {
    vector <double> S;
    for (int i = 0; i < t.size(); i++) {
        S.push_back(0);
    }
    for (int i = 0; i < t.size(); i++)
        for (int n = 0; n < number_of_members; n++)
            S[i] += (2.0 * amplitude / pi) * pow(-1.0, (double)n) * sin(((double)n + 1.0) * 2.0 * pi / period * t[i]);
    return S;
}

// Возвращает массив значений графика треугольного сигнала
vector <double> create_triangle_signal(double amplitude, double period, vector <double> t) {
    vector <double> M;
    for (int i = 0; i < t.size(); i++) {
        int n = (t[i] + period / 2) / period * 2;
        if ((t[i] + period / 2) >= 0 && n % 2 == 0)
            M.push_back(((t[i] + period / 2) - (double)n * period / 2) * amplitude * 4 / period - amplitude);
        else if ((t[i] + period / 2) >= 0 && n % 2 != 0)
            M.push_back(((t[i] + period / 2) - (double)n * period / 2) * (-amplitude * 4 / period) + amplitude);
        else if ((t[i] + period / 2) < 0 && n % 2 == 0)
            M.push_back(((t[i] + period / 2) - (double)n * period / 2) * (-amplitude * 4 / period) - amplitude);
        else
            M.push_back(((t[i] + period / 2) - (double)n * period / 2) * amplitude * 4 / period + amplitude);
    }
    return M;
}

// Возвращает массив значений графика ряда Фурье для треугольного сигнала
vector <double> Fourier_series_triangle_signal(double amplitude, double period, int number_of_members, vector <double> t) {
    vector <double> S;
    for (int i = 0; i < t.size(); i++) {
        S.push_back(0);
    }
    for (int i = 0; i < t.size(); i++)
        for (int n = 0; n < number_of_members; n++)
            S[i] += (8.0 * amplitude / pow(pi, 2)) * 1 / ((double)n * 2.0 + 1.0) * cos(((double)n * 2.0 + 1.0) * 2.0 * pi / period * t[i]);
    return S;
}

// Возвращает массив значений графика симметричного треугольного импульса
vector <double> symmetric_triangle_signal(double amplitude, double period, vector <double> t) {
    vector <double> M;

    for (int i = 0; i < t.size(); i++) {
        if (abs(t[i]) <= period)
            M.push_back(amplitude * (1 - abs(t[i]) / period));
        else
            M.push_back(0);
    }

    return M;
}

// Возвращает массив значений графика спектра симметричного теугольного импульса
vector <double> symmetric_triangle_signal_spectrum(double amplitude, double period, vector <double> t) {
    vector <double> S;
    for (int i = 0; i < t.size(); i++)
        S.push_back(amplitude * period * pow(sin(t[i] * period / 2), 2) / pow(t[i] * period / 2, 2));
    return S;
}

// Возвращает массив значений графика несимметричного треугольного импульса
vector <double> asymmetrical_triangle_signal(double amplitude, double period, vector <double> t) {
    vector <double> M;

    for (int i = 0; i < t.size(); i++) {
        if (t[i] >= 0 && t[i] <= period)
            M.push_back(amplitude * t[i] / period);
        else
            M.push_back(0);
    }

    return M;
}

// Возвращает массив значений графика спектра несимметричного теугольного импульса
vector <double> asymmetrical_triangle_signal_spectrum(double amplitude, double period, vector <double> t) {
    vector <double> S;
    complex<double> j(0, 1);
    for (int i = 0; i < t.size(); i++)
        if (abs(t[i]) > 0.05)
            S.push_back(abs(amplitude / pow((j * t[i]), 2) / period - amplitude / pow((j * t[i]), 2) / period * pow(e, -j * t[i] * period) - amplitude / j / t[i] * pow(e, -j * t[i] * period)));
        else
            S.push_back(amplitude * period / 2);
    return S;
}

// Выводит два графика на одной координатной плоскости по заданным точкам
int drawGraph(vector <double> t, vector <double> M, vector <double> S) {
    // Запуск Gnuplot
    FILE* gnuplot = _popen("gnuplot -persistent", "w");
    if (gnuplot == nullptr) {
        std::cerr << "Ошибка при запуске Gnuplot." << std::endl;
        return 1;
    }

    // Установка параметров Gnuplot
    fprintf(gnuplot, "set title 'Fourier Series for Selected Signal'\n");
    fprintf(gnuplot, "set xlabel 't'\n");
    fprintf(gnuplot, "set ylabel 's(t)'\n");
    fprintf(gnuplot, "set xrange [%g:%g]\n", *min_element(t.begin(), t.end()), *max_element(t.begin(), t.end()));
    fprintf(gnuplot, "set yrange [%g:%g]\n", *min_element(S.begin(), S.end()) - 0.5,
        *max_element(S.begin(), S.end()) + 0.5);

    // Установка стиля линии и создание функции для отображения
    fprintf(gnuplot, "plot '-' with lines lc rgb '0x000000' lw 5 title 'original function', '-' with lines lc rgb '0xC00D0D' lw 2 title 'fourier series'\n");

    // Генерация точек для исходной функции
    for (int i = 0; i < t.size(); ++i) {
        fprintf(gnuplot, "%g %g\n", t[i], M[i]);
    }
    fprintf(gnuplot, "e\n");

    // Генерация точек для ряда Фурье
    for (int i = 0; i < t.size(); ++i) {
        fprintf(gnuplot, "%g %g\n", t[i], S[i]);
    }
    fprintf(gnuplot, "e\n");



    // Закрытие Gnuplot
    _pclose(gnuplot);
    return 0;
}

// Выводит два графика, каждый в своем окне
int drawGraph2(vector<double> t, vector<double> M, vector<double> S, double amplitude, double period) {
    // Первая фигура: исходная функция
    FILE* gnuplot1 = _popen("gnuplot -persistent", "w");
    if (gnuplot1 == nullptr) {
        std::cerr << "Ошибка при запуске Gnuplot (первая фигура)." << std::endl;
        return 1;
    }

    fprintf(gnuplot1, "set title 'Original Function'\n");
    fprintf(gnuplot1, "set xlabel 't'\n");
    fprintf(gnuplot1, "set ylabel 'S(t)'\n");
    fprintf(gnuplot1, "set xrange [%g:%g]\n", *min_element(t.begin(), t.end()), *max_element(t.begin(), t.end()));
    fprintf(gnuplot1, "set yrange [%g:%g]\n", *min_element(M.begin(), M.end()) - 0.5,
        *max_element(M.begin(), M.end()) + 0.5);
    fprintf(gnuplot1, "plot '-' with lines lc rgb '0x0000FF' lw 2 title 'Original Function'\n");
    for (int i = 0; i < t.size(); ++i) {
        fprintf(gnuplot1, "%g %g\n", t[i], M[i]);
    }
    fprintf(gnuplot1, "e\n");
    _pclose(gnuplot1);

    // Вторая фигура: ряд Фурье
    FILE* gnuplot2 = _popen("gnuplot -persistent", "w");
    if (gnuplot2 == nullptr) {
        std::cerr << "Ошибка при запуске Gnuplot (вторая фигура)." << std::endl;
        return 1;
    }

    fprintf(gnuplot2, "set title 'Specter'\n");
    fprintf(gnuplot2, "set xlabel 'w'\n");
    fprintf(gnuplot2, "set ylabel 'S(w)'\n");
    fprintf(gnuplot2, "set xrange [%g:%g]\n", *min_element(t.begin(), t.end()), *max_element(t.begin(), t.end()));
    // fprintf(gnuplot2, "set yrange [%g:%g]\n", *min_element(S.begin(), S.end()) - 0.5,
    //                                          *max_element(S.begin(), S.end()) + 0.5);
    fprintf(gnuplot2, "set yrange [%g:%g]\n", *min_element(S.begin(), S.end()) - 0.5, amplitude * period + 0.5);

    fprintf(gnuplot2, "plot '-' with lines lc rgb '0xFF0000' lw 2 title 'Specter'\n");
    for (int i = 0; i < t.size(); ++i) {
        fprintf(gnuplot2, "%g %g\n", t[i], S[i]);
    }
    fprintf(gnuplot2, "e\n");
    _pclose(gnuplot2);

    return 0;
}
