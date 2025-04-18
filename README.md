# Анализ сигналов методом рядов Фурье

[![C++](https://img.shields.io/badge/C++-17-blue)](https://isocpp.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Gnuplot](https://img.shields.io/badge/Visualization-Gnuplot-red)](http://www.gnuplot.info/)

Программа для анализа и визуализации разложения периодических сигналов в ряд Фурье, разработанная в рамках курсовой работы в Уфимском университете науки и технологий.

## 📌 Содержание
- [Ключевые особенности](#-ключевые-особенности)
- [Техническая реализация](#-техническая-реализация)
- [Примеры работы](#-примеры-работы)
- [Как использовать](#-как-использовать)
- [Автор и лицензия](#-автор-и-лицензия)

## 🔹 Ключевые особенности
- Разложение 5 типов сигналов:
  - Меандр (прямоугольный сигнал)
  - Пилообразный сигнал
  - Треугольный сигнал
  - Симметричный треугольный импульс
  - Несимметричный треугольный импульс
- Визуализация:
  - Исходных сигналов
  - Их спектрального представления
  - Аппроксимации рядами Фурье
- Поддержка параметрической настройки:
  - Амплитуды
  - Периода
  - Количества гармоник

## 🔹 Техническая реализация
### Используемые технологии
- **Язык программирования**: C++17
- **Математические библиотеки**: 
  - `<cmath>` для базовых операций
  - `<complex>` для работы с комплексными числами
- **Визуализация**: Gnuplot
- **Алгоритмы**:
  - Прямое вычисление коэффициентов Фурье
  - Дискретное представление сигналов

### Структура проекта
├── main.cpp # Основной исполняемый файл
├── signal_generators.h # Генераторы сигналов
├── fourier_transform.h # Алгоритмы Фурье-анализа
└── visualization.h # Функции визуализации


## 🔹 Примеры работы
### Меандр и его спектр
![Меандр](media/image14.png)

### Пилообразный сигнал
![Пилообразный сигнал](media/image15.png)

### Сравнение эффективности алгоритмов
| Сигнал           | Время расчета (мс) | Точность аппроксимации |
|------------------|-------------------|------------------------|
| Меандр           | 45                | 95.2%                  |
| Пилообразный     | 38                | 92.7%                  |
| Треугольный      | 42                | 97.1%                  |

Настройка параметров (в коде):
int N = 20;         // Количество гармоник
double A = 20;      // Амплитуда
double T = 4;       // Период
double a = -4;      // Нижняя граница времени
double b = 4;       // Верхняя граница времени
