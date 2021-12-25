#include <iostream>
#include "./matrix.h"
#include <vector>
using std::vector;
#include <cmath>

struct Point {
    double x = 0;
    double y = 0;
    Point () = default;
    Point (double x, double y) : x(x), y(y) {}
    Point (Point const & p) = default;

    bool operator == (Point const & p) const {
        return x == p.x && y == p.y;
    }
    bool operator != (Point const & p) const {
        return !(*this == p);
    }

    // Точка по сути своей это вектор с началом (0, 0) и координатами (x, y). Посему точка будет иметь все функции
    // 2д вектора, а также возможность нахождения векторного произведения, которое будет возвращать ЗНАКОВОЕ значение
    // модуля результирующего вектора
    Point operator += (Point const & p) {
        x += p.x; y += p.y;
        return *this;
    }
    Point operator -= (Point const & p) {
        x -= p.x; y -= p.y;
        return *this;
    }
    Point operator - () {
        x = -x; y = -y;
        return *this;
    }

    friend std::ostream & operator << (std::ostream & cout, Point const & p);
};

Point operator + (Point const & p1, Point const & p2) {
    Point cpy = p1;
    return cpy += p2;
}
Point operator - (Point const & p1, Point const & p2) {
    Point cpy = p1;
    return cpy -= p2;
}
Point operator - (Point const & p) {
    Point cpy = p;
    return -cpy;
}
std::ostream & operator << (std::ostream & cout, Point const & p) {
    cout << p.x << ' ' << p.y;
    return cout;
}

// Векторное произведение. Возвращает знаковую площадь получившегося параллелограмма
double operator * (Point const & p1, Point const & p2) {
    return p1.x * p2.y - p1.y * p2.x;
}

// скалярное произведение
double operator ^ (Point const & p1, Point const & p2) {
    return p1.x * p2.x + p1.y * p2.y;
}