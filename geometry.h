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
    Point & operator += (Point const & p) {
        x += p.x; y += p.y;
        return *this;
    }
    Point & operator -= (Point const & p) {
        x -= p.x; y -= p.y;
        return *this;
    }
    Point & operator - () {
        x = -x; y = -y;
        return *this;
    }
    Point & operator *= (double l) {
        x *= l; y *= l;
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
Point operator * (Point const & p, double l) {
    Point cpy(p);
    return cpy *= l;
}
Point operator * (double l, Point const & p) {
    Point cpy(p);
    return cpy *= l;
}
std::ostream & operator << (std::ostream & cout, Point const & p) {
    cout << p.x << ' ' << p.y << '\n';
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

class Line {
public:
    Line () = default;
    Line (Point const & p0, Point const & p1) :
            A(p0.y-p1.y), B(p1.x - p0.x), C(p0.x * p1.y - p0.y * p1.x) {}
    Line (double k, double b) : A(-k), B(1), C(-b) {}

    [[nodiscard]] bool belongs (Point const & p) const {
        return A * p.x + B * p.y + C == 0;
    }

    bool operator == (Line const & right) const {
        return Matrix<2, 3>{{A, B, C}, {right.A, right.B, right.C}}.rank() != 2;
    }

    bool operator != (Line const & right) const { return !(*this == right); }

    // линия хранится в виде Ax+By+C=0
    double A = -1;
    double B = 1;
    double C = 0;

    friend std::ostream & operator << (std::ostream & cout, Line const & l);
};

std::ostream & operator << (std::ostream & cout, Line const & l) {
    cout << l.A << ' ' << l.B << ' ' << l.C << '\n';
    return cout;
}

