#include <iostream>
#include "./matrix.h"
#include <vector>
using std::vector;
#include <cmath>
#include <stdarg.h>

struct Point {
    double x = 0;
    double y = 0;
    Point () = default;
    Point (double x, double y) : x(x), y(y) {}
    Point (Point const & p) = default;
    Point & operator = (Point const & p) = default;
    Point & operator = (Matrix<1, 3> const & M) {
        x = M[0][0]; y = M[0][1]; return *this;
    }
    explicit operator Matrix <1, 3> () {
        return Matrix<1, 3> ({{x, y, 1}});
    }

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

    static double distance (Point const & p1, Point const & p2);
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

double Point::distance(const Point &p1, const Point &p2)  {
    Point diff = p1 - p2;
    return sqrt(diff.x * diff.x + diff.y * diff.y);
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

class Shape {
public:
    virtual ~Shape () = default;
    virtual void rotate (const Point &, double) = 0;
    virtual void reflect (const Point &) = 0;
    virtual void reflect (const Line &) = 0;
    virtual void scale (const Point &, double) = 0;
    virtual bool operator == (Shape const &) const = 0;
    virtual bool operator != (Shape const &) const = 0;
    [[nodiscard]] virtual double perimeter () const = 0;
    [[nodiscard]] virtual double area () const = 0;
    [[nodiscard]] virtual bool isCongruentTo (Shape const &) const = 0;
    [[nodiscard]] virtual bool isSimilarTo (Shape const &) const = 0;
    [[nodiscard]] virtual bool containsPoint (Shape const &) const = 0;
};

class Polygon : public Shape {
public:
    Polygon () = default;
    ~Polygon () override = default;

    Polygon (Point ver1, Point ver2, Point ver3, ...) : vertices({ver1, ver2, ver3}) {
        va_list points;
        va_start(points, ver3);

    }

    explicit Polygon (vector<Point> const & vertices) : vertices(vertices) {}
    [[nodiscard]] size_t verticesCount () const {
        return vertices.size();
    }
    [[nodiscard]] bool isConvex () const;
    void rotate (const Point & center, double angle) override {
        applyMatrix(Matrix<3, 3>{ { cos(angle), sin(angle), 0 },
                                      { -sin(angle), cos(angle), 0 },
                                      { -center.x * (cos(angle) - 1) + center.y * sin(angle),
                                        -center.x * sin(angle) - center.y * (cos(angle) - 1), 1 } } );
    }
    void reflect (const Point & center) override {

    }
    void reflect (const Line &) override {

    }
    void scale (const Point &, double) override {

    }
    bool operator == (Polygon const & s) const {
        if (vertices.size() != s.vertices.size()) return false;
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t j = i;
            /*for (; j < i + vertices.size(); ++j) {
                if (vertices)
            }*/
        }
    }
    bool operator != (Shape const &) const override {

    }
    [[nodiscard]] double perimeter () const override {
        Point prev = *--vertices.end();
        double res = 0;
        for (auto vertex : vertices) {
            res += Point::distance(vertex, prev);
            prev = vertex;
        }
        return res;
    }
    [[nodiscard]] double area () const override {
        double res = 0;
        for (size_t i = 0; i < vertices.size(); i += 2) {
            res += (vertices[i + 1] - vertices[i]) * (vertices[i + 2] - vertices[i]) * 0.5f;
        }
        return fabs(res);
    }
    [[nodiscard]] bool isCongruentTo (Shape const &) const override {

    }
    [[nodiscard]] bool isSimilarTo (Shape const &) const override {

    }
    [[nodiscard]] bool containsPoint (Shape const &) const override {

    }
private:
    vector <Point> vertices = {};
    void applyMatrix (Matrix <3, 3> const & M) {
        for (auto & vertex : vertices) {
            vertex = static_cast<Matrix<1, 3>>(vertex) * M;
        }
    }
};

class Rectangle : public Polygon {

};