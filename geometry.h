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
    Point & operator /= (double l) {
        x /= l; y /= l;
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
Point operator / (Point const & p, double l) {
    Point cpy(p);
    return cpy /= l;
}
Point operator / (double l, Point const & p) {
    Point cpy(p);
    return cpy /= l;
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
    Line (double A, double B, double C) : A(A), B(B), C(C) {}
    Line (Point const & p0, Point const & p1) :
            A(p0.y-p1.y), B(p1.x - p0.x), C(p0.x * p1.y - p0.y * p1.x) {}
    Line (double k, double b) : A(-k), B(1), C(-b) {}
    Line (Point const & p0, double k) : A(-k), B(1), C(k * p0.x - p0.y) {}

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
    static Point intersection (Line const & l1, Line const & l2) {
        double det = l1.A * l2.B - l1.B * l2.A;
        if (det == 0) return {1 / 0.0, 1 / 0.0};
        return { (l1.C * l2.B - l2.C * l1.B) / det, (l1.C * l2.A - l2.C * l1.A) / det };
    }
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
    virtual void reflect(const Line &) = 0;
    virtual void scale (const Point &, double) = 0;
    virtual bool operator == (Shape const &) const {
        return false;
    }
    virtual bool operator != (Shape const &) const {
        return true;
    }
    [[nodiscard]] virtual double perimeter () const = 0;
    [[nodiscard]] virtual double area () const = 0;
    [[nodiscard]] virtual bool isSimilarTo (Shape const &) const {
        return false;
    }
    [[nodiscard]] virtual bool isCongruentTo (Shape const &) const {
        return false;
    }
    [[nodiscard]] virtual bool containsPoint (Point const &) const {
        return false;
    }
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
    [[nodiscard]] bool isConvex () const {
        return true;
    }
    void rotate (const Point & center, double angle) override {
        applyMatrix(Matrix<3, 3>{ { cos(angle), sin(angle), 0 },
                                      { -sin(angle), cos(angle), 0 },
                                      { -center.x * (cos(angle) - 1) + center.y * sin(angle),
                                        -center.x * sin(angle) - center.y * (cos(angle) - 1), 1 } } );
    }
    void reflect (const Point & center) override {
        applyMatrix(Matrix<3, 3>{ { -1, 0, 0 },
                                      { 0, -1, 0 },
                                      { 2 * center.x, 2 * center.y, 1 } });
    }
    void reflect(const Line & axis) override {
        double k = axis.A * axis.A + axis.B * axis.B;
        applyMatrix(Matrix<3, 3>{
            { (axis.B * axis.B - axis.A * axis.A) / k, -2 * axis.A * axis.B / k, 0 },
            { -2 * axis.A * axis.B / k, (axis.A * axis.A - axis.B * axis.B) / k, 0 },
            { -2 * axis.C * axis.A / k, -2 * axis.C * axis.B / k, 1 } });
    }
    void scale (const Point & center, double coefficient) override {
        applyMatrix(Matrix<3, 3>{
                { coefficient, 0, 0 },
                { 0, coefficient, 0 },
                { (-coefficient + 1) * center.x, (-coefficient + 1) * center.y, 1 }
        });
    }
    bool operator == (Polygon const & right) const {
        if (vertices.size() != right.vertices.size()) return false;
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t j = i;
            for (; j < i + vertices.size(); ++j) {
                if (vertices[j % vertices.size()] != right.vertices[j - i])
                    break;
            }
            if (j == i + vertices.size())
                return true;
        }
        return false;
    }
    bool operator != (Polygon const & right) const {
        return !(*this == right);
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
    [[nodiscard]] bool isSimilarTo (Polygon const & right) const {
        if (vertices.size() != right.vertices.size()) return false;
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t j = i;
            Point diff = vertices[j % vertices.size()] - right.vertices[j - i];
            for (++j; j < i + vertices.size(); ++j) {
                if (vertices[j % vertices.size()] - right.vertices[j - i] != diff)
                    break;
            }
            if (j == i + vertices.size())
                return true;
        }
        return false;
    }
    [[nodiscard]] bool isCongruentTo (Polygon const & right) const {
        if (vertices.size() != right.vertices.size()) return false;
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t j = i;
            Point center = Line::intersection(
                    Line(vertices[j % vertices.size()], right.vertices[j - i]),
                    Line(vertices[(j + 1) % vertices.size()], right.vertices[j - i + 1]));
            if (center == Point{1/0.0, 1/0.0}) break;
            for (++j; j < i + vertices.size(); ++j) {
                if (!Line(vertices[j % vertices.size()], right.vertices[j - i]).belongs(center))
                    break;
            }
            if (j == i + vertices.size())
                return true;
        }
        return false;
    }
    [[nodiscard]] bool containsPoint (Point const & p) const override {
        double sign = (vertices[vertices.size() - 1] - vertices[0]) * p;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if ((vertices[i] - vertices[i + 1]) * p * sign < 0)
                return false;
        }
        return true;
    }

    friend std::ostream & operator << (std::ostream & cout, Polygon const & p);
protected:
    vector <Point> vertices = {};
    void applyMatrix (Matrix <3, 3> const & M) {
        for (auto & vertex : vertices) {
            vertex = Matrix<1, 3>{{vertex.x, vertex.y, 1}} * M;
        }
    }
};

std::ostream & operator << (std::ostream & cout, Polygon const & p) {
    for (auto x : p.vertices)
        cout << x << ' ';
    cout << '\n';
    return cout;
}

class Rectangle : public Polygon {
public:
    Rectangle (Point ver1, Point ver2, double coefficient) {

    }
    Point center () const {
        return vertices[0] - vertices[2];
    }
    std::pair <Line, Line> diagonals () const {
        return {{vertices[0], vertices[2]}, {vertices[1], vertices[3]}};
    }

};

class Triangle : public Polygon {
public:
    Triangle (vector<Point> vertices) : Polygon (vertices) {}
    [[nodiscard]] bool isCongruentTo (Triangle const & right) const {
        const Matrix<6, 1> raw_tr = Matrix<6, 6>{{vertices[0].x, vertices[0].y, 1, 0, 0, 0},
                     {0, 0, 0, vertices[0].x, vertices[0].y, 1},
                     {vertices[1].x, vertices[1].y, 1, 0, 0, 0},
                     {0, 0, 0, vertices[1].x, vertices[1].y, 1},
                     {vertices[2].x, vertices[2].y, 1, 0, 0, 0},
                     {0, 0, 0, vertices[2].x, vertices[2].y, 1}}.inverted()
                     * Matrix<6, 1>{{right.vertices[0].x},
                                    {right.vertices[0].y},
                                    {right.vertices[1].x},
                                    {right.vertices[1].y},
                                    {right.vertices[2].x},
                                    {right.vertices[2].y}};
        const Matrix<2, 3> tr = {{raw_tr[0][0], raw_tr[1][0], raw_tr[2][0]},
                           {raw_tr[3][0], raw_tr[4][0], raw_tr[5][0]}};
        const Point center = (vertices[0] + vertices[1] + vertices[2]) / 3;
        const Point right_center = (right.vertices[0] + right.vertices[1] + right.vertices[2]) / 3;
        std::cout << tr << center << right_center;
        std::cout << ((tr * Matrix<3, 1>{{center.x}, {center.y}, {1}})[1][0]);

        const Matrix<2, 1> res = tr * Matrix<3, 1>{{center.x}, {center.y}, {1}}
                        - Matrix<2, 1>{{right_center.x}, {right_center.y}};
        const double eps = 1e-6;
        return res[0][0] < eps && res[1][0] < eps;
    }
};

class Ellipse : public Shape {
public:
    ~Ellipse () = default;
    Ellipse (Point focus1, Point focus2, double axis) : focus1(focus1), focus2(focus2), axis(axis) {}
    void rotate (const Point &, double) {
        return;
    }
    void reflect (const Point & center) {
        applyMatrix(Matrix<3, 3>{ { -1, 0, 0 },
                      { 0, -1, 0 },
                      { 2 * center.x, 2 * center.y, 1 } });
    }
    void reflect (const Line & a) {
        double k = a.A * a.A + a.B * a.B;
        applyMatrix(Matrix<3, 3>{
                {(a.B * a.B - a.A * a.A) / k, -2 * a.A * a.B / k,          0 },
                {-2 * a.A * a.B / k,          (a.A * a.A - a.B * a.B) / k, 0 },
                {-2 * a.C * a.A / k,          -2 * a.C * a.B / k,          1 } });
    }
    void scale (const Point & center, double coefficient) {
        applyMatrix(Matrix<3, 3>{
                { coefficient, 0, 0 },
                { 0, coefficient, 0 },
                { (-coefficient + 1) * center.x, (-coefficient + 1) * center.y, 1 }
        });
        axis *= coefficient;
    }
    bool operator == (Ellipse const & right) const {
        return focus1 == right.focus1 && focus2 == right.focus2 && axis == right.axis;
    }
    bool operator != (Ellipse const & right) const {
        return !(*this == right);
    }
    [[nodiscard]] double perimeter () const {
        double b = saxis();
        return M_PI * (3 * (axis + b) - sqrt((3 * axis + b) * (axis + 3 * b)));
    }
    [[nodiscard]] double area () const {
        return M_PI * axis * saxis();
    }
    [[nodiscard]] bool isSimilarTo (Ellipse const & right) const {
        return axis == right.axis && saxis() == right.saxis();
    }
    [[nodiscard]] bool isCongruentTo (Ellipse const &) const {
        return false;
    }
    [[nodiscard]] bool containsPoint (Point const & p) const {
        Point center = (focus1 + focus2) / 2;
        double b = saxis();
        return (p.x - center.x) * (p.x - center.x) / axis / axis +
                (p.y - center.y) * (p.y - center.y) / b / b <= 1;
    }

    std::pair<Point, Point> focuses () const {
        return {focus1, focus2 };
    }
    std::pair<Line, Line> directrices () const {
        double c = Point::distance(focus1, focus2) / 2;
        Line a(focus1, focus2);
        Line b((focus2 + focus1) / 2, -a.A / a.B);
        return {Line{b.A, b.B, b.C + c}, Line{b.A, b.B, b.C - c}};
    }
protected:
    [[nodiscard]] double saxis () const {
        double c = Point::distance(focus1, focus2) / 2;
        return sqrt(axis * axis - c * c);
    }
    void applyMatrix (Matrix<3, 3> const & M) {
        auto f1 = Matrix<1, 3>{{focus1.x, focus1.y, 1}} * M;
        focus1.x = f1[0][0]; focus1.y = f1[0][1];
        auto f2 = Matrix<1, 3>{{focus2.x, focus2.y, 1}} * M;
        focus2.x = f2[0][0]; focus2.y = f2[0][1];
    }
    Point focus1, focus2;
    double axis; // большая полуось
};

class Circle : public Ellipse {
public:
    Circle (Point Center, double radius) : Ellipse(Center, Center, radius) {}
    double radius () const {
        return axis;
    }
};

class Square : public Rectangle {
public:
    Square (Point ver1, Point ver2) : Rectangle(ver1, ver2, 1) {}
    Circle circumscribedCircle () const {
        return Circle{center(), Point::distance(vertices[0], vertices[2]) / 2};
    }
    Circle inscribedCircle () const {
        return Circle{center(), Point::distance(vertices[0], vertices[1]) / 2};
    }
};