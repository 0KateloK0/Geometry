#include <iostream>
#include <vector>
#include <cassert>
using std::vector;

template <size_t N, size_t M, typename Field = double>
class Matrix {
public:
    explicit Matrix (vector<vector<Field>> const & v) {
        assert(v.size() == N);
        destroyBody();
        body = new Field*[N];
        for (size_t i = 0; i < N; ++i) {
            assert(v[i].size() == M);
            body[i] = new Field[M];
            for (size_t j = 0; j < M; ++j)
                body[i][j] = v[i][j];
        }
    }

    explicit Matrix (vector<vector<int>> const & v) {
        assert(v.size() == N);
        destroyBody();
        body = new Field*[N];
        for (size_t i = 0; i < N; ++i) {
            assert(v[i].size() == M);
            body[i] = new Field[M];
            for (size_t j = 0; j < M; ++j)
                body[i][j] = static_cast<Field>(v[i][j]);
        }
    }

    Matrix (std::initializer_list<std::initializer_list<double>> const & v) {
        assert(v.size() == N);
        destroyBody();
        body = new Field*[N];
        size_t i = 0;
        for (auto & x : v) {
            assert(x.size() == M);
            body[i] = new Field[M];
            size_t j = 0;
            for (auto & y : x) {
                body[i][j] = static_cast<Field>(y);
                ++j;
            }
            ++i;
        }
    }

    Matrix (Matrix const & X) {
        body = new Field*[N];
        for (size_t i = 0; i < N; ++i) {
            body[i] = new Field[M];
            for (size_t j = 0; j < M; ++j)
                body[i][j] = X.body[i][j];
        }
    }

    Matrix () {
        createBody();
    }

    template <size_t A = N, size_t B = M, typename = std::enable_if<A == B>>
    Matrix () {
        createBody();
        for (size_t i = 0; i < N; ++i)
            body[i][i] = 1;
    }

    ~Matrix () {
        destroyBody();
    }

    Matrix & operator = (Matrix const & right) {
        if (&right == this) return *this;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = right.body[i][j];
            }
        }
        return *this;
    }
    
    Matrix & operator += (Matrix const & right) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] += right.body[i][j];
            }
        }
        return *this;
    }

    Matrix & operator -= (Matrix const & right) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] -= right.body[i][j];
            }
        }
        return *this;
    }

    template <size_t K>
    Matrix<N, K, Field> operator * (Matrix <M, K, Field> const & right) const {
        Matrix<N, K, Field> res;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < K; ++j) {
                res.body[i][j] = static_cast<Field>(0);
                for (size_t m = 0; m < M; ++m) {
                    res.body[i][j] += body[i][m] * right.body[m][j];
                }
            }
        }
        return res;
    }

    Matrix & operator *= (Field const & right) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] *= right;
            }
        }
        return *this;
    }

    Matrix & operator *= (Matrix const & right) {
        static_assert(N == M);
        return *this = *this * right;
    }

    Field det () const {
        static_assert(N == M);
        if (N == 2) {
            return body[0][0] * body[1][1] - body[0][1] * body[1][0];
        }
        const Matrix<N, M, Field> ret(makeEchelonForm(*this));
        auto res = static_cast<Field>(1);
        for (size_t i = 0; i < N; ++i) {
            res *= ret.body[i][i];
        }
        return res;
    }
    
    Matrix<M, N, Field> transposed () const {
        Matrix <M, N, Field> ret;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                ret.body[j][i] = body[i][j];
            }
        }
        return ret;
    }

    [[nodiscard]] size_t rank () const {
        const auto X = echelon_<N, M, Field>::f(*this);
        size_t count = 0;
        for (size_t i = 0; i < std::min(N, M); ++i) {
            if (X.body[i][i] == static_cast<Field>(0)) {
                size_t j = i + 1;
                while (j < std::max(N, M) && X.body[i][j] == static_cast<Field>(0)) ++j;
                if (j == std::max(N, M)) return count;
                else ++count;
            }
            else ++count;
        }
        return count;
    }

    Field trace () const {
        static_assert(N == M);
        Field ret{0};
        for (size_t i = 0; i < N; ++i) {
            ret += body[i][i];
        }
        return ret;
    }

    Matrix inverted () const {
        Matrix<N, M, Field> ret = *this;
        ret.invert();
        return ret;
    }

    Matrix & invert () {
        // Используется метод Гаусса. Изначально создается матрица, рядом с которой записывается единичная матрица,
        // для которой повторяются все операции
        static_assert(N == M);
        Matrix <N, 2 * M, Field> P;

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                P.body[i][j] = body[i][j];
            }
            for (size_t j = M; j < 2 * M; ++j) {
                P.body[i][j] = static_cast<Field>(i == j - M);
            }
        }

        P = makeEchelonForm(P, true);

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                body[i][j] = P.body[i][j + M];
            }
        }

        return *this;
    }
    
    vector<Field> getRow (size_t ind) const {
        return vector<Field>{body[ind], body[ind] + M};
    }

    vector<Field> getColumn (size_t ind) const {
        auto ret = vector<Field>(N);
        for (size_t i = 0; i < N; ++i)
            ret[i] = body[i][ind];
        return ret;
    }

    Field* operator [](const size_t ind) {
        return body[ind];
    }

    vector<Field> operator [](const size_t ind) const {
        return vector<Field>{body[ind], body[ind] + M};
    }

    template <size_t A, size_t B, typename T>
    bool operator == (Matrix<A, B, T> const & right) const {
        if (A != N || B != M) return false;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j)
                if (body[i][j] != right.body[i][j]) return false;
        }
        return true;
    }

    template <size_t A, size_t B, typename T>
    bool operator != (Matrix<A, B, T> const & right) const {
        return !(*this == right);
    }

    // для того чтобы можно было использовать приватные поля матриц любого размера
    template <size_t A, size_t B, typename T>
    friend class Matrix;
    
private:
    // не сделал через вектор векторов, потому что как таковой функционал вектора не используется
    Field** body = nullptr;

    void createBody () {
        destroyBody();
        body = new Field*[N];
        for (size_t i = 0; i < N; ++i) {
            body[i] = new Field[M];
        }
    }

    void destroyBody () {
        if (body != nullptr) {
            for (size_t i = 0; i < N; ++i) {
                delete [] body[i];
            }
            delete [] body;
        }
    }

    // эти структуры необходимы, чтобы можно было выполнять приведение к треугольному виду матрицы любого размера
    // использование: echelon_< размеры матрицы и тип элемента, 2 >::f(матрица)
    // без них попытки сделать похожий функционал приводили к CE.
    template <size_t A, size_t B, typename T, int v = 2>
    struct echelon_ {
        static Matrix<std::min(A, B), std::max(A, B), T> f (Matrix<A, B, T> const & X) {
            return echelon_<std::min(A, B), std::max(A, B), T, (A > B)>::f(X);
        }
    };

    template <size_t A, size_t B, typename T>
    struct echelon_ <A, B, T, 0> {
        static Matrix<A, B, T> f (Matrix<A, B, T> const & X) {
            return makeEchelonForm(X);
        }
    };

    template <size_t A, size_t B, typename T>
    struct echelon_ <A, B, T, 1> {
        static Matrix<A, B, T> f (Matrix <B, A, T> const & X) {
            return makeEchelonForm(X.transposed());
        }
    };

    // Предполагает что дана прямоугольная матрица размеров N*M, где M >= N.
    // матрица передается не по ссылке потому что в любом случае создавалась бы копия
    // преобразовывает матрицу X к ступенчатому виду методом Гаусса.
    // Если определитель не 0, может также выполнить обратный ход
    template <size_t A, size_t B, typename T>
    static Matrix<A, B, T> makeEchelonForm (Matrix<A, B, T> X, bool rev = false) {
        static_assert(B >= A);
        T zero = T(0); // чтобы не создавать триллион раз нули статик кастами
        size_t j = 0;
        for (size_t i = 0; i + 1 < A; ++i) {
            if (X.body[i][j] == zero) {
                //ищет первую строку в которой ненулевой body[k][j] число и складывает ее с i строкой
                // после чего продолжает метод Гаусса
                size_t k = i + 1;
                while (k < A && X.body[k][j] == zero) ++k;
                if (k == A) {
                    j += j + 1 < B;
                    --i;
                    continue;
                }
                for (size_t w = j; w < B; ++w) {
                    X.body[i][w] += X.body[k][w];
                }
                --i; // цикл выполнится заново, зайдя во вторую ветку ифа. Должно корректно работать и с i == 0;
            }
            else {
                for (size_t k = i + 1; k < A; ++k) {
                    const T q = X.body[k][j] / X.body[i][j];
                    for (size_t w = j; w < B; ++w) {
                        X.body[k][w] -= X.body[i][w] * q;
                    }
                }
                j += j + 1 < B;
            }
        }
        // выполняет обратный метод Гаусса, если rev = true
        if (rev) { // работает только если определитель больше нуля, иначе залезет в бесконечный цикл
            for (size_t i = A; i > 0; --i) {
                for (size_t k = i - 1; k > 0; --k) {
                    const T q = X.body[k - 1][i - 1] / X.body[i - 1][i - 1];
                    for (size_t w = i - 1; w < B; ++w) {
                        X.body[k - 1][w] -= X.body[i - 1][w] * q;
                    }
                }
                // приводит к единичному виду
                const T q = X.body[i - 1][i - 1];
                for (size_t k = 0; k < B; ++k) {
                    X.body[i - 1][k] /= q;
                }
            }
        }
        return X;
    }
};

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator + (Matrix<N, M, Field> const & left, Matrix<N, M, Field> const & right) {
    Matrix<N, M, Field> cpy = left;
    return cpy += right;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator - (Matrix<N, M, Field> const & left, Matrix<N, M, Field> const & right) {
    Matrix<N, M, Field> cpy = left;
    return cpy -= right;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator * (Matrix <N, M, Field> const & left, Field const & right) {
    Matrix<N, M, Field> cpy = left;
    return cpy *= right;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator * (Field const & left, Matrix <N, M, Field> const & right) {
    Matrix<N, M, Field> cpy = right;
    return cpy *= left;
}

template <size_t N, typename Field = double>
using SquareMatrix = Matrix <N, N, Field>;
