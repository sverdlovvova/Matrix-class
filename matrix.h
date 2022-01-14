#include <vector>
#include <cmath>
#include "biginteger.h"
#include "rational.h"
#include "residue.h"

template<size_t N, size_t M, typename Field = Rational>
class Matrix {
private:
    std::vector<std::vector<Field>> matrix;

public:
    Matrix() : matrix(N, std::vector<Field>(M, Field(0))) {
        if (N != M) {
            return;
        }
        for (size_t i = 0; i < N; ++i) {
            matrix[i][i] = Field(1);
        }
    }

    Matrix(const std::vector<std::vector<Field>> &values) : matrix(values) {}

    Matrix(const std::initializer_list<std::initializer_list<int>> &values) {
        matrix.resize(N);
        auto row = values.begin();
        for (size_t i = 0; i < N; ++i) {
            auto elem = (*row).begin();
            for (size_t j = 0; j < M; ++j) {
                matrix[i].push_back(Field(*elem));
                ++elem;
            }
            ++row;
        }
    }

    template<size_t K, size_t L>
    bool operator==(const Matrix<K, L, Field> &another) const {
        return matrix == another.matrix;
    }

    template<size_t K, size_t L>
    bool operator!=(const Matrix<K, L, Field> &another) const {
        return !(*this == another);
    }

    template<size_t K, size_t L>
    Matrix &operator+=(const Matrix<K, L, Field> &another) {
        static_assert(N == K && M == L);
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> rowAnother = another[i];
            for (size_t j = 0; j < M; ++j) {
                matrix[i][j] += rowAnother[j];
            }
        }
        return *this;
    }

    Matrix &operator*=(const Field &multiplier) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                matrix[i][j] *= multiplier;
            }
        }
        return *this;
    }

    template<size_t K, size_t L>
    Matrix &operator-=(const Matrix<K, L, Field> &another) {
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> rowAnother = another[i];
            for (size_t j = 0; j < M; ++j) {
                matrix[i][j] -= rowAnother[j];
            }
        }
        return *this;
    }

    template<size_t K, size_t L>
    Matrix operator+(const Matrix<K, L, Field> &another) const {
        Matrix result = *this;
        result += another;
        return result;
    }

    template<size_t K, size_t L>
    Matrix operator-(const Matrix<K, L, Field> &another) const {
        Matrix result = *this;
        result -= another;
        return result;
    }

    template<size_t K, size_t L>
    Matrix operator*(const Field &multiplier) const {
        Matrix result = *this;
        result *= multiplier;
        return result;
    }

    template<size_t K, size_t L>
    Matrix<N, L, Field> operator*(const Matrix<K, L, Field> &another) const {
        static_assert(M == K);
        Matrix<N, L, Field> result;
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> row(L, Field(0));
            for (size_t j = 0; j < L; ++j) {
                std::vector<Field> anotherColumn = another.getColumn(j);
                for (size_t k = 0; k < M; ++k) {
                    row[j] += matrix[i][k] * anotherColumn[k];
                }
            }
            result[i] = row;
        }
        return result;
    }

    template<size_t K, size_t L>
    Matrix<N, L, Field> &operator*=(const Matrix<K, L, Field> &another) {
        static_assert(M == K && N == M);
        *this = *this * another;
        return *this;
    }

    std::vector<Field> &operator[](int i) {
        return matrix[i];
    }

    std::vector<Field> operator[](int i) const {
        return matrix[i];
    }

    Matrix gauss(bool forInverting = false) const {
        Matrix copy = *this;
        size_t k = 0;
        size_t untilColumn = M;
        if (forInverting) {
            untilColumn /= 2;
        }
        for (size_t i = 0; i < untilColumn; ++i) {
            size_t j = k;
            while (j < N && copy.matrix[j][i] == Field(0)) {
                ++j;
            }
            if (j == N) {
                continue;
            }
            swap(copy.matrix[j], copy.matrix[k]);
            for (size_t t = k + 1; t < N; ++t) {
                Field coefficient = copy.matrix[t][i] / copy.matrix[k][i];
                for (size_t s = i; s < M; ++s) {
                    copy.matrix[t][s] -= copy.matrix[i][s] * coefficient;
                }
            }
            ++k;
        }
        return copy;
    }

    Matrix invertedGauss() const {
        Matrix result = *this;
        for (size_t c = M / 2 - 1; c + 1 != 0; --c) {
            for (size_t i = 0; i < c; ++i) {
                Field coefficient = result.matrix[i][c] / result.matrix[c][c];
                for (size_t j = c; j < M; ++j) {
                    result.matrix[i][j] -= coefficient * result.matrix[c][j];
                }
            }
        }
        return result;
    }

    Field det() const {
        static_assert(N == M);
        Matrix copy = gauss();
        Field det = Field(1);
        for (size_t i = 0; i < N; ++i) {
            det *= copy.matrix[i][i];
        }
        return det;
    }

    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> result;
        for (size_t i = 0; i < M; ++i) {
            result[i] = getColumn(i);
        }
        return result;
    }

    size_t rank() const {
        Matrix copy = gauss();
        size_t j = 0;
        size_t answer = 0;
        for (size_t i = 0; i < N; ++i) {
            while (j < M && copy.matrix[i][j] == Field(0)) {
                ++j;
            }
            if (j >= M) {
                break;
            }
            ++answer;
        }
        return answer;
    }

    Field trace() const {
        Field result = Field(0);
        for (size_t i = 0; i < std::min(N, M); ++i) {
            result += matrix[i][i];
        }
        return result;
    }

    Matrix inverted() const {
        static_assert(N == M);
        Matrix<N, 2 * N, Field> copy;
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> row = matrix[i];
            for (size_t j = 0; j < N; ++j) {
                row.push_back(Field(0));
                if (j == i) {
                    row.back() = Field(1);
                }
            }
            copy[i] = row;
        }
        copy = copy.gauss();
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> row = copy[i];
            for (size_t j = i + 1; j < 2 * N; ++j) {
                row[j] = row[j] / row[i];
            }
            row[i] = Field(1);
            copy[i] = row;
        }
        copy = copy.invertedGauss();
        Matrix result;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                result[i][j] = copy[i][j + N];
            }
        }
        return result;
    }

    Matrix &invert() {
        *this = inverted();
        return *this;
    }

    std::vector<Field> getRow(unsigned rowNumber) const {
        return matrix[rowNumber];
    }

    std::vector<Field> getColumn(unsigned columnNumber) const {
        std::vector<Field> result;
        for (size_t i = 0; i < N; ++i) {
            result.push_back(matrix[i][columnNumber]);
        }
        return result;
    }
};

template<size_t N, size_t M, typename Field = Rational>
std::ostream &operator<<(std::ostream &out, const Matrix<N, M, Field> &matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            out << int(matrix[i][j]) << ' ';
        }
        out << '\n';
    }
    return out;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field> operator*(const Field &multiplier, const Matrix<N, M, Field> &matrix) {
    Matrix<N, M, Field> result = matrix;
    result *= multiplier;
    return result;
}