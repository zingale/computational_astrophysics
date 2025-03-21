#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <cassert>

// a 2-d matrix with contiguous storage
// here the data is stored in row-major order in a 1-d memory space
// managed as a vector.  We overload () to allow us to index this as
// a(irow, icol)

class Matrix {

    std::size_t _rows;
    std::size_t _cols;
    std::vector<double> _data;

public:

    Matrix (std::size_t rows, std::size_t cols, double val=0.0)
        : _rows{rows},
          _cols{cols},
          _data(rows * cols, val)
    {
        assert (rows > 0 && cols > 0);
    }

    // Constructor for list initialization
    // e.g., Matrix A{{1, 2}, {3, 4}};

    Matrix(std::initializer_list<std::initializer_list<double>> values)
        : _rows(values.size()),
          _cols(values.begin()->size()),
          _data()
    {
        assert(_rows > 0);
        for (const auto& row : values) {
            // Ensure all rows have the same number of columns
            assert(row.size() == _cols);
            _data.insert(_data.end(), row.begin(), row.end());
        }
    }

    inline std::vector<double>& flat() {return _data;}

    // note the "const" after the argument list here -- this means
    // that this can be called on a const Matrix

    [[nodiscard]] inline std::size_t ncols() const { return _cols;}
    [[nodiscard]] inline std::size_t nrows() const { return _rows;}

    inline double& operator()(std::size_t row, std::size_t col) {
        // size_t is unsigned, so will always be >= 0
        assert (row < _rows);
        assert (col < _cols);
        return _data[row*_cols + col];
    }

    inline const double& operator()(std::size_t row, std::size_t col) const {
        // size_t is unsigned, so will always be >= 0
        assert (row < _rows);
        assert (col < _cols);
        return _data[row*_cols + col];
    }

    [[nodiscard]] Matrix transpose() const {
        Matrix A_T(_cols, _rows);
        for (std::size_t irow = 0; irow < _rows; ++irow) {
            for (std::size_t jcol = 0; jcol < _cols; ++jcol) {
                A_T(jcol, irow) = (*this)(irow, jcol);
            }
        }
        return A_T;
    }


    // matrix-vector multiplication
    // b = A x -- we make the result a column-matrix

    Matrix
    operator* (const std::vector<double>& x) {
        assert (_cols == x.size());
        Matrix b(_rows, 1);

        for (std::size_t irow = 0; irow < _rows; ++irow) {
            for (std::size_t k = 0; k < _cols; ++k) {
                b(irow, 0) += (*this)(irow, k) * x[k];
            }
        }

        return b;
    }

    // matrix-matrix multiplication,
    // C = A B

    Matrix
    operator* (const Matrix& B) {
        assert (_cols == B.nrows());
        Matrix C(_rows, B.ncols());

        for (std::size_t irow = 0; irow < C.nrows(); ++irow) {
            for (std::size_t jcol = 0; jcol < C.ncols(); ++jcol) {
                for (std::size_t k = 0; k < _cols; ++k) {
                    C(irow, jcol) += (*this)(irow, k) * B(k, jcol);
                }
            }
        }

        return C;
    }

    // multiply by a scalar, C = A * a

    Matrix
    operator* (double a) {
        Matrix C(_rows, _cols);
        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = a * _data[n];
        }
        return C;
    }

    friend Matrix
    operator* (double a, const Matrix& A);

    // element-by-element multiplication

    Matrix
    operator% (const Matrix& B) {
        assert (_rows == B.nrows() && _cols == B.ncols());
        Matrix C(_rows, _cols);

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = _data[n] * B._data[n];
        }

        return C;
    }

    // addition

    Matrix
    operator+ (const Matrix& B) {
        assert (_rows == B.nrows() && _cols == B.ncols());
        Matrix C(_rows, _cols);

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = _data[n] + B._data[n];
        }

        return C;
    }

    Matrix
    operator+ (double a) {
        Matrix C(_rows, _cols);

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = _data[n] + a;
        }

        return C;
    }

    friend Matrix
    operator+ (double a, const Matrix& A);

    // subtraction

    Matrix
    operator- (const Matrix& B) {
        assert (_rows == B.nrows() && _cols == B.ncols());
        Matrix C(_rows, _cols);

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = _data[n] - B._data[n];
        }

        return C;
    }

    Matrix
    operator- (double a) {
        Matrix C(_rows, _cols);

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            C._data[n] = _data[n] - a;
        }

        return C;
    }

    friend Matrix
    operator- (double a, const Matrix& A);

    // compound operators

    Matrix&
    operator+= (const Matrix& B) {
        assert(_rows == B.nrows() && _cols == B.ncols());

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            _data[n] += B._data[n];
        }
        return *this;
    }

    Matrix&
    operator-= (const Matrix& B) {
        assert(_rows == B.nrows() && _cols == B.ncols());

        for (std::size_t n = 0; n < _rows * _cols; ++n) {
            _data[n] -= B._data[n];
        }
        return *this;
    }

    // apply a function on element-by-element basis

    void apply(const std::function<double(double)>& f) {
        for (auto& e : _data) {
            e = f(e);
        }
    }

};

// right-multiply by a scalar, C = a * A


inline Matrix
operator* (double a, const Matrix& A) {
    Matrix C(A._rows, A._cols);
    for (std::size_t n = 0; n < A._rows * A._cols; ++n) {
        C._data[n] = a * A._data[n];
    }
    return C;
}

inline Matrix
operator+ (double a, const Matrix& A) {
    Matrix C(A.nrows(), A.ncols());

    for (std::size_t n = 0; n < A.nrows() * A.ncols(); ++n) {
        C._data[n] = a + A._data[n];
    }

    return C;
}

inline Matrix
operator- (double a, const Matrix& A) {
    Matrix C(A.nrows(), A.ncols());

    for (std::size_t n = 0; n < A.nrows() * A.ncols(); ++n) {
        C._data[n] = a - A._data[n];
    }

    return C;
}

// the << operator is not part of the of the class, so it is not a
// member

inline
std::ostream& operator<< (std::ostream& os, const Matrix& a) {

    for (std::size_t row = 0; row < a.nrows(); ++row) {
        for (std::size_t col = 0; col < a.ncols(); ++col) {
            os << a(row, col) << " ";
        }
        os << std::endl;
    }
    return os;
}

#endif
