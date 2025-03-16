#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <cassert>

// a 2-d matrix with contiguous storage
// here the data is stored in row-major order in a 1-d memory space
// managed as a vector.  We overload () to allow us to index this as
// a(irow, icol)

struct Matrix {

    std::size_t _rows;
    std::size_t _cols;
    std::vector<double> _data;

    Matrix (std::size_t rows, std::size_t cols, double val=0.0)
        : _rows{rows},
          _cols{cols},
          _data(rows * cols, val)
    {
        assert (rows > 0 && cols > 0);
    }

    // note the "const" after the argument list here -- this means
    // that this can be called on a const Matrix

    inline std::size_t ncols() const { return _cols;}
    inline std::size_t nrows() const { return _rows;}

    inline double& operator()(int row, int col) {
        assert (row >= 0 && row < static_cast<int>(_rows));
        assert (col >= 0 && col < static_cast<int>(_cols));
        return _data[row*_cols + col];
    }

    inline const double& operator()(int row, int col) const {
        assert (row >= 0 && row < static_cast<int>(_rows));
        assert (col >= 0 && col < static_cast<int>(_cols));
        return _data[row*_cols + col];
    }

    Matrix transpose() {
        Matrix A_T(_cols, _rows);
        for (std::size_t irow = 0; irow < _rows; ++irow) {
            for (std::size_t jcol = 0; jcol < _cols; ++jcol) {
                A_T(jcol, irow) = (*this)(irow, jcol);
            }
        }
        return A_T;
    }

    std::vector<double>
    operator* (const std::vector<double>& x) {
        // do b = A x
        assert (_cols == x.size());

        std::vector<double> b(_rows, 0.0);

        for (std::size_t irow = 0; irow < _rows; ++irow) {
            for (std::size_t k = 0; k < _cols; ++k) {
                b[irow] += (*this)(irow, k) * x[k];
            }
        }

        return b;
    }

    Matrix
    operator* (const Matrix& B) {
        // do C = A B
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

};

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
