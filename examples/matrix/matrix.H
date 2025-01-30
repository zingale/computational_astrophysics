#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <cassert>

// a simple matrix class that supports matrix-matrix multiplication
// addition of matrices, and multiplication by a scalar.
//
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
        // we do the asserts here after the initialization of _data
        // in the initialization list, but if the size is zero, we'll
        // abort here anyway.

        assert (rows > 0);
        assert (cols > 0);

    }

    // constructor that allows us to initialize Matrix
    // via an initializer list of rows

    explicit Matrix (const std::vector<std::vector<double>>& v)
        : _rows{v.size()},
          _cols{v[0].size()},
          _data(_rows * _cols, 0.0)
    {
        int idx = 0;
        for (std::size_t i = 0; i < _rows; ++i) {
            assert (v[i].size() == _cols);
            for (std::size_t j = 0; j < _cols; ++j) {
                _data[idx] = v[i][j];
                ++idx;
            }
        }
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

    inline Matrix operator* (double c) {
        Matrix A(*this);
        for (std::size_t irow = 0; irow < nrows(); ++irow) {
            for (std::size_t icol = 0; icol < ncols(); ++icol) {
                A(irow, icol) *= c;
            }
        }
        return A;
    }

    inline Matrix operator/ (double c) {
        Matrix A(*this);
        for (std::size_t irow = 0; irow < nrows(); ++irow) {
            for (std::size_t icol = 0; icol < ncols(); ++icol) {
                A(irow, icol) /= c;
            }
        }
        return A;
    }

    inline Matrix operator+ (const Matrix& B) {
        // A = this + B
        assert(B.nrows() == nrows() && B.ncols() == ncols());

        Matrix A(*this);
        for (std::size_t irow = 0; irow < nrows(); ++irow) {
            for (std::size_t icol = 0; icol < ncols(); ++icol) {
                A(irow, icol) = (*this)(irow, icol) + B(irow, icol);
            }
        }
        return A;
    }

    inline Matrix operator- (const Matrix& B) {
        // A = this - B
        assert(B.nrows() == nrows() && B.ncols() == ncols());

        Matrix A(*this);
        for (std::size_t irow = 0; irow < nrows(); ++irow) {
            for (std::size_t icol = 0; icol < ncols(); ++icol) {
                A(irow, icol) = (*this)(irow, icol) - B(irow, icol);
            }
        }
        return A;
    }

};

// the << operator is not part of the of the class, so it is not a member

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
