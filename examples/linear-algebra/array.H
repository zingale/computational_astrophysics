#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include <iostream>
#include <cassert>

// a contiguous 2-d array
// here the data is stored in row-major order in a 1-d memory space
// managed as a vector.  We overload () to allow us to index this as
// a(irow, icol)

struct Array {

    int _rows;
    int _cols;
    std::vector<double> _data;

    Array (int rows, int cols, double val=0.0)
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

    // constructor that allows us to initialize array
    // via an initializer list of rows

    Array (std::vector<std::vector<double>> v)
        : _rows{static_cast<int>(v.size())},
          _cols{static_cast<int>(v[0].size())},
          _data(_rows * _cols, 0.0)
    {
        int idx = 0;
        for (int i = 0; i < _rows; ++i) {
            assert (static_cast<int>(v[i].size()) == _cols);
            for (int j = 0; j < _cols; ++j) {
                _data[idx] = v[i][j];
                ++idx;
            }
        }
    }

    // note the "const" after the argument list here -- this means
    // that this can be called on a const Array

    inline int ncols() const { return _cols;}
    inline int nrows() const { return _rows;}

    inline double& operator()(int row, int col) {
        assert (row >= 0 && row < _rows);
        assert (col >= 0 && col < _cols);
        return _data[row*_cols + col];
    }

    inline const double& operator()(int row, int col) const {
        assert (row >= 0 && row < _rows);
        assert (col >= 0 && col < _cols);
        return _data[row*_cols + col];
    }

};

// the << operator is not part of the of the class, so it is not a member

std::ostream& operator<< (std::ostream& os, const Array& a) {

    for (int row = 0; row < a.nrows(); ++row) {
        for (int col = 0; col < a.ncols(); ++col) {
            os << a(row, col) << " ";
        }
        os << std::endl;
    }

    return os;
}

#endif
