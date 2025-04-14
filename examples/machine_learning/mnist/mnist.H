#ifndef MNIST_H
#define MNIST_H

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <ranges>
#include <vector>

#include "matrix.H"

namespace mnist {

    constexpr std::size_t DIGIT_ROWS{28};
    constexpr std::size_t DIGIT_COLS{28};


    // a single digit from MNIST

    class MNISTDigit {

    public:

        Matrix scaled;
        Matrix out;
        int num{-1};

        MNISTDigit(const std::string& raw_string)
            : scaled(DIGIT_ROWS * DIGIT_COLS, 1), out(10, 1, 0.01)
        {
            // loop over the raw string and split it up.  The first number
            // is the answer, the rest are the digit

            std::stringstream ss(raw_string);
            std::string value;
            int n = 0;
            while (std::getline(ss, value, ',')) {
                if (n == 0) {
                    num = std::stoi(value);
                } else {
                    scaled(n-1, 0) = std::stod(value) / 255.0 * 0.99 + 0.01;
                }
                n++;
            }

            // set the categorical data
            out(num, 0) = 0.99;

        }

        void display() {

            for (auto [n, e] : scaled.flat() | std::views::enumerate) {

                if (n % DIGIT_COLS == 0) {
                    std::cout << std::endl;
                }

                if (e > 0.8) {
                    std::cout << "█";
                } else if (e > 0.6) {
                    std::cout << "▓";
                } else if (e > 0.4) {
                    std::cout << "▒";
                } else if (e > 0.2) {
                    std::cout << "░";
                } else {
                    std::cout << " ";
                }

            }
            std::cout << std::endl;
        }

        // given a vector of 10 doubles, take the index with the largest
        // value as the guess for the digit and compare to the correct
        // answer

        int interpret(const Matrix& categorical_guess) {
            assert (categorical_guess.nrows() == 10);

            double max_val{std::numeric_limits<double>::lowest()};
            int idx{-1};
            for (std::size_t n = 0; n < categorical_guess.nrows(); ++n) {
                if (categorical_guess(n, 0) > max_val) {
                    max_val = categorical_guess(n, 0);
                    idx = static_cast<int>(n);
                }
            }
            return idx;

        }

        bool validate(const Matrix& categorical_guess) {
            return interpret(categorical_guess) == num;
        }

    };

    inline
    std::vector<MNISTDigit> read_training_set() {

        std::vector<MNISTDigit> digits;
        std::ifstream mf("mnist_train.csv");
        if (mf.fail()) {
            std::cout << "Error: training file does not exist" << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(mf, line)) {
            digits.emplace_back(line);
        }
        return digits;
    }

    inline
    std::vector<MNISTDigit> read_test_set() {

        std::vector<MNISTDigit> digits;

        std::ifstream mf("mnist_test.csv");
        if (mf.fail()) {
            std::cout << "Error: test file does not exist" << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(mf, line)) {
            digits.emplace_back(line);
        }
        return digits;
    }
}
#endif
