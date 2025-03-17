#ifndef MNIST_H
#define MNIST_H

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ranges>
#include <vector>

constexpr int DIGIT_ROWS{28};
constexpr int DIGIT_COLS{28};


// a single digit from MNIST

class MNISTDigit {

public:

    std::vector<double> scaled;
    std::vector<double> out;
    int num{-1};

    MNISTDigit(const std::string& raw_string)
        : scaled(DIGIT_ROWS * DIGIT_COLS, 0.0), out(10, 0.01)
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
                scaled[n-1] = std::stod(value) / 255.0 * 0.99 + 0.01;
            }
            n++;
        }

        // set the categorical data
        out[num] = 0.99;

    }

    void display() {

        for (auto [n, e] : scaled | std::views::enumerate) {

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

    bool validate(const std::vector<double>& categorical_guess) {
        assert (categorical_guess.size() == 10);

        auto max_iter = std::max_element(categorical_guess.cbegin(),
                                         categorical_guess.cend());
        std::size_t index = std::distance(categorical_guess.cbegin(), max_iter);

        return static_cast<int>(index) == num;
    }
};

std::vector<MNISTDigit> read_training_set() {

    std::vector<MNISTDigit> digits;

    std::ifstream mf("mnist_train.csv");
    std::string line;
    while (std::getline(mf, line)) {
        digits.emplace_back(line);
    }
    return digits;
}

std::vector<MNISTDigit> read_test_set() {

    std::vector<MNISTDigit> digits;

    std::ifstream mf("mnist_test.csv");
    std::string line;
    while (std::getline(mf, line)) {
        digits.emplace_back(line);
    }
    return digits;
}

#endif
