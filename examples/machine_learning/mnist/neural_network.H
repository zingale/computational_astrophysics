#ifndef NEURAL_NETWORK_H
#define NEURAL_NETWORK_H

#include <algorithm>
#include <cmath>
#include <random>
#include <ranges>

#include "matrix.H"
#include "mnist.H"

// our sigmoid function

double sigmoid(double xi) {
    // the sigmoid activation function
    return 1.0 / (1.0 + std::exp(-xi));
}

double sigmoid_deriv(double z) {
    // this is the derivative of the sigmoid in terms of the
    // sigmoid, z = sigmoid(xi)
    return z * (1.0 - z);
}

class NeuralNetwork {

    int nin{};
    int nout{};
    int nhidden{};

    Matrix A;
    Matrix B;

public:

    NeuralNetwork(int input_size, int output_size, int hidden_layer_size)
        : nin(input_size), nout(output_size), nhidden(hidden_layer_size),
          A(nout, nhidden), B(nhidden, nin)
    {

        // initialize matrix elements to random values
        {
            std::random_device rd;
            std::mt19937 generator(rd());
            std::normal_distribution<double> randn(0.0, 1.0 / std::sqrt(nhidden));
            std::generate(A.flat().begin(), A.flat().end(),
                          [&]() -> double {return randn(generator);});
        }
        {
            std::random_device rd;
            std::mt19937 generator(rd());
            std::normal_distribution<double> randn(0.0, 1.0 / std::sqrt(nin));
            std::generate(B.flat().begin(), B.flat().end(),
                          [&]() -> double {return randn(generator);});
        }

    };

    void train(const std::vector<mnist::MNISTDigit>& training_data,
               const std::vector<mnist::MNISTDigit>& test_data,
               int n_epochs, double learning_rate, bool verbose=true) {

        std::random_device rd;
        std::mt19937 gen(rd());

        // vector of indices for randomly iterating over data
        std::vector<std::size_t> indices(training_data.size());
        std::iota(indices.begin(), indices.end(), 0);

        for (int i = 0; i < n_epochs; ++i) {
            if (verbose) {
                std::cout << "epoch " << i << " ... ";
            }
            std::shuffle(indices.begin(), indices.end(), gen);

            for (auto i : indices) {
                const auto& model = training_data[i];
                auto z_tilde = B * model.scaled;
                z_tilde.apply_inplace(sigmoid);

                auto z = A * z_tilde;
                z.apply_inplace(sigmoid);

                auto e = z - model.out;
                auto e_tilde = A.transpose() * (e % z.apply_new(sigmoid_deriv));

                // corrections
                auto dA = ((-2 * learning_rate * e) % z % (1.0 - z)) *
                    z_tilde.transpose();
                auto dB = ((-2 * learning_rate * e_tilde) % z_tilde.apply_new(sigmoid_deriv)) * model.scaled.transpose();

                A += dA;
                B += dB;
            }

            // finished an epoch -- see how we are doing so far
            int n_correct{0};
            for (auto model : test_data) {
                auto res = predict(model);
                if (model.validate(res)) {
                    n_correct += 1;
                }
            }
            if (verbose) {
                std::cout << "accuracy = "
                          << static_cast<double>(n_correct) /
                             static_cast<double>(test_data.size())
                          << std::endl;
            }
        }

    }

    Matrix
    predict(const mnist::MNISTDigit& model) {
        auto ztilde = B * model.scaled;
        ztilde.apply_inplace(sigmoid);
        auto z = A * ztilde;
        z.apply_inplace(sigmoid);

        return z;
    }
};

#endif
