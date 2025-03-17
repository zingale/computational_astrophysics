#ifndef NEURAL_NETWORK_H
#define NEURAL_NETWORK_H

#include <algorithm>
#include <cmath>
#include <random>
#include <ranges>

#include "matrix.H"
#include "mnist.H"

// our sigmoid function

double g(double xi) {
    return 1.0 / (1.0 + std::exp(-xi));
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
            std::generate(A._data.begin(), A._data.end(),
                          [&]() -> double {return randn(generator);});
        }
        {
            std::random_device rd;
            std::mt19937 generator(rd());
            std::normal_distribution<double> randn(0.0, 1.0 / std::sqrt(nin));
            std::generate(B._data.begin(), B._data.end(),
                          [&]() -> double {return randn(generator);});
        }

    };

    void train(const std::vector<MNISTDigit>& training_data,
               int n_epochs, double learning_rate=0.1) {

        std::random_device rd;
        std::mt19937 gen(rd());

        for (int i = 0; i < n_epochs; ++i) {
            std::cout << "epoch " << i << std::endl;

            for (const auto& model : training_data) { //std::ranges::shuffle(training_data, gen)) {
                auto z_tilde = B * model.scaled;
                z_tilde.apply(g);

                auto z = A * z_tilde;
                z.apply(g);

                auto e = z - model.out;
                auto e_tilde = A.transpose() * (e % z % (1.0 - z));

                // corrections
                auto dA = ((-2 * learning_rate * e) % z % (1.0 - z)) *
                    z_tilde.transpose();
                auto dB = ((-2 * learning_rate * e_tilde) % z_tilde % (1.0 - z_tilde)) * model.scaled.transpose();

                A += dA;
                B += dB;
            }
        }
    }

    Matrix
    predict(const MNISTDigit& model) {
        auto ztilde = B * model.scaled;
        ztilde.apply(g);
        auto z = A * ztilde;
        z.apply(g);

        return z;
    }
};

#endif
