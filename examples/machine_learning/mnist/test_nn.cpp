#include <iostream>

#include "mnist.H"
#include "neural_network.H"


int main() {

    // train the network using the training set

    std::cout << "reading in the data...";

    auto training_set = mnist::read_training_set();
    auto test_set = mnist::read_test_set();

    std::cout << "done\n" << std::endl;

    constexpr int hidden_layer_size{50};
    constexpr int n_epochs{10};
    constexpr double learning_rate{0.1};

    auto n = NeuralNetwork(mnist::DIGIT_ROWS * mnist::DIGIT_COLS,
                           mnist::DIGIT_CATEGORIES, hidden_layer_size);
    n.train(training_set, test_set, n_epochs, learning_rate);

    // now assess how well we did using the test set

    int n_correct{0};
    for (auto model : test_set) {
        auto res = n.predict(model);
        if (model.validate(res)) {
            n_correct += 1;
        }
    }

    std::cout << "finished training!" << std::endl;
    std::cout << "final test set accuracy = "
              << static_cast<double>(n_correct) /
                 static_cast<double>(test_set.size())
              << std::endl << std::endl;

#if 0
    // output the first 10 we get wrong

    int n_wrong{0};
    for (auto model : test_set) {
        auto res = n.predict(model);
        if (!model.validate(res)) {
            std::cout << "correct answer = " << model.num << "\n"
                      << "prediction = " << model.interpret(res) << std::endl;
            model.display();
            n_wrong += 1;
        }
        if (n_wrong >= 10) {
            break;
        }
    }
#endif

    std::cout << "Exploring hidden layer size:" << std::endl;

    // now explore the size of the hidden layer
    for (auto hl_size : {25, 50, 100, 150}) {
        auto n = NeuralNetwork(mnist::DIGIT_ROWS * mnist::DIGIT_COLS,
                               mnist::DIGIT_CATEGORIES, hl_size);
        constexpr bool verbose{false};
        n.train(training_set, test_set, n_epochs,
                learning_rate, verbose);

        // now assess how well we did using the test set

        int n_correct{0};
        for (auto model : test_set) {
            auto res = n.predict(model);
            if (model.validate(res)) {
                n_correct += 1;
            }
        }

        std::cout << "hidden layer size = " << hl_size
                  << " accuracy = "
                  << static_cast<double>(n_correct) /
                     static_cast<double>(test_set.size())
                  << std::endl;
    }

}
