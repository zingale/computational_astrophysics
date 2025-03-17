#include <iostream>

#include "mnist.H"
#include "neural_network.H"


int main() {

    auto training_set = read_training_set();

    auto n = NeuralNetwork(784, 10, 50);
    n.train(training_set, 5);

    auto test_set = read_test_set();

    int n_correct{0};
    for (auto model : test_set) {
        auto res = n.predict(model);
        if (model.validate(res)) {
            n_correct += 1;
        }
    }

    std::cout << "accuracy = " << static_cast<double>(n_correct) / test_set.size() << std::endl;


}
