#include <iostream>

#include "mnist.H"
#include "neural_network.H"


int main() {

    // train the network using the training set

    auto training_set = read_training_set();

    auto n = NeuralNetwork(784, 10, 50);
    n.train(training_set, 5);

    // now assess how well we did using the test set

    auto test_set = read_test_set();

    int n_correct{0};
    for (auto model : test_set) {
        auto res = n.predict(model);
        if (model.validate(res)) {
            n_correct += 1;
        }
    }

    std::cout << "accuracy = " << static_cast<double>(n_correct) / test_set.size() << std::endl;

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

}
