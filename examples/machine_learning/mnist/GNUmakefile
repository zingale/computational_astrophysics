ALL: test_nn

HEADERS := $(wildcard *.H)

%.o: %.cpp $(HEADERS)
	g++ -Wall -Wextra -Wpedantic -std=c++23 -DNDEBUG -c -O3 $<

test_nn: test_nn.o $(HEADERS)
	g++ -o $@ $<
