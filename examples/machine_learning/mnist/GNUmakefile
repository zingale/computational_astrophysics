ALL: test_nn

HEADERS := $(wildcard *.H)

%.o: %.cpp
	g++ -std=c++23 -DNDEBUG -c -O3 $<

test_nn: test_nn.o $(HEADERS)
	g++ -o $@ $<
