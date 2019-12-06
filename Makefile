CC = g++
CFLAGS = -Wall

all: main ./test/tests

main: ./build/cluster.o ./build/utilities.o ./build/classification.o ./build/Grid.o ./build/LSH_Structure.o
	$(CC) -o ./bin/cluster ./build/cluster.o ./build/utilities.o ./build/classification.o ./build/Grid.o ./build/LSH_Structure.o

./test/tests: ./tests.o ./build/utilities.o ./build/classification.o ./build/Grid.o ./build/LSH_Structure.o
	g++ ./tests.o ./build/utilities.o ./build/classification.o ./build/Grid.o ./build/LSH_Structure.o -o ./test/tests -Wall -g -pthread -lgtest_main -lgtest -lpthread

./tests.o: ./test/tests.cpp
	g++ -Wall -g -pthread -lgtest_main -lgtest -lpthread -std=c++11 -c ./test/tests.cpp

# Main code
./build/cluster.o: ./src/cluster.cpp
	$(CC) -o ./build/cluster.o -c ./src/cluster.cpp

# Libraries code
./build/utilities.o: ./src/utilities.cpp
	$(CC) ${CFLAGS} -c -o ./build/utilities.o ./src/utilities.cpp

./build/classification.o: ./src/classification.cpp
	$(CC) $(CFLAGS) -c -o ./build/classification.o ./src/classification.cpp

./build/LSH_Structure.o: ./src/LSH_Structure.cpp
	$(CC) $(CFLAGS) -c -o ./build/LSH_Structure.o ./src/LSH_Structure.cpp

./build/Grid.o: ./src/Grid.cpp
	$(CC) $(CFLAGS) -c -o ./build/Grid.o ./src/Grid.cpp

# Clean
clean:
	rm ./bin/* ./test/tests ./tests.o