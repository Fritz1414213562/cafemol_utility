CC=g++
CFLAG= -std=c++14 -pedantic -Wall -g
TEST= test.o CSVParser.o
MAKEDATA= main.o CSVParser.o 

test: $(TEST)
	$(CC) -o test $(TEST)

makedata: $(MAKEDATA)
	$(CC) -o makedata $(MAKEDATA)

clean:
	rm makedata main.o test.o CSVParser.o 
