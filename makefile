compiler = mpic++

# All of them
all:
	make run; make generate; make parser; make clean

# Main calculational unit in this project

run: main.o matrix.o InitialSpectra.o node.o inhomogeneities.o segmentation.o io.o 
	$(compiler) main.o node.o InitialSpectra.o matrix.o inhomogeneities.o segmentation.o io.o -o run

main.o: ./source/main.cpp node.o
	$(compiler) -c ./source/main.cpp -lm

node.o: io.o matrix.o inhomogeneities.o segmentation.o ./source/node.cpp 
	$(compiler) -c ./source/node.cpp -lm

matrix.o: ./source/matrix.cpp
	$(compiler) -c ./source/matrix.cpp -lm

InitialSpectra.o: ./source/InitialSpectra.cpp
	$(compiler) -c ./source/InitialSpectra.cpp -lm

inhomogeneities.o: ./source/inhomogeneities.cpp
	$(compiler) -c ./source/inhomogeneities.cpp -lm

io.o: ./source/io.cpp segmentation.o
	$(compiler) -c ./source/io.cpp -lm

segmentation.o: ./source/segmentation.cpp
	$(compiler) -c ./source/segmentation.cpp -lm

# Executable which generates for us initial data, saved in rec.bin

generate: InitialSpectra.o generate.o inhomogeneities.o
	$(compiler) InitialSpectra.o generate.o inhomogeneities.o  -o generate

generate.o: ./source/generate.cpp
	$(compiler) -c ./source/generate.cpp -lm

# Additional executable which is used for processing data of full bin-files

parser: InitialSpectra.o parser.o
	$(compiler) InitialSpectra.o parser.o -o parser

parser.o: ./source/parser.cpp
	$(compiler) -c ./source/parser.cpp -lm

clean:
	rm -rf *.o