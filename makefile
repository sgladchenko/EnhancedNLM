compiler = mpic++

# All of them
all:
	make run; make generate; make parser; make clean

# Test section:

test: Layer.o Vector.o Unit.o Z_test.o
	$(compiler) Unit.o Vector.o Layer.o Z_test.o -o test_containers; make clean

Z_test.o:
	$(compiler) -c ./source/Z_test.cpp -lm

# Main calculational unit in this project

run: Main.o Vector.o InitialSpectra.o Inhomogeneities.o Segmentation.o StandardScheme.o BaseBin.o Matrix.o StandardRK.o AdiabaticScheme.o Log.o Layer.o Unit.o AdiabaticFDM.o
	$(compiler) Main.o InitialSpectra.o Segmentation.o StandardScheme.o BaseBin.o Matrix.o StandardRK.o Inhomogeneities.o AdiabaticScheme.o Log.o Layer.o Vector.o Unit.o AdiabaticFDM.o -o run 

Main.o: ./source/Main.cpp ./source/Constants.h
	$(compiler) -c ./source/Main.cpp -lm

StandardScheme.o: ./source/StandardScheme.cpp ./source/StandardScheme.h ./source/BaseScheme.h ./source/Constants.h ./source/Log.h
	$(compiler) -c ./source/StandardScheme.cpp -lm

AdiabaticScheme.o: ./source/AdiabaticScheme.cpp ./source/AdiabaticScheme.h ./source/BaseScheme.h ./source/Constants.h ./source/Log.h
	$(compiler) -c ./source/AdiabaticScheme.cpp -lm	

InitialSpectra.o: ./source/InitialSpectra.cpp ./source/InitialSpectra.h ./source/PhysicalConstants.h
	$(compiler) -c ./source/InitialSpectra.cpp -lm

Segmentation.o: ./source/Segmentation.cpp ./source/Segmentation.h ./source/Constants.h
	$(compiler) -c ./source/Segmentation.cpp -lm

BaseBin.o: ./source/BaseBin.cpp ./source/BaseBin.h ./source/Constants.h
	$(compiler) -c ./source/BaseBin.cpp -lm

Matrix.o: ./source/Matrix.cpp ./source/Matrix.h ./source/Constants.h
	$(compiler) -c ./source/Matrix.cpp -lm

StandardRK.o: ./source/StandardRK.cpp ./source/StandardRK.h ./source/Constants.h
	$(compiler) -c ./source/StandardRK.cpp -lm

Inhomogeneities.o: ./source/Inhomogeneities.cpp ./source/Inhomogeneities.h ./source/Constants.h
	$(compiler) -c ./source/Inhomogeneities.cpp -lm

Log.o: ./source/Log.cpp ./source/Log.h ./source/Constants.h
	$(compiler) -c ./source/Log.cpp -lm

Unit.o: ./source/Unit.cpp ./source/Containers.h ./source/Constants.h
	$(compiler) -c ./source/Unit.cpp -lm

Vector.o: ./source/Vector.cpp ./source/Containers.h ./source/Constants.h
	$(compiler) -c ./source/Vector.cpp -lm	

Layer.o: ./source/Layer.cpp ./source/Containers.h ./source/Constants.h
	$(compiler) -c ./source/Layer.cpp -lm

AdiabaticFDM.o: ./source/AdiabaticFDM.cpp ./source/Containers.h ./source/Constants.h ./source/AdiabaticScheme.h
	$(compiler) -c ./source/AdiabaticFDM.cpp -lm

# Executable which generates for us initial data, saved in rec.bin

generate: InitialSpectra.o Generate.o Inhomogeneities.o ./source/Constants.h ./source/BaseBin.h
	$(compiler) InitialSpectra.o Generate.o Inhomogeneities.o  -o generate

Generate.o: ./source/Generate.cpp
	$(compiler) -c ./source/Generate.cpp -lm

# Additional executable which is used for processing data of full bin-files

parser: InitialSpectra.o Parser.o
	$(compiler) InitialSpectra.o Parser.o -o parser

Parser.o: ./source/Parser.cpp
	$(compiler) -c ./source/Parser.cpp -lm

# Clean all the object files

clean:
	rm -rf *.o
