all:
	g++ -Ofast -I ../Eigen -march=native -std=c++20 -o main.exe main.cpp
	main.exe
parallel:
	g++ -fopenmp -O3 -I ../Eigen -march=native -std=c++20 -o main.exe main.cpp 
	main.exe
debug:
	g++ -O1 -g3 -I ../Eigen -march=native -std=c++20 -o main.exe main.cpp -ftree-vectorize -funsafe-math-optimizations -fopenmp
	main.exe