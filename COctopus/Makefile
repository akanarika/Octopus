all: main

main: main.o phyxel.o objloader.o
	g++ main.o phyxel.o objloader.o -o main -framework OpenGL -framework GLUT -lGLEW

main.o: main.cpp
	g++ -c main.cpp

phyxel.o: phyxel.cpp
	g++ -c phyxel.cpp

objloader.o: objloader.cpp
	g++ -c objloader.cpp
	
