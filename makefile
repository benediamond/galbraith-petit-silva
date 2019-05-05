CXX = g++
CXXFLAGS = -isystem/usr/local/include -std="c++11" -g2 # -L/usr/local/lib -lLiDIA -lgmp -lm -lcryptopp

main: main.o gps.o walker.o order.o basis.o step.o extension.o
	$(CXX) $(CXXFLAGS) -o main.o gps.o walker.o order.o basis.o step.o extension.o
main.o: main.cc gps.h walker.h step.h order.h basis.h extension.h
	$(CXX) $(CXXFLAGS) -c main.cc
gps.o: gps.h walker.h step.h order.h basis.h extension.h

walker.o: walker.h step.h order.h basis.h extension.h

order.o: order.h step.h basis.h extension.h

step.o: step.h extension.h

basis.o: basis.h extension.h

extension.o: extension.h
