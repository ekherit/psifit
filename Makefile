all : psifit draw

LIBS = `root-config --libs` -lMinuit -L$(HOME)/local/lib  -lboost_program_options
CXXFLAGS = `root-config --cflags` -I$(HOME)/local/include -I$(HOME)/work/ -std=c++11

CC=g++

BINDIR=$(HOME)/work/bin
psifit : FitOniumRCompactLib.o psifit.o
		$(CC) -o $(BINDIR)/psifit  FitOniumRCompactLib.o psifit.o  $(LIBS)

psifit.o : psifit.cpp
		$(CC) -o psifit.o  $(CXXFLAGS) -c psifit.cpp

draw : FitOniumRCompactLib.o draw.o 
		$(CC) -o $(BINDIR)/draw FitOniumRCompactLib.o  draw.o   $(LIBS)

FitOniumRCompactLib.o : FitOniumRCompactLib.cc
		$(CC) -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc


#interference.o :
#		g++ -o interference.o  $(CXXFLAGS) -c ../interference.cpp
	
	
clean :
	rm $(BINDIR)/psifit *.o $(BINDIR)/draw
