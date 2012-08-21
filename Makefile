all : psifit draw

LIBS = `root-config --libs` -lMinuit -L$(HOME)/local/lib  -lboost_program_options
CXXFLAGS = `root-config --cflags` -I$(HOME)/local/include -I$(HOME)/work/ -std=c++0x 

BINDIR=$(HOME)/work/bin
psifit : FitOniumRCompactLib.o psifit.o
		g++ -o $(BINDIR)/psifit  FitOniumRCompactLib.o psifit.o  $(LIBS)

draw : FitOniumRCompactLib.o draw.o 
		g++ -o $(BINDIR)/draw FitOniumRCompactLib.o  draw.o   $(LIBS)

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

psifit.o : FitScan.cc
		g++ -o psifit.o  $(CXXFLAGS) -c FitScan.cc 

#interference.o :
#		g++ -o interference.o  $(CXXFLAGS) -c ../interference.cpp
	
	
clean :
	rm $(BINDIR)/psifit *.o $(BINDIR)/draw
