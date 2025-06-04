CC = g++
LIBPATH =

#CFLAGS = -O3 -Wall -g -D HAVE_PTHREADS
#CFLAGS = -O3 -std=c++11 -g -Wall
CFLAGS = -O3 -std=c++11 -Wall

all: scist


Utils.o : Utils.cpp  Utils.h
	$(CC)  $(CFLAGS) -c Utils.cpp

Utils2.o : Utils2.cpp  Utils2.h Utils.h
	$(CC)  $(CFLAGS) -c Utils2.cpp

Utils3.o : Utils3.cpp  Utils3.h Utils.h
	$(CC)  $(CFLAGS) -c Utils3.cpp

Utils4.o : Utils4.cpp  Utils4.h Utils.h
	$(CC)  $(CFLAGS) -c Utils4.cpp

UtilsNumerical.o : UtilsNumerical.cpp  UtilsNumerical.h
	$(CC)  $(CFLAGS) -c UtilsNumerical.cpp

#Utils-basic.o : Utils-basic.cpp  Utils-basic.h
#	$(CC)  $(CFLAGS) -c Utils-basic.cpp

#GraphUtils.o : GraphUtils.cpp  GraphUtils.h
#	$(CC)  $(CFLAGS) -c GraphUtils.cpp

RerootTreeUtils.o : RerootTreeUtils.cpp  RerootTreeUtils.h
	$(CC)  $(CFLAGS) -c RerootTreeUtils.cpp

TreeBuilder.o : TreeBuilder.cpp  TreeBuilder.h
	$(CC)  $(CFLAGS) -c TreeBuilder.cpp

UnWeightedGraph.o : UnWeightedGraph.cpp  UnWeightedGraph.h
	$(CC)  $(CFLAGS) -c UnWeightedGraph.cpp

MarginalTree.o : MarginalTree.cpp MarginalTree.h
	$(CC)  $(CFLAGS) -c MarginalTree.cpp

RBT.o : RBT.cpp RBT.h
	$(CC)  $(CFLAGS) -c RBT.cpp

PhylogenyTreeBasic.o : PhylogenyTreeBasic.cpp PhylogenyTreeBasic.h
	$(CC)  $(CFLAGS) -c PhylogenyTreeBasic.cpp

PhylogenyTree.o : PhylogenyTree.cpp PhylogenyTree.h
	$(CC)  $(CFLAGS) -c PhylogenyTree.cpp

BioSequenceMatrix.o : BioSequenceMatrix.cpp BioSequenceMatrix.h
	$(CC)  $(CFLAGS) -c BioSequenceMatrix.cpp

BinaryMatrix.o : BinaryMatrix.cpp BinaryMatrix.h
	$(CC)  $(CFLAGS) -c BinaryMatrix.cpp

GenotypeMatrix.o : GenotypeMatrix.cpp GenotypeMatrix.h
	$(CC)  $(CFLAGS) -c GenotypeMatrix.cpp

ScistGenotype.o : ScistGenotype.cpp ScistGenotype.hpp
	$(CC)  $(CFLAGS) -c ScistGenotype.cpp

ScistPerfPhyUtils.o : ScistPerfPhyUtils.cpp ScistPerfPhyUtils.hpp
	$(CC)  $(CFLAGS) -c ScistPerfPhyUtils.cpp

ScistPerfPhyImp.o : ScistPerfPhyImp.cpp ScistPerfPhyImp.hpp
	$(CC)  $(CFLAGS) -c ScistPerfPhyImp.cpp

ScistDoublet.o : ScistDoublet.cpp ScistDoublet.hpp
	$(CC)  $(CFLAGS) -c ScistDoublet.cpp

ScistErrRateInf.o : ScistErrRateInf.cpp ScistErrRateInf.hpp
	$(CC)  $(CFLAGS) -c ScistErrRateInf.cpp

ScistLocalSearch.o : ScistLocalSearch.cpp ScistLocalSearch.hpp
	$(CC)  $(CFLAGS) -c ScistLocalSearch.cpp

main.o : main.cpp
	$(CC)  $(CFLAGS) -c main.cpp

# Now final link

scist: Utils.o Utils2.o  Utils3.o Utils4.o UtilsNumerical.o UnWeightedGraph.o TreeBuilder.o RerootTreeUtils.o \
    BioSequenceMatrix.o BinaryMatrix.o GenotypeMatrix.o PhylogenyTreeBasic.o PhylogenyTree.o MarginalTree.o RBT.o \
    ScistGenotype.o ScistPerfPhyUtils.o ScistDoublet.o ScistErrRateInf.o ScistPerfPhyImp.o ScistLocalSearch.o main.o
	$(CC) -o scistree Utils.o Utils2.o  Utils3.o Utils4.o UtilsNumerical.o UnWeightedGraph.o TreeBuilder.o RerootTreeUtils.o BioSequenceMatrix.o  \
       BinaryMatrix.o GenotypeMatrix.o PhylogenyTreeBasic.o PhylogenyTree.o MarginalTree.o RBT.o \
       ScistGenotype.o ScistPerfPhyUtils.o ScistPerfPhyImp.o ScistDoublet.o ScistErrRateInf.o ScistLocalSearch.o main.o \
        -lm -lpthread
