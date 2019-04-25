#DEFINE MACROS:
CC       =  g++
FLAGS    = -std=c++11 -fopenmp
OPTIMIZE =  -unroll -O0 -g -Wno-deprecated -Wno-write-strings 

LIBS = -lgsl -lgslcblas -lgomp -lm

CFITSIO_LIB_DIR = ../Lib/cfitsio
HEALPIX_INC_DIR = ../Lib/Healpix_3.31/src/cxx/generic_gcc/include
HEALPIX_LIB_DIR = ../Lib/Healpix_3.31/src/cxx/generic_gcc/lib

LIBS_H =  -L$(HEALPIX_LIB_DIR) -L$(CFITSIO_LIB_DIR)  -lhealpix_cxx -lcxxsupport -lsharp -lc_utils -lfftpack -lcfitsio 
FLAGS_INC = -I$(HEALPIX_INC_DIR)	
#########################################################################
#ADD CODES:
#########################################################################
SOURCE1   = clgal.cpp 
OBJS1    = $(SOURCE1:.cpp=)
EXEC1     = $(SOURCE1:.cpp=.o)
$(EXEC1): $(OBJS1) 
$(OBJS1) : $(SOURCE1)
	$(CC) $(FLAGS) $(SOURCE1) -o $(EXEC1) $(FLAGS_INC) $(LIBS) $(LIBS_H) $(OPTIMIZE)  
.PHONY : clean
clean:
	rm -f  core.* *o *~
