UNAME=$(shell uname -s)
LIBNAME:=libnetdecode.a


CXXFLAGS=-g3 -std=c++11 -pthread -fopenmp -w -m64 -fPIC -march=native -O4 -DUSE_BLAS_LIB -DAXPBY -DINT_64BITS -DDEBUG
ifeq ($(UNAME),Linux)
	CXX=g++
	LINALG=-lopenblas -llapack
	#MKLROOT=/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl
	#LINALG=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lm -ldl
	#LINALG=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -ldl
endif	
ifeq ($(UNAME),Darwin)
	CXX=/usr/local/opt/llvm/bin/clang++
	LINALG=-framework accelerate
	CXXFLAGS+=-DACCELERATE
endif

LIB_FLAGS=-lstdc++ ${LINALG} -lpthread
INCLUDE=-I./include -I./include/armadillo/include

SRC=src/NetDECODE.cc
OBJ=$(SRC:.cc=.o)

MATLAB_SRC=matlab/asssessCoactivity.cc matlab/constructKstarNN.cc matlab/predictActivityScores.cc
MATLAB_MEX=$(MATLAB_SRC:.cc=.mexa64)

MATLAB_FLAGS=-DUSE_BLAS_LIB -DAXPBY -DINT_64BITS -DNDEBUG -largeArrayDims

MATLAB = $(shell matlab -e | sed -n 's/MATLAB=//p')
MEX = $(MATLAB)/bin/mex

PROGRAM=NetDECODE
	
all: $(LIBNAME) message
	
src/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) ${INCLUDE} -c -o $@  $<


$(LIBNAME): $(OBJ)
	ar -rcs $@ $(OBJ)

NetDECODE: $(LIBNAME) src/main.o
	$(CXX) $(CXXFLAGS)  -o $@ src/main.o -L. ${LIBNAME} $(LIB_FLAGS) 

	
matlab/%.mexa64: matlab/%.cc
	$(MEX) $< -outdir matlab/ $(LIBNAME) $(INCLUDE) LDFLAGS="${LIB_FLAGS} -fopenmp" CXXFLAGS="-fopenmp -std=c++11 -fPIC" 	

matlab: $(LIBNAME) $(MATLAB_MEX)
		


R: $(LIBNAME)
	cp $(LIBNAME) Rpackage/bin/	
	cp include/NetDECODE.h Rpackage/inst/include	
	R CMD INSTALL --preclean Rpackage/

				
.PHONY: clean ar

clean:
	rm -f $(PROGRAM)  src/*.o src/*~ include/*~ *~ $(LIBNAME) Rpackage/bin/$(LIBNAME) matlab/*.mexa64

install:
	cp $(LIBNAME) /usr/lib/
	cp include/decode.h /usr/include
	@echo done

archive:
	make clean
	tar -czvf ../$(PROGRAM)"(`date`)".tar.gz *

message:
	@echo "$(PROGRAM) library have been created"	
