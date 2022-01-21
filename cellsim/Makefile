FLAGS = -fPIC

# the python interface through swig
PYTHONI =  `pkg-config python3 --cflags-only-I`
PYTHONL = -Xlinker -export-dynamic 
CFLAGS = -I/usr/include  -Wall -Wextra -march=native  -O3
LIBFLAGS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm
# default super-target
#
all: POP3p2x2  Rm 

Rm: 
	g++ $(CFLAGS)  -march=native -O3 -fPIC -c $@.c -o $@.o
	swig -c++ -python -o $@_wrap.cxx $@.i 
	g++ $(FLAGS) $(PYTHONI)  -c $@_wrap.cxx -o $@_wrap.o
	g++ $@.o $@_wrap.o $(PYTHONL) $(LIBFLAGS) -shared  -o _$@.so


POP3p2x2: 
	g++ $(CFLAGS)  -march=native -O3 -fPIC -c $@.c -o $@.o
	swig -c++ -python -o $@_wrap.cxx $@.i 
	g++ $(FLAGS) $(PYTHONI)  -c $@_wrap.cxx -o $@_wrap.o
	g++ $@.o $@_wrap.o $(PYTHONL) $(LIBFLAGS) -shared  -o _$@.so


