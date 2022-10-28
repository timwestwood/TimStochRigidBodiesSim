LIBS = -larmadillo -lopenblas -llapack -DARMA_DONT_USE_WRAPPER -L/home/taw16/sprng/lib -I/home/taw16/sprng/include -lsprng -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm
FASTOPTS = -DARMA_NO_DEBUG -O3 -w
STD = -std=c++11
CPP = rcross.cpp parallel_random_numbers.cpp quaternion.cpp particle.cpp linked_list.cpp collection_of_rigid_bodies.cpp fft_solver.cpp fcm_fluid_solver.cpp gmres_solver.cpp main.cpp

default:

	mpic++ -o rigid_bodies $(CPP) $(STD) $(LIBS)

clean:
	-rm rigid_bodies

fast:

	mpic++ -o rigid_bodies $(CPP) $(STD) $(LIBS) $(FASTOPTS)
