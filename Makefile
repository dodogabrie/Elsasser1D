MAIN = main.f90
MDIR = ./modules/
MODULES_names := physicals.f90 numericals.f90 lele_filter.f90 derivative.f90 IO.f90
MODULES := $(addprefix $(MDIR), $(MODULES_names))
COMPILE_BASE := gfortran -J $(MDIR) $(MODULES)

help:
	@echo "Available tasks:"
	@echo "compile  -> Compile fortran main and modules"
	@echo "run      -> Run the compiled codes"
	@echo "all      -> Call compile, run and plot results"
	@echo "clean    -> Remove *.mod *.so *.c *.o *.html build"

all: compile run plot

compile:
	$(COMPILE_BASE) $(MAIN) -Ofast -o run.out

.PHONY: run
run:
	./run.out

plot:
	python plots/main.py

# Phony targets for cleanup and similar uses
.PHONY: clean
clean:
	rm -rf *~ *.mod $(MDIR)/*.mod *.out *.o *.html build __pycache__

clean-data:
	rm -rf data/nlin/* data/line/*

test-penta: 
	$(COMPILE_BASE) tests/test_penta.f90 -llapack -lblas -o tests/test_penta.out
	clear
	./tests/test_penta.out
test-tridiag: 
	$(COMPILE_BASE) tests/test_tridiag.f90 -llapack -lblas -o tests/test_tridiag.out
	clear
	./tests/test_tridiag.out
test-derivative: 
	$(COMPILE_BASE) tests/test_derivative.f90 -o tests/test_derivative.out
	clear
	./tests/test_derivative.out
	python tests/plot_der.py
test-filters: 
	$(COMPILE_BASE) tests/test_lele.f90 -o tests/test_lele.out
	clear
	./tests/test_lele.out
	python tests/plot_lele.py
	
	

# 2D extension (Orszag-Tang MHD turbulence, see main2d.f90)
compile2d:
	gfortran -J $(MDIR) $(MDIR)physicals.f90 $(MDIR)numericals.f90 $(MDIR)lele_filter.f90 $(MDIR)derivative.f90 $(MDIR)derivative2d.f90 $(MDIR)poisson2d.f90 main2d.f90 -Ofast -fno-stack-arrays -o run2d.out
run2d:
	mkdir -p data2d && ./run2d.out
