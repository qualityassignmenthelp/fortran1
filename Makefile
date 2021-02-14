compiler = intel

###=====================================================================================================================
ifeq ($(compiler), intel)
    FC = ifort

#    SPEC_FOPT1 = -mp1 -pc80 -O0 -debug all -check all -warn all
    SPEC_FOPT1 = -fast

    SPEC_FOPT = $(SPEC_FOPT1) -traceback -assume buffered_io -static-intel
endif

ifeq ($(compiler), gnu)
    FC = gfortran

#    SPEC_FOPT1 = -O0 -Wall -fcheck=all -Wno-unused-dummy-argument
    SPEC_FOPT1 = -Ofast


    SPEC_FOPT = $(SPEC_FOPT1) -fbacktrace -ffree-line-length-none -fimplicit-none
endif

###=====================================================================================================================
MAKE_FILE  = Makefile

SRC_MAIN   = solution.f90
EXE        = $(SRC_MAIN:.f90=.exe)

$(EXE): $(SRC_MAIN) $(MAKE_FILE)
	$(FC) -o $(EXE)  $(SPEC_FOPT) $(SRC_MAIN)
#


