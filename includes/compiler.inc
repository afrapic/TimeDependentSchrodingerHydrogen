#-----------------------------------------------------------------------
# 	Makefile for Intel Fortran compiler with MPI
#	Linux architecture
#-----------------------------------------------------------------------
#Suffixes to the source files (all posibilities, F77 or F90)
.SUFFIXES: .f .f90 .F .o .mod

#-----------------------------------------------------------------------
# fortran compiler and flags
#-----------------------------------------------------------------------

#sequential
FC=ifort
# Extensive runtime checks
FFLAGS = -save -mcmodel=large
MODULES = -module $(MDIR)

#compiling commands

$(ODIR)/%.o : %.f
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o : %.for
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(MODULES) -I $(INCDIR) -c $< -o $(ODIR)/$*.o

