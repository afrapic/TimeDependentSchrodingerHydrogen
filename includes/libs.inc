#
#-----------------------------------------------------------------------
# the following lines specify the position of SCALAPACK and BLACS
# with OPENMPI in the clusters hmem, lemaitre2 and manneback.
# It works for compilation with Intel MKL Composer XE 2011 + Open MPI 1.4.5
#-----------------------------------------------------------------------
#
HOST=$(shell hostname)

ifeq ($(HOST),ubuntu-VirtualBox)
MKLROOT = /opt/intel/compilers_and_libraries_2020.1.217/linux/mkl
LAPACK_LIB = -L$(MKLROOT)/include -mkl=sequential -lpthread -lm -ldl
endif

ifeq ($(HOST),ihpc01.ifisur-conicet.gob.ar)
MKLROOT=/share/localhost.apps/intel/mkl
LAPACK_LIB=  -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),afrapic-oficina)
MKLROOT=/opt/intel/composerxe-2011.4.191/mkl
LAPACK_LIB=  -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),afrapic)
MKLROOT=/opt/intel/composer_xe_2011_sp1/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),manneback.cism.ucl.ac.be)
MKLROOT=/opt/intel/compilerpro-12.0.0.084/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),usuario-Z68MX-UD2H-B3)
MKLROOT=/opt/intel/composer_xe_2011_sp1/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),ihpc01.ifisur-conicet.gob.ar)
MKLROOT=/opt/intel/composer_xe_2013/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),pauli.pamo.ucl.ac.be)
MKLROOT=/opt/intel/compilerpro-12.0.0.084/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),pauli03.pamo.ucl.ac.be)
MKLROOT=/opt/intel/compilerpro-12.0.0.084/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

ifeq ($(HOST),pauli02.pamo.ucl.ac.be)
MKLROOT=/opt/intel/compilerpro-12.0.0.084/mkl
LAPACK_LIB= -L$(MKLROOT)/include -mkl=sequential
endif

LIBS = $(LAPACK_LIB) 
