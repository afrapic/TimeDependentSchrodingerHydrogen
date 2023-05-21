#-----------------------------------------------------------------------
#  Local Makefile for stages
#  For Intel Fortran Composer and MKL Library
#-----------------------------------------------------------------------
SHELL   = /bin/bash
#
#Includes files for directories and compiler
#
include ../../includes/root.inc
include $(LIDIR)/dirs.inc
include $(LIDIR)/compiler.inc
#
#Points path to common and modules source folders
#
COMMON = ../common
MOD = ../modules
#
# Include file for path to Lapack+Blas and ScaLapack
#
include $(LIDIR)/libs.inc
#
# Stage name (local folder name)
#
NAME   = $(notdir $(PWD))
#
# Specifies where the compiled (.mod) module files will be searched
#
INC = -module $(MDIR) 
#
# ALL libraries
#
LIB = $(LIBS) 
#
#searches the directories listed in VPATH for a files, besides the local folder
#
VPATH = $(MOD) $(COMMON) 
#
#Local objects
#
LSRC   :=$(wildcard  *.f *.for *.f90)
LSRCS   = $(notdir $(LSRC))
LBASENAME = $(basename $(LSRCS))
LOBJS     = $(addsuffix .o, $(LBASENAME))
LOBJECTS  = $(addprefix $(ODIR)/, $(LOBJS))
#
#Common objects
#
CSRC   :=$(wildcard  $(COMMON)/*.f $(COMMON)/*.for $(COMMON)/*.f90)
CSRCS   = $(notdir $(CSRC))
CBASENAME = $(basename $(CSRCS))
COBJS     = $(addsuffix .o, $(CBASENAME))
COBJECTS  = $(addprefix $(ODIR)/, $(COBJS))
#
#Module objects
#
MSRC   :=$(wildcard  $(MOD)/*.f $(MOD)/*.for $(MOD)/*.f90)
MSRCS   = $(notdir $(MSRC))
MBASENAME = $(basename $(MSRCS))
MOBJS     = $(addsuffix .o, $(MBASENAME))
MOBJECTS  = $(addprefix $(ODIR)/, $(MOBJS))
#
#ALL objects (.o compiled objects)
#
OBJECTS = $(MOBJECTS) $(COBJECTS) $(LOBJECTS)
#
#Name and location of the executable 
#
EXE = $(STHOME)/$(NAME).x
#
#Directives to build this executable
#
all:  exe

exe: $(OBJECTS) 
	$(FC) $(FFLAGS) $(OBJECTS) $(LIB) -o $(EXE)

clean: 
	rm -f *.*~ makefile~
