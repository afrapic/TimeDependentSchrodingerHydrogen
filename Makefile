#----------------------------------------------------------------------------------------------------------
# 	Driver makefile for the code
#	The individual makefiles with stage specific variables are
#	inside each corresponding folder stage in src/(stage name)
# 	This model makefile works for compilation with Intel MKL Composer XE 2011 + Open MPI 1.4.5
#	Model Makefile (and makefiles in each stage), and includes provided by F Colavecchia,
#	modified by A L Frapiccini
#	report link or compilation errors to afrapic@uns.edu.ar
#----------------------------------------------------------------------------------------------------------
SHELL = /bin/bash
PROGNAME = TDFullGS
ALLFILES = src includes Makefile *.input 
SCRATCHDIR=scratch
VER=1.3

include includes/libs.inc
include includes/compiler.inc

#Tests whether the SCRATCHDIR and MKLROOT directories provided exists
locdir = $(shell test -d $(SCRATCHDIR) && echo 'ok' || echo 'no')
loclib = $(shell test -d $(MKLROOT) && echo 'ok' || echo 'no')

install: deploy testcomp testlib testdir

#Creates necessary directories, stores in includes/root.inc
deploy:
	@echo -e "******** Creating root includes at" $(PWD)
	@echo "STHOME = " $(PWD) > includes/root.inc
	@echo "LIDIR  = " $(PWD)/includes >> includes/root.inc
	@echo "SCRATCHDIR  = " $(SCRATCHDIR) >> includes/root.inc
	@cat includes/root.inc
	@echo -e "******** Creating directories objs modules and output"
	@mkdir -p  objs modules output scratch

#Prints message to warn about the compiler requiered
testcomp:
	@echo -e "******** Compiler requiered is FC=" $(FC) " check if available"

#Prints message if the location of library ScaLapack is ok or not	
testlib:
ifeq ($(loclib), ok)
	@echo -e "The location of MKLROOT for Lapack linking is OK"
else
	@echo -e "******** The location MKLROOT=' $(MKLROOT) ' in includes/libs.inc for MKLROOT for Lapack linking is NOT CORRECT"
endif

#Prints message if the location for SCRATCHDIR is ok or not. If correct, creates scratchdir_module with the fortran variable requiered
testdir: 
ifeq ($(locdir), ok)
	@echo -e '******** The location of SCRATCHDIR is OK:' $(SCRATCHDIR)
	$(shell make  filedir)
else
	@echo -e '******** The location SCRATCHDIR=' $(SCRATCHDIR) ' does NOT EXISTS'
endif

#Creates fortran module to pass the location of scratch as variable
filedir:
	@echo "module scratchdir_module" > src/modules/scratchdir_module.f90 
	@echo "implicit none" >> src/modules/scratchdir_module.f90 
	@echo "character*40,parameter :: dir_scratch='"$(SCRATCHDIR)"/'" >> src/modules/scratchdir_module.f90
	@echo "end module" >> src/modules/scratchdir_module.f90

target:
	sudo ./includes/script.sh $(INPUT)

#Executables
all: allclean basis matrices bstates tdse eigval spect

basis:
	(cd src/basis ; $(MAKE) exe) 

matrices:
	(cd src/matrices ; $(MAKE) exe)

bstates:
	(cd src/bstates ; $(MAKE) exe)  

tdse:
	(cd src/tdse ; $(MAKE) exe) 

eigval:
	(cd src/eigval ; $(MAKE) exe) 

spect:
	(cd src/spect ; $(MAKE) exe) 


#Cleans the object and modules files in objs and modules folders. Erases temporary files due to modification of source files
allclean:
	rm -f objs/* modules/* makefile~ src/common/*.*~ src/modules/*.*~ includes/*.*~ *~
	(cd src/basis ; $(MAKE) clean)
	(cd src/matrices ; $(MAKE) clean)
	(cd src/bstates ; $(MAKE) clean)
	(cd src/tdse ; $(MAKE) clean)
	(cd src/eigval ; $(MAKE) clean)
	(cd src/spect ; $(MAKE) clean)

#Cleans the output folder 
outclean:
	rm -f output/*.* slurm-* nohup* st*.o* st*.e*

#Cleans the scratch folder
scratchclean:
	ionice -c 3 rm -f $(SCRATCHDIR)/*.*

#Makes a tar.gz compressed folder with the sources, includes, slurm, Makefile and input files
zip:
	@$(MAKE) allclean ; mkdir $(PROGNAME)-$(VER) ; cp -pR $(ALLFILES) $(PROGNAME)-$(VER) ; \
	tar -czvf $(PROGNAME)-$(VER).tar.gz $(PROGNAME)-$(VER); rm -rf $(PROGNAME)-$(VER)

