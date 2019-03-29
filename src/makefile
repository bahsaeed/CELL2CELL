.SUFFIXES: .f90
################# Compiler flag #####################
FC=gfortran        # IBM-FORTRAN compiler! If you are 
                # using another comppiler, FC should
                # be editted to the availabe one.
################ Processor options ##################
PREFLAGS=
################ Compiler options ###################
FFLAGS= -O3 
################ Link flags #########################
LFLAGS=-L/Users/admin/Documents/Utility/lapack-3.6.1 -llapack -lrefblas
################ Directory path for EXE file ########
LPATH= ../bin 
                # You need to change the above path 
                # to your own path address
#####################################################

################ 
WKDIR= .
################ Profram files ######################
MODUL= $(WKDIR)/fun.f90    
INPUT=  
PUBLIC=
OUTPUT=
TRANS= 
MAIN=  $(WKDIR)/cell2cell.f90
################ Special files ######################
OTHER= $(WKDIR)/INPUT       
#####################################################
ALL= $(MODUL) $(INPUT) $(PUBLIC) $(OUTPUT) $(TRANS) $(MAIN)
################ Execuation #########################
.f90.o:
	$(FC) $(PREFLAGS) $(FFLAGS) -c $*.f90
#	mv $(*F).o $(*D)
#	mv $(*F).mod $(*D)
 ALLOBJ=${ALL:.f90=.o}
 ALLMOD=${ALL:.f90=.mod}
################ Dependencies #######################

############# Create Exe file #######################
cell2cell:  $(ALLOBJ)
	$(FC) $(ALLOBJ) $(LFLAGS)  $(FFLAGS) -o $@
	chmod ugo+r *
	chmod ugo+x cell2cell
############ Clean ##################################
clean:
	rm -f $(WKDIR)/*.mod $(WKDIR)/*.o $(WKDIR)/*.omp \
              $(WKDIR)/*.lst $(WKDIR)/*.log $(WKDIR)/*.L cell2cell
