.SUFFIXES:
.SUFFIXES: .f .f90 .o


#sgi
#FL = f90
#CMPLFLG = -c -O3
#LNKFLG = -o $(TARGET)  

#alpha
#FL = fort
#CMPLFLG = -c -fast 
#LNKFLG = -o $(TARGET) -lcxml
#-warn argument_checking -C

#P4
#FL = ifort
#FL = /opt/mpich/intel/fc/mpich-1.2.7p1/bin/mpif90

#Raj
FL = gfortran -fbounds-check -ffree-line-length-none

CMPLFLG = -c -O3 #-lmpi
#CMPLFLG = -C  
LNKFLG = -o $(TARGET) #-L/opt/intel/l_mkl_enh_p_9.1.021/lib/32/ -Vaxlib -lguide -lmkl_lapack -lmkl_ia32  -i-static
TARGET = IMLS.x

OBJS = dynamic_parameters_4D_vdW.o index.o PES_CO-CO.o PES_4D_subroutine_vdW.o diag.o
$(TARGET) : $(OBJS)
	$(FL) $(OBJS) $(LNKFLG)

# this is a freeform

$(OBJS) :

#equivalent to
#%.f : %.o
.f90.o:
	$(FL) $(CMPLFLG) $<

.f.o:
	$(FL) $(CMPLFLG) $<
