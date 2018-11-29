
LIB = -lm 

CC = gcc
CFLAGS = -O2
EXEC = massfunction.x

OBJ_HOD = header.o main.o utility.o sigmac.o transfnc.o transfunc_file.o \
	nonlinear_power_spectrum.o \
	input_params.o \
	tf_eisenstein_hu.o \
	dark_matter_statistics.o halo_mass_function.o mstar.o growthfactor.o

OBJ_NR = nrutil.o qromo.o midpnt.o midinf.o polint.o splint.o spline.o \
	qtrap.o zbrent.o trapzd.o

OBJS01 = $(OBJ_HOD) $(OBJ_NR)

$(EXEC): $(OBJS01)
	$(CC) -o $@ $(OBJS01) $(LIB)
$(OBJS01):	header.h nrutil.h


clean:
	rm -f *.o
