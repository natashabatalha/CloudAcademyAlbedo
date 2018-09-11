LDFLAGS= 

OBJECTS = pangle_data.o egpalb.o subs_egp.o raman2.o trist.o

#FFLAGS= -tpp7 -extend_source -w90 -w95 -r8 -132
#F77 = ifort
#FFLAGS=  -frecord-marker=4  -ffixed-line-length-132 -O3 -fdefault-real-8 -fdefault-double-8
FFLAGS=  -frecord-marker=4  -ffixed-line-length-132 -O3 -fdefault-real-8 -fdefault-double-8 -fno-align-commons
F77=gfortran

.SUFFIXES: .o .f .c

run/geom: $(OBJECTS)
	$(F77) $(FFLAGS) -o $@ $(OBJECTS)

.f.o:
	$(F77) $(FFLAGS) -c $<

pangle_data.o: pangle_data.f
egpalb.o: egpalb.f 
subs_egp.o: subs_egp.f 
raman2.o: raman2.f 
trist.o: trist.f

clean: 
	rm pangle_data.o egpalb.o subs_egp.o raman2.o trist.o
