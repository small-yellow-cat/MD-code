
a.out: read100.o te.o cut.o
	ifort -o a.out te.o read100.o cut.o


cut.o cut.mod: cut.f90
	ifort -c cut.f90


te.o te.mod: te.f90 cut.mod
	ifort -c te.f90

read100.o: read100.f90 te.mod
	ifort -c read100.f90

clean:
	rm  a.out  cut.o cut.mod te.o read100.o te.mod
