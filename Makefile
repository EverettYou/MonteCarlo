EXE = MC
OBJ = ./build
FC = gfortran
TOUCH = touch
RM = rm
CFLAGS = -fimplicit-none -ffree-line-length-none
LFLAGS = -framework Accelerate 

all: $(EXE)
$(EXE): $(OBJ)/$(EXE).o
	@$(FC) -o $(EXE) $(OBJ)/$(EXE).o $(LFLAGS) -I $(OBJ)
$(OBJ)/$(EXE).o: $(EXE).f90
	@$(FC) -c $(EXE).f90 -o $(OBJ)/$(EXE).o $(CFLAGS) -J $(OBJ)
clean: 
	$(RM) build/*
	$(RM) $(EXE)
