# Compilador Fortran
FC = gfortran

# Arquivo fonte: (módulos)
COMMON_SOURCES = quant_module.f90

# Lista de arquivos objeto
OBJECTS = $(COMMON_SOURCES:.f90=.o) $(S:.f90=.o)

# Flags para vincular as bibliotecas BLAS e LAPACK
LFLAGS = -lblas -llapack

EX = a.out

# Comando para compilar
compile: $(EX)

$(EX): $(OBJECTS)12
	$(FC) $^ -o $@ $(LFLAGS)

# Regra para criar arquivos objeto
%.o: %.f90
	$(FC) -c $< -o $@

# Comando para limpar (opcional)
clean:
	rm -f $(EX) $(OBJECTS)

# Para compilar: make S=filename.f90 EX=name_exec.out
