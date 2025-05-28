# Définit le compilateur Fortran
CC      = gfortran
FFLAGS  = -fopenmp -O3

# Cibles par défaut
all: genetic_docking para_genetic_docking

# ------------------------------------------------------------------
# Règles de construction des exécutables
# ------------------------------------------------------------------



genetic_docking: genetic_docking.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+

para_genetic_docking: para_genetic_docking.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+


# ------------------------------------------------------------------
# Règle générique pour compiler les .f90 en .o
# ------------------------------------------------------------------

%.o: %.f90
	@echo ""
	@echo "Generating $@"
	$(CC) $(FFLAGS) -c $<

# Liste des exécutables à générer
EXEC = genetic_docking para_genetic_docking

# ------------------------------------------------------------------
# Nettoyage
# ------------------------------------------------------------------

clean: clean_all
	@rm -rf *.o

clean_all: clean cleanSource
	@rm -rf $(EXEC)

cleanSource:
	@find . \( -name "*~" -o -name "*.old" -o -name "#*" \) -print -exec rm {} \;
	@find . -name "*.mod" -print -exec rm {} \;

.PHONY: all clean clean_all cleanSource
