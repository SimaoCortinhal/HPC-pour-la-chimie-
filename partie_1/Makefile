# Définit le compilateur Fortran
CC      = gfortran
FFLAGS  = -fopenmp -O3

# Cibles par défaut
all: chargeur_covalence lecteur_mol2 affiche_topologie complete_mol2 genetic_docking

# ------------------------------------------------------------------
# Règles de construction des exécutables
# ------------------------------------------------------------------

chargeur_covalence: chargeur_covalence.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+

lecteur_mol2: lecteur_mol2.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+

affiche_topologie: affiche_topologie.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+

complete_mol2: complete_mol2.o
	@echo ""
	@echo "Linking $(@)"
	$(CC) $(FFLAGS) -o $@ $+

genetic_docking: genetic_docking.o
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
EXEC = chargeur_covalence \
       lecteur_mol2 \
       affiche_topologie \
       complete_mol2 \
       genetic_docking

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
