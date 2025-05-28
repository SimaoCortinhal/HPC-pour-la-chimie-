# CHPS1002 – MOL2, Covalence and co.
**Simao CORTINHAL — 2024/2025**

## Présentation

Ce projet a pour objectif la reconstruction de fichiers moléculaires **SYBYL mol2** valides à partir d’informations atomiques, l’identification des liaisons chimiques, puis la mise en place d’un algorithme génétique pour l’amarrage moléculaire (docking). Le tout développé en Fortran 90, avec des versions séquentielles et parallèles (OpenMP).

## Organisation du projet

Le projet est structuré en deux grandes parties :
- **partie_1/** : Reconnaissance atomique, détermination des liaisons (topologie), génération de fichiers MOL2 complets et optimisation/parallélisation.
- **partie_2/** : Implémentation d’un algorithme génétique d’amarrage moléculaire basé sur les fichiers générés, avec analyse de performance et tests.

Chaque partie possède ses propres scripts Fortran, Makefile, exemples de fichiers mol2, et fichiers de paramètres.

---

## Détail des étapes

### Partie 1 – 

1. **Lecture d’un fichier mol2**
   - `lecteur_mol2.f90` : Code qui lit un fichier mol2 passé en argument, et affiche pour chaque atome son symbole chimique et ses coordonnées 3D.
   - Utilisation de fichiers `.mol2` ne contenant que les coordonnées atomiques pour les tests.

2. **Chargement des rayons de covalence**
   - `chargeur_covalence.f90` : Permet de lire le fichier de rayons de covalence, d’interroger ces valeurs pour différents éléments, et de vérifier leur bonne intégration dans le programme.

3. **Détection de la topologie (liaisons chimiques)**
   - `affiche_topologie.f90` : Combine les données des rayons de covalence et des coordonnées atomiques pour proposer une topologie complète de la molécule, en testant différents seuils de tolérance (delta) pour relier les atomes.
   - Génère la liste des liaisons (simple, double, triple) pour chaque paire d’atomes.

4. **Génération de fichiers mol2 complets**
   - `complete_mol2.f90` : Ajoute la section BOND dans un nouveau fichier mol2 (sans écraser l’original), pour produire une molécule prête à être visualisée (ex : sous UCSF Chimera, VMD…).

5. **Optimisation et parallélisation**
   - `para_complete_mol2.f90` : Version parallèle (OpenMP) de la recherche de topologie, accélérant le traitement sur les grosses molécules.

---

### Partie 2 

6. **Amarrage moléculaire par algorithme génétique**
   - `genetic_docking.f90` : Implémente l’algorithme génétique pour faire évoluer des solutions représentant la position/orientation du ligand.
   - Les solutions sont évaluées sur le **nombre de ponts hydrogène** formés avec le site cible (distance 2.2–4 Å, angle 90–150°).
   - Inclut génération aléatoire initiale, sélection, recombinaison, mutation, et suivi de la population.

7. **Parallélisation**
   - `para_genetic_docking.f90` : Parallélise certaines étapes de l’algorithme génétique (ex : évaluation de la population, recombinaison, etc.) via OpenMP.

8. **Protocole de test et résultats**
   - Exécution répétée 10 fois avec population=100, 500 générations max, taux de recombinaison=95%, mutation=5%.
   - Sauvegarde de la population finale (mol2) et statistiques d’évolution (CSV).
   - Rapport des gains de performance obtenus avec la parallélisation.

---

## Structure des dossiers

- **partie_1/**
  - Fichiers `.mol2` exemples, rayons de covalence, scripts Fortran (lecteur, chargeur, topologie, complétion, parallèle), Makefile, etc.
- **partie_2/**
  - Scripts Fortran pour le docking, version parallèle, Makefile, fichiers ligand/site, fichiers résultats.

- **CHPS1002_Simao_Cortinhal.pdf** : Rapport technique détaillé, réponses aux questions, organisation du travail et difficultés rencontrées.
