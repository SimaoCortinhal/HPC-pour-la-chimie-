program chargeur_covalence
  implicit none

  ! Paramètre pour double précision
  integer, parameter :: dp = kind(1.0d0)

  ! Type pour stocker les rayons de covalence en double précision
  type :: info_covalence
    character(len=2) :: element            ! Symbole chimique
    real(kind=dp)    :: rayon_simple
    real(kind=dp)    :: rayon_double
    real(kind=dp)    :: rayon_triple
  end type info_covalence

  integer                    :: statut_IO, unite_num
  integer                    :: nb_elements, idx, jdx, compteur
  integer                    :: statut_IO_temp
  character(len=256)         :: nom_fichier, ligne, ligne_trimmee
  type(info_covalence), allocatable :: covalences(:)
  character(len=2)           :: element_str
  character(len=3)           :: rayon_str
  real(kind=dp)              :: val_tmp
  logical                    :: trouve
  character(len=10)          :: element_saisi

  ! Récupération du nom du fichier en argument
  call get_command_argument(1, nom_fichier)
  if (len_trim(nom_fichier) == 0) then
    print *, "Usage: chargeur_covalence <fichier_rayons.txt>"
    stop
  end if

  ! Ouverture du fichier de rayons
  open(newunit=unite_num, file=trim(nom_fichier), status="old", action="read", iostat=statut_IO)
  if (statut_IO /= 0) then
    print *, "Erreur lors de l'ouverture du fichier :", trim(nom_fichier)
    stop
  end if

  ! compter les lignes exploitables
  nb_elements = 0
  do
    read(unite_num, '(A)', iostat=statut_IO) ligne
    if (statut_IO < 0) exit                   ! fin de fichier
    if (statut_IO > 0) then
      print *, "Erreur de lecture (1ère passe)."
      stop
    end if

    ligne_trimmee = trim(ligne)
    if (len_trim(ligne_trimmee) == 0) cycle   ! ignorer ligne vide
    if (ligne_trimmee(1:1) == "#") cycle      ! ignorer commentaire
    if (len(ligne_trimmee) < 5) cycle         ! trop court
    nb_elements = nb_elements + 1
  end do

  rewind(unite_num)
  if (nb_elements == 0) then
    print *, "Fichier vide ou aucune ligne exploitable."
    stop
  end if

  ! Allocation du tableau de données
  allocate(covalences(nb_elements))
  compteur = 0

  ! lecture et parsing des valeurs
  do
    read(unite_num, '(A)', iostat=statut_IO) ligne
    if (statut_IO < 0) exit                   ! fin de fichier
    if (statut_IO > 0) then
      print *, "Erreur de lecture (2ème passe)."
      stop
    end if

    ligne_trimmee = trim(ligne)
    if (len_trim(ligne_trimmee) == 0) cycle
    if (ligne_trimmee(1:1) == "#") cycle
    if (len(ligne_trimmee) < 5) cycle

    if (compteur == nb_elements) exit
    compteur = compteur + 1

    !  Symbole (colonnes 4-5)
    element_str = ligne_trimmee(4:5)

    ! Rayon simple (colonnes 6-8) -> division par 100 pour Å
    if (len(ligne_trimmee) >= 8) then
      rayon_str = ligne_trimmee(6:8)
      read(rayon_str, *, iostat=statut_IO_temp) val_tmp
      if (statut_IO_temp == 0) then
        covalences(compteur)%rayon_simple = val_tmp / 100.0_dp
      else
        covalences(compteur)%rayon_simple = 0.0_dp
      end if
    else
      covalences(compteur)%rayon_simple = 0.0_dp
    end if

    ! Rayon double (colonnes 9-11)
    if (len(ligne_trimmee) >= 11) then
      rayon_str = ligne_trimmee(9:11)
      read(rayon_str, *, iostat=statut_IO_temp) val_tmp
      if (statut_IO_temp == 0) then
        covalences(compteur)%rayon_double = val_tmp / 100.0_dp
      else
        covalences(compteur)%rayon_double = 0.0_dp
      end if
    else
      covalences(compteur)%rayon_double = 0.0_dp
    end if

    ! Rayon triple (colonnes 12-14)
    if (len(ligne_trimmee) >= 14) then
      rayon_str = ligne_trimmee(12:14)
      read(rayon_str, *, iostat=statut_IO_temp) val_tmp
      if (statut_IO_temp == 0) then
        covalences(compteur)%rayon_triple = val_tmp / 100.0_dp
      else
        covalences(compteur)%rayon_triple = 0.0_dp
      end if
    else
      covalences(compteur)%rayon_triple = 0.0_dp
    end if

    covalences(compteur)%element = trim(element_str)
  end do

  close(unite_num)
  if (compteur < nb_elements) nb_elements = compteur

  ! Interrogation utilisateur
  do idx = 1, 3
    print *, "Entrez un symbole d'élément chimique (ex: H, He, Li...) :"
    read(*, *) element_saisi
    trouve = .false.
    do jdx = 1, nb_elements
      if (trim(element_saisi) == trim(covalences(jdx)%element)) then
        print *, "Rayons de covalence pour", trim(covalences(jdx)%element), "(en Å) :"
        print '(A, F6.2)', "  Simple :",  covalences(jdx)%rayon_simple
        print '(A, F6.2)', "  Double :",  covalences(jdx)%rayon_double
        print '(A, F6.2)', "  Triple :",  covalences(jdx)%rayon_triple
        trouve = .true.
        exit
      end if
    end do
    if (.not. trouve) then
      print *, "Élément inconnu."
    end if
  end do

  ! Libération de la mémoire
  deallocate(covalences)

end program chargeur_covalence
