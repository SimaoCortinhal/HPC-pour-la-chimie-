program lecteur_mol2
  implicit none

  ! type pour stocker les informations d'un atome
  type :: info_atome
    character(len=8) :: element    ! Le symbole de l'atome, ex: N, C, O, etc.
    real :: x, y, z     ! Coordonnées
  end type info_atome

  integer           :: statut_IO, unite_num, nb_atomes, id_atome, idx
  character(len=256):: ligne, nom_fichier
  logical           :: dans_section_atome
  type(info_atome), allocatable :: atomes(:)

  ! Variables temporaires pour la lecture
  character(len=8)  :: nom_atome_tmp
  character(len=8)  :: type_atome_tmp
  real :: x_tmp, y_tmp, z_tmp
  integer :: entier_a_ignorer
  character(len=20) :: chaine_a_ignorer
  integer           :: pos_point

  dans_section_atome = .false.
  nb_atomes = 0

  ! nom du fichier en argument
  call get_command_argument(1, nom_fichier)
  if (len_trim(nom_fichier) == 0) then
    print *, "Usage  : lecteur_mol2 <nom_du_fichier.mol2>"
    stop
  end if

  ! Ouverture du fichier mol2
  open(newunit=unite_num, file=trim(nom_fichier), status="old", action="read", iostat=statut_IO)
  if (statut_IO /= 0) then
    print *, "Erreur lors de l'ouverture du fichier ", trim(nom_fichier)
    stop
  end if

  ! Parcours du fichier pour trouver les sections MOLECULE et ATOM
  do
    read(unite_num,'(A)', iostat=statut_IO) ligne
    if (statut_IO /= 0) exit

    if (trim(ligne) == "@<TRIPOS>MOLECULE") then
      read(unite_num, '(A)', iostat=statut_IO) ligne
      read(unite_num, *) nb_atomes ! Nombre d'atomes
      cycle
    end if

    if (trim(ligne) == "@<TRIPOS>ATOM") then
      dans_section_atome = .true.
      exit
    end if
  end do

  if (nb_atomes <= 0) then
    print *, "Erreur : nombre d'atomes non trouvé dans la section MOLECULE"
    stop
  end if

  allocate(atomes(nb_atomes))

  ! Lecture de chaque atome
  do idx = 1, nb_atomes
    read(unite_num, '(A)', iostat=statut_IO) ligne
    if (statut_IO /= 0) exit

    read(ligne,*) id_atome, nom_atome_tmp, x_tmp, y_tmp, z_tmp, &
                   type_atome_tmp, entier_a_ignorer, chaine_a_ignorer

    ! On retire tout ce qui est après le point, s'il y a un point
    pos_point = index(type_atome_tmp, ".")
    if (pos_point >0) then
       type_atome_tmp = type_atome_tmp(1:pos_point-1)
    end if

    atomes(idx)%element = trim(type_atome_tmp)
    atomes(idx)%x = x_tmp
    atomes(idx)%y = y_tmp
    atomes(idx)%z = z_tmp
  end do
  close(unite_num)

  ! Affichage de contrôle
  print *, "Nb atomes lus =", nb_atomes
  do idx = 1, nb_atomes
    print '(I3,1X,A8,3F10.4)', idx, trim(atomes(idx)%element), &
          atomes(idx)%x, atomes(idx)%y, atomes(idx)%z
  end do
  deallocate(atomes)

end program lecteur_mol2
