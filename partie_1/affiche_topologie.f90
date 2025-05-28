program affiche_topologie
  implicit none

  
  !  Paramètre double précision 
  integer, parameter :: dp = kind(1.0d0)

  
  ! Type pour stocker les rayons de covalence en double précision
  type :: covalence_data
    character(len=2) :: element
    real(kind=dp)    :: rayon_simple
    real(kind=dp)    :: rayon_double
    real(kind=dp)    :: rayon_triple
  end type covalence_data

  ! type pour stocker les informations d'un atome
  type :: atom_info
    character(len=8) :: elem
    real(kind=dp)    :: x, y, z
  end type atom_info

  type :: bond
    integer          :: a1, a2
    character(len=10):: bond_type
  end type bond

  character(len=256)             :: fichier_mol2, fichier_covalence
  type(covalence_data), allocatable :: covalences(:)
  type(atom_info),       allocatable :: atomes(:)
  type(bond),            allocatable :: liaisons(:)
  integer,               allocatable :: connectes(:)

  integer :: nb_atomes, nb_covalence
  integer :: pourcentage_delta, idx, jdx, kdx, nombre_liaisons
  real(kind=dp) :: delta_frac
  logical :: tous_connectes
  real(kind=dp) :: distance_ij
  character(len=2) :: element_i, element_j

  
  ! Début du programme
  call get_command_argument(1, fichier_mol2)
  call get_command_argument(2, fichier_covalence)

  if (len_trim(fichier_mol2) == 0 .or. len_trim(fichier_covalence) == 0) then
    print *, "Usage: affiche_topologie <fichier.mol2> <fichier_rayons.txt>"
    stop
  end if

  ! Lecture du fichier de covalence 
  call read_covalence(fichier_covalence, covalences, nb_covalence)

  ! Affichage des covalences lues
  print *, "--------------------------------------------------"
  print *, "Liste des éléments dans le fichier de covalence :"
  print *, "(Symbole, rayon_simple, rayon_double, rayon_triple) en Angstroms"
  do idx = 1, nb_covalence
    write(*,'(A3,2X,F6.3,2X,F6.3,2X,F6.3)') trim(covalences(idx)%element), &
         covalences(idx)%rayon_simple, covalences(idx)%rayon_double, covalences(idx)%rayon_triple
  end do
  print *, "--------------------------------------------------"

  ! lecture du fichier mol2
  call read_mol2(fichier_mol2, atomes, nb_atomes)

  ! Affichage des atomes lus
  print *, "Liste des atomes dans le fichier mol2 :"
  print *, "(Index, Symbole, x, y, z)"
  do idx = 1, nb_atomes
    write(*,'(I3,1X,A8,3F10.4)') idx, trim(atomes(idx)%elem), &
         atomes(idx)%x, atomes(idx)%y, atomes(idx)%z
  end do
  print *, "--------------------------------------------------"

  ! Vérifi que chaque atome a un rayon de covalence défini
  call check_elements(atomes, nb_atomes, covalences, nb_covalence)

  ! Boucle sur delta (10%..35% par pas de 5%)
  do pourcentage_delta = 10, 35, 5
    delta_frac = real(pourcentage_delta, kind=dp) / 100.0_dp

    ! On réinitialise la liste de liaisons
    if (allocated(liaisons)) deallocate(liaisons)
    allocate(liaisons(nb_atomes*(nb_atomes-1)/2))
    nombre_liaisons = 0

    do idx = 1, nb_atomes - 1
      do jdx = idx + 1, nb_atomes
        call calculate_distance(atomes, idx, jdx, distance_ij)

        ! Récupére les symboles
        element_i = trim(atomes(idx)%elem)
        element_j = trim(atomes(jdx)%elem)

        ! On teste triple, double, single
        if (try_bond_type('triple', element_i, element_j, distance_ij, delta_frac, &
             idx, jdx, covalences, nb_covalence, nombre_liaisons, liaisons)) cycle

        if (try_bond_type('double', element_i, element_j, distance_ij, delta_frac, &
             idx, jdx, covalences, nb_covalence, nombre_liaisons, liaisons)) cycle

        if (try_bond_type('simple', element_i, element_j, distance_ij, delta_frac, &
             idx, jdx, covalences, nb_covalence, nombre_liaisons, liaisons)) cycle
      end do
    end do

    call check_connectivity(nb_atomes, liaisons, nombre_liaisons, tous_connectes)

    if (tous_connectes) then
      print *, "Topologie trouvée avec un delta de", pourcentage_delta, "%"
      do kdx = 1, nombre_liaisons
        print '(A,I0,A,I0,A,A)', "Liaison entre atome ", liaisons(kdx)%a1, &
                                " et ", liaisons(kdx)%a2, " : ", trim(liaisons(kdx)%bond_type)
      end do
      stop
    end if

    deallocate(liaisons)
  end do

  print *, "Impossible de déterminer la topologie avec un delta maximum de 35%."

contains

  subroutine read_covalence(fichier, covs, nb_cov)
    character(len=*), intent(in) :: fichier
    type(covalence_data), allocatable, intent(out) :: covs(:)
    integer, intent(out) :: nb_cov

    integer :: unite, statut_IO, compteur, statut_temp
    character(len=256) :: ligne, ligne_trimmee
    character(len=2) :: element_str
    character(len=3) :: rayon_str
    real(kind=dp) :: val_tmp

    open(newunit=unite, file=fichier, status="old", action="read", iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture covalence :", trim(fichier)
      stop
    end if

    nb_cov = 0
    do
      read(unite,'(A)',iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      if (statut_IO > 0) then
        print *, "Erreur lecture (1ère passe) dans", trim(fichier)
        stop
      end if
      ligne_trimmee = trim(ligne)
      if (len_trim(ligne_trimmee)==0 .or. ligne_trimmee(1:1)=="#" .or. len(ligne_trimmee)<5) cycle
      nb_cov = nb_cov + 1
    end do

    rewind(unite)
    if (nb_cov == 0) then
      print *, "Aucune ligne exploitable dans", trim(fichier)
      stop
    end if

    allocate(covs(nb_cov))
    compteur = 0

    do
      read(unite,'(A)',iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      if (statut_IO > 0) then
        print *, "Erreur lecture (2ème passe) dans", trim(fichier)
        stop
      end if
      ligne_trimmee = trim(ligne)
      if (len_trim(ligne_trimmee)==0 .or. ligne_trimmee(1:1)=="#" .or. len(ligne_trimmee)<5) cycle

      if (compteur == nb_cov) exit
      compteur = compteur + 1

      element_str = ligne_trimmee(4:5)

      if (len(ligne_trimmee) >= 8) then
        rayon_str = ligne_trimmee(6:8)
        read(rayon_str,*,iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covs(compteur)%rayon_simple = val_tmp / 100.0_dp
        else
          covs(compteur)%rayon_simple = 0.0_dp
        end if
      else
        covs(compteur)%rayon_simple = 0.0_dp
      end if

      if (len(ligne_trimmee) >= 11) then
        rayon_str = ligne_trimmee(9:11)
        read(rayon_str,*,iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covs(compteur)%rayon_double = val_tmp / 100.0_dp
        else
          covs(compteur)%rayon_double = 0.0_dp
        end if
      else
        covs(compteur)%rayon_double = 0.0_dp
      end if

      if (len(ligne_trimmee) >= 14) then
        rayon_str = ligne_trimmee(12:14)
        read(rayon_str,*,iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covs(compteur)%rayon_triple = val_tmp / 100.0_dp
        else
          covs(compteur)%rayon_triple = 0.0_dp
        end if
      else
        covs(compteur)%rayon_triple = 0.0_dp
      end if

      covs(compteur)%element = trim(element_str)
    end do

    close(unite)
    if (compteur < nb_cov) nb_cov = compteur
  end subroutine read_covalence


  subroutine read_mol2(fichier, atomes, nb_atomes)
    character(len=*), intent(in) :: fichier
    type(atom_info), allocatable, intent(out) :: atomes(:)
    integer, intent(out) :: nb_atomes

    integer :: unite, statut_IO, idx, id_dummy
    character(len=256) :: ligne
    character(len=8)   :: nom_tmp, type_tmp
    real(kind=dp) :: x_tmp, y_tmp, z_tmp
    integer :: pos_point

    open(newunit=unite, file=fichier, status="old", action="read", iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture mol2 :", trim(fichier)
      stop
    end if

    nb_atomes = 0
    do
      read(unite,'(A)',iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      if (statut_IO > 0) then
        print *, "Erreur lecture (1ère passe mol2)."
        stop
      end if
      if (trim(ligne) == "@<TRIPOS>MOLECULE") then
        read(unite,'(A)',iostat=statut_IO) ligne
        if (statut_IO < 0 .or. statut_IO > 0) cycle
        read(unite,*) nb_atomes
        exit
      end if
    end do

    if (nb_atomes <= 0) then
      print *, "Impossible de lire le nombre d'atomes."
      stop
    end if

    do
      read(unite,'(A)',iostat=statut_IO) ligne
      if (statut_IO < 0) then
        print *, "Pas de section ATOM."
        stop
      end if
      if (statut_IO > 0) then
        print *, "Erreur lecture (recherche ATOM)."
        stop
      end if
      if (trim(ligne) == "@<TRIPOS>ATOM") exit
    end do

    allocate(atomes(nb_atomes))

    do idx = 1, nb_atomes
      read(unite,'(A)',iostat=statut_IO) ligne
      if (statut_IO < 0) then
        print *, "EOF inattendu."
        stop
      else if (statut_IO > 0) then
        print *, "Erreur lecture atomes."
        stop
      end if

      read(ligne,*,iostat=statut_IO) id_dummy, nom_tmp, x_tmp, y_tmp, z_tmp, type_tmp
      if (statut_IO /= 0) then
        print *, "Erreur parsing ATOM :"; print *, ligne; stop
      end if
      pos_point = index(type_tmp, '.')
      if (pos_point > 0) type_tmp = type_tmp(1:pos_point-1)
      call uppercase(type_tmp)

      atomes(idx)%elem = trim(type_tmp)
      atomes(idx)%x    = x_tmp
      atomes(idx)%y    = y_tmp
      atomes(idx)%z    = z_tmp
    end do

    close(unite)
  end subroutine read_mol2


  subroutine check_elements(atomes, n, covs, nb_cov)
    type(atom_info), intent(in) :: atomes(:)
    integer, intent(in) :: n, nb_cov
    type(covalence_data), intent(in) :: covs(:)

    integer :: i, j
    logical :: trouve
    do i = 1, n
      trouve = .false.
      do j = 1, nb_cov
        if (trim(atomes(i)%elem) == trim(covs(j)%element)) then
          trouve = .true.; exit
        end if
      end do
      if (.not. trouve) then
        print *, "Élément non trouvé :", trim(atomes(i)%elem)
        stop
      end if
    end do
  end subroutine check_elements


  subroutine calculate_distance(atomes, i, j, distance_ij)
    type(atom_info), intent(in) :: atomes(:)
    integer, intent(in)         :: i, j
    real(kind=dp), intent(out)  :: distance_ij

    real(kind=dp) :: dx, dy, dz
    dx = atomes(i)%x - atomes(j)%x
    dy = atomes(i)%y - atomes(j)%y
    dz = atomes(i)%z - atomes(j)%z
    distance_ij = sqrt(dx*dx + dy*dy + dz*dz)
  end subroutine calculate_distance

  logical function try_bond_type(type_liaison, e1, e2, distance_ij, delta_frac, &
                                i, j, covs, nb_cov, nombre_liaisons, liaisons)
    character(len=*), intent(in) :: type_liaison, e1, e2
    real(kind=dp), intent(in)    :: distance_ij, delta_frac
    integer, intent(in)          :: i, j, nb_cov
    type(covalence_data), intent(in) :: covs(:)
    integer, intent(inout)       :: nombre_liaisons
    type(bond), intent(inout)    :: liaisons(:)

    real(kind=dp) :: r1, r2, somme_r, borne_inf, borne_sup

    r1 = get_radius(e1, type_liaison, covs, nb_cov)
    r2 = get_radius(e2, type_liaison, covs, nb_cov)
    somme_r = r1 + r2

    if (somme_r <= 0.0_dp) then
      try_bond_type = .false.; return
    end if

    borne_inf = somme_r * (1.0_dp - delta_frac)
    borne_sup = somme_r * (1.0_dp + delta_frac)

    if (distance_ij >= borne_inf .and. distance_ij <= borne_sup) then
      nombre_liaisons = nombre_liaisons + 1
      liaisons(nombre_liaisons)%a1       = i
      liaisons(nombre_liaisons)%a2       = j
      liaisons(nombre_liaisons)%bond_type = type_liaison
      try_bond_type = .true.
    else
      try_bond_type = .false.
    end if
  end function try_bond_type

 
  real(kind=dp) function get_radius(element, type_liaison, covs, nb_cov)
    character(len=*), intent(in) :: element, type_liaison
    type(covalence_data), intent(in) :: covs(:)
    integer, intent(in) :: nb_cov
    integer :: idx_k

    get_radius = 0.0_dp
    do idx_k = 1, nb_cov
      if (trim(covs(idx_k)%element) == trim(element)) then
        select case(type_liaison)
        case('simple')
          get_radius = covs(idx_k)%rayon_simple
        case('double')
          get_radius = covs(idx_k)%rayon_double
        case('triple')
          get_radius = covs(idx_k)%rayon_triple
        end select
        return
      end if
    end do
  end function get_radius

  subroutine check_connectivity(nb_atomes, liaisons, nombre_liaisons, tous_connectes)
    integer, intent(in)         :: nb_atomes, nombre_liaisons
    type(bond), intent(in)      :: liaisons(:)
    logical, intent(out)        :: tous_connectes
    integer, allocatable :: connectes(:)
    integer :: idx_k

    allocate(connectes(nb_atomes))
    connectes = 0

    do idx_k = 1, nombre_liaisons
      connectes(liaisons(idx_k)%a1) = connectes(liaisons(idx_k)%a1) + 1
      connectes(liaisons(idx_k)%a2) = connectes(liaisons(idx_k)%a2) + 1
    end do

    tous_connectes = all(connectes > 0)
    deallocate(connectes)
  end subroutine check_connectivity

  subroutine uppercase(str)
    character(len=*), intent(inout) :: str
    integer :: i, c
    do i = 1, len(str)
      c = ichar(str(i:i))
      if (c >= ichar('a') .and. c <= ichar('z')) then
        str(i:i) = char(c - 32)
      end if
    end do
  end subroutine uppercase

end program affiche_topologie
