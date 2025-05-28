! complete_mol2.f90
program complete_mol2
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! Déclarations des types utilisés
  type :: covalence_data
    character(len=2) :: element
    real(kind=dp)    :: rayon_simple
    real(kind=dp)    :: rayon_double
    real(kind=dp)    :: rayon_triple
  end type covalence_data

  type :: atom_info
    character(len=8) :: elem
    real(kind=dp)    :: x, y, z
  end type atom_info

  type :: bond
    integer          :: a1, a2
    character(len=10):: bond_type
  end type bond

  ! Variables pour les fichiers et les données
  character(len=256) :: fichier_mol2, fichier_covalence, fichier_sortie
  type(covalence_data), allocatable :: covalences(:)
  type(atom_info), allocatable      :: atomes(:)
  type(bond), allocatable           :: liaisons(:)
  integer :: nb_atomes, nb_covalence
  integer :: nombre_liaisons
  logical :: tous_connectes

  real :: temps_debut, temps_fin
  ! Variables pour la gestion de la topologie
  integer :: pourcentage_delta, idx, jdx, kdx
  real(kind=dp) :: delta, distance_ij
  character(len=2) :: element_i, element_j

  ! Variables pour la lecture du fichier MOL2 complet
  character(len=256), allocatable :: lignes_fichier(:)
  integer :: nb_lignes, idx_ligne, ligne_section_atome, premiere_apres_atome

  ! Variables de lecture de ligne
  character(len=256) :: ligne

  
  ! Récupération des paramètres
  call cpu_time(temps_debut)
  call get_command_argument(1, fichier_mol2)
  call get_command_argument(2, fichier_covalence)

  if (len_trim(fichier_mol2) == 0 .or. len_trim(fichier_covalence) == 0) then
    print *, "Usage: complete_mol2 <fichier.mol2> <fichier_rayons.txt>"
    stop
  end if

  ! Calcul du nom du fichier de sortie (ex : 1QSN_WITH_BONDS.mol2)
  call compute_output_filename(trim(fichier_mol2), fichier_sortie)

  
  ! Lecture du fichier de covalence
  call read_covalence(fichier_covalence, covalences, nb_covalence)

  ! Affichage  de la liste des éléments
  print *, "Eléments du fichier de covalence lus :"
  do idx = 1, nb_covalence
    write(*,'(A3,2X,F6.3,2X,F6.3,2X,F6.3)') trim(covalences(idx)%element), &
         covalences(idx)%rayon_simple, covalences(idx)%rayon_double, covalences(idx)%rayon_triple
  end do
  print *, "-----------"

  
  !Lecture du fichier MOL2 pour récupérer la liste des atomes
  call read_mol2(fichier_mol2, atomes, nb_atomes)

  ! Affichage des atomes
  print *, "Nombre d'atomes lus :", nb_atomes

  ! Vérification que chaque atome a un rayon défini
  call check_elements(atomes, nb_atomes, covalences, nb_covalence)

  
  ! Calcul des liaisons (topologie)
  ! On essaie différents delta de tolérance
  do pourcentage_delta = 10, 35, 5
    delta = real(pourcentage_delta, kind=dp) / 100.0_dp

    ! Allouer la table des liaisons (taille max : n*(n-1)/2)
    if (allocated(liaisons)) deallocate(liaisons)
    allocate(liaisons(nb_atomes*(nb_atomes-1)/2))
    nombre_liaisons = 0

    do idx = 1, nb_atomes - 1
      do jdx = idx + 1, nb_atomes
        call calculate_distance(atomes, idx, jdx, distance_ij)
        element_i = trim(atomes(idx)%elem)
        element_j = trim(atomes(jdx)%elem)

        ! Teste triple, double, puis simple
        if (try_bond_type('triple', element_i, element_j, distance_ij, delta, idx, jdx, &
                           covalences, nb_covalence, nombre_liaisons, liaisons)) cycle
        if (try_bond_type('double', element_i, element_j, distance_ij, delta, idx, jdx, &
                           covalences, nb_covalence, nombre_liaisons, liaisons)) cycle
        if (try_bond_type('single', element_i, element_j, distance_ij, delta, idx, jdx, &
                           covalences, nb_covalence, nombre_liaisons, liaisons)) cycle
      end do
    end do

    call check_connectivity(nb_atomes, liaisons, nombre_liaisons, tous_connectes)
    if (tous_connectes) then
      print *, "Topologie trouvée avec un delta de", pourcentage_delta, "%"
      exit
    end if
    deallocate(liaisons)
  end do

  if (.not. tous_connectes) then
    print *, "Impossible de déterminer une topologie connectée avec delta<=35%."
    stop
  end if

  
  ! Lecture complète du fichier MOL2 dans une table de chaînes
  call count_file_lines(fichier_mol2, nb_lignes)
  allocate(lignes_fichier(nb_lignes))
  call read_file_lines(fichier_mol2, lignes_fichier)

  ! Recherche de la ligne contenant "@<TRIPOS>ATOM"
  ligne_section_atome = 0
  do idx_ligne = 1, nb_lignes
    if (trim(lignes_fichier(idx_ligne)) == "@<TRIPOS>ATOM") then
      ligne_section_atome = idx_ligne
      exit
    end if
  end do
  if (ligne_section_atome == 0) then
    print *, "Section ATOM non trouvée dans le fichier MOL2."
    stop
  end if

  ! La section ATOM contient exactement nb_atomes lignes, donc :
  premiere_apres_atome = ligne_section_atome + nb_atomes

  
  ! Écriture du nouveau fichier MOL2 complet
  open(unit=10, file=fichier_sortie, status='replace', action='write')
  ! Recopie des lignes depuis le début jusqu'à la fin de la section ATOM
  do idx_ligne = 1, premiere_apres_atome
    write(10, '(A)') trim(lignes_fichier(idx_ligne))
  end do

  ! Insertion de la section BOND
  write(10, *)
  write(10, '(A)') "@<TRIPOS>BOND"
  do kdx = 1, nombre_liaisons
    ! Format standard pour MOL2 : bond_id, atome origine, atome destination, type de liaison
    write(10, '(I4,2X,I4,2X,I4,2X,A10)') kdx, liaisons(kdx)%a1, liaisons(kdx)%a2, trim(liaisons(kdx)%bond_type)
  end do

  ! Recopie du reste du fichier 
  if (premiere_apres_atome < nb_lignes) then
    do idx_ligne = premiere_apres_atome+1, nb_lignes
      write(10, '(A)') trim(lignes_fichier(idx_ligne))
    end do
  end if
  close(10)
  print *, "Fichier complet créé :", trim(fichier_sortie)

  call cpu_time(temps_fin)
  print *, "Temps d'exécution :", temps_fin - temps_debut, "secondes"

contains

  
  subroutine compute_output_filename(nom_entree, nom_sortie)
    character(len=*), intent(in)  :: nom_entree
    character(len=256), intent(out) :: nom_sortie
    integer :: pos

    ! Si le nom contient "_NO_BOND", on le remplace par "_WITH_BONDS"
    pos = index(nom_entree, '_NO_BOND')
    if (pos > 0) then
      nom_sortie = nom_entree(1:pos-1) // "_WITH_BONDS" // nom_entree(pos+7:)
    else
      ! Sinon, on insère _WITH_BONDS avant l'extension
      pos = index(nom_entree, '.mol2')
      if (pos > 0) then
        nom_sortie = nom_entree(1:pos-1) // "_WITH_BONDS.mol2"
      else
        nom_sortie = nom_entree // "_WITH_BONDS"
      end if
    end if
  end subroutine compute_output_filename

  
  subroutine count_file_lines(nom_fichier, nb_lignes)
    character(len=*), intent(in) :: nom_fichier
    integer, intent(out) :: nb_lignes
    integer :: unite, statut_IO
    character(len=256) :: ligne_temp

    nb_lignes = 0
    open(newunit=unite, file=nom_fichier, status='old', action='read', iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture fichier :", trim(nom_fichier)
      stop
    end if
    do
      read(unite, '(A)', iostat=statut_IO) ligne_temp
      if (statut_IO /= 0) exit
      nb_lignes = nb_lignes + 1
    end do
    close(unite)
  end subroutine count_file_lines

  
  subroutine read_file_lines(nom_fichier, lignes)
    character(len=*), intent(in) :: nom_fichier
    character(len=256), intent(out) :: lignes(:)
    integer :: unite, idx, statut_IO

    open(newunit=unite, file=nom_fichier, status='old', action='read', iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture fichier :", trim(nom_fichier)
      stop
    end if
    idx = 0
    do
      idx = idx + 1
      read(unite, '(A)', iostat=statut_IO) lignes(idx)
      if (statut_IO /= 0) exit
    end do
    close(unite)
  end subroutine read_file_lines

  
  subroutine read_covalence(nom_fichier, covalences, nb_covalence)
    character(len=*), intent(in) :: nom_fichier
    type(covalence_data), allocatable, intent(out) :: covalences(:)
    integer, intent(out) :: nb_covalence

    integer :: unite, statut_IO, compteur, statut_temp
    character(len=256) :: ligne, ligne_trimmee
    character(len=2) :: element_str
    character(len=3) :: rayon_str
    real(kind=dp) :: val_tmp

    open(newunit=unite, file=nom_fichier, status='old', action='read', iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture fichier de covalence :", trim(nom_fichier)
      stop
    end if

    nb_covalence = 0
    do
      read(unite, '(A)', iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      if (statut_IO > 0) then
        print *, "Erreur lecture dans", trim(nom_fichier)
        stop
      end if
      ligne_trimmee = trim(ligne)
      if (len_trim(ligne_trimmee) == 0) cycle
      if (ligne_trimmee(1:1) == "#") cycle
      if (len(ligne_trimmee) < 5) cycle
      nb_covalence = nb_covalence + 1
    end do

    rewind(unite)
    allocate(covalences(nb_covalence))
    compteur = 0
    do
      read(unite, '(A)', iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      ligne_trimmee = trim(ligne)
      if (len_trim(ligne_trimmee) == 0) cycle
      if (ligne_trimmee(1:1) == "#") cycle
      if (len(ligne_trimmee) < 5) cycle
      compteur = compteur + 1
      element_str = ligne_trimmee(4:5)
      if (len(ligne_trimmee) >= 8) then
        rayon_str = ligne_trimmee(6:8)
        read(rayon_str, *, iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covalences(compteur)%rayon_simple = val_tmp / 100.0_dp
        else
          covalences(compteur)%rayon_simple = 0.0_dp
        end if
      else
        covalences(compteur)%rayon_simple = 0.0_dp
      end if
      if (len(ligne_trimmee) >= 11) then
        rayon_str = ligne_trimmee(9:11)
        read(rayon_str, *, iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covalences(compteur)%rayon_double = val_tmp / 100.0_dp
        else
          covalences(compteur)%rayon_double = 0.0_dp
        end if
      else
        covalences(compteur)%rayon_double = 0.0_dp
      end if
      if (len(ligne_trimmee) >= 14) then
        rayon_str = ligne_trimmee(12:14)
        read(rayon_str, *, iostat=statut_temp) val_tmp
        if (statut_temp == 0) then
          covalences(compteur)%rayon_triple = val_tmp / 100.0_dp
        else
          covalences(compteur)%rayon_triple = 0.0_dp
        end if
      else
        covalences(compteur)%rayon_triple = 0.0_dp
      end if
      covalences(compteur)%element = trim(element_str)
    end do
    close(unite)
    if (compteur < nb_covalence) nb_covalence = compteur
  end subroutine read_covalence

  
  subroutine read_mol2(nom_fichier, atomes, nb_atomes)
    character(len=*), intent(in) :: nom_fichier
    type(atom_info), allocatable, intent(out) :: atomes(:)
    integer, intent(out) :: nb_atomes

    integer :: unite, statut_IO, idx, id_dummy, pos_point
    character(len=256) :: ligne
    character(len=8)   :: nom_tmp, type_tmp
    real(kind=dp) :: x_tmp, y_tmp, z_tmp

    open(newunit=unite, file=nom_fichier, status='old', action='read', iostat=statut_IO)
    if (statut_IO /= 0) then
      print *, "Erreur ouverture MOL2 file :", trim(nom_fichier)
      stop
    end if

    nb_atomes = 0
    do
      read(unite, '(A)', iostat=statut_IO) ligne
      if (statut_IO < 0) exit
      if (trim(ligne) == "@<TRIPOS>MOLECULE") then
        read(unite, '(A)', iostat=statut_IO) ligne
        read(unite, *) nb_atomes
        exit
      end if
    end do

    if (nb_atomes <= 0) then
      print *, "Erreur: nombre d'atomes non lu dans la section MOLECULE."
      stop
    end if

    do
      read(unite, '(A)', iostat=statut_IO) ligne
      if (statut_IO < 0) then
        print *, "Section @<TRIPOS>ATOM non trouvée."
        stop
      end if
      if (trim(ligne) == "@<TRIPOS>ATOM") exit
    end do

    allocate(atomes(nb_atomes))
    do idx = 1, nb_atomes
      read(unite, '(A)', iostat=statut_IO) ligne
      if (statut_IO /= 0) then
        print *, "Erreur lecture atome :", ligne
        stop
      end if
      read(ligne, *, iostat=statut_IO) id_dummy, nom_tmp, x_tmp, y_tmp, z_tmp, type_tmp
      if (statut_IO /= 0) then
        print *, "Erreur de parsing pour l'atome :", ligne
        stop
      end if
      pos_point = index(type_tmp, '.')
      if (pos_point > 0) then
        type_tmp = type_tmp(:pos_point-1)
      end if
      call uppercase(type_tmp)
      atomes(idx)%elem = trim(type_tmp)
      atomes(idx)%x = x_tmp
      atomes(idx)%y = y_tmp
      atomes(idx)%z = z_tmp
    end do
    close(unite)
  end subroutine read_mol2

  
  subroutine check_elements(atomes, n, covalences, nb_covalence)
    type(atom_info), intent(in) :: atomes(:)
    integer, intent(in) :: n
    type(covalence_data), intent(in) :: covalences(:)
    integer, intent(in) :: nb_covalence
    integer :: idx, jdx
    logical :: trouve

    do idx = 1, n
      trouve = .false.
      do jdx = 1, nb_covalence
        if (trim(atomes(idx)%elem) == trim(covalences(jdx)%element)) then
          trouve = .true.
          exit
        end if
      end do
      if (.not. trouve) then
        print *, "Element inconnu dans le fichier de covalence :", trim(atomes(idx)%elem)
        stop
      end if
    end do
  end subroutine check_elements

  
  subroutine calculate_distance(atomes, idx, jdx, distance_ij)
    type(atom_info), intent(in) :: atomes(:)
    integer, intent(in) :: idx, jdx
    real(kind=dp), intent(out) :: distance_ij
    real(kind=dp) :: dx, dy, dz

    dx = atomes(idx)%x - atomes(jdx)%x
    dy = atomes(idx)%y - atomes(jdx)%y
    dz = atomes(idx)%z - atomes(jdx)%z
    distance_ij = sqrt(dx*dx + dy*dy + dz*dz)
  end subroutine calculate_distance

  
  logical function try_bond_type(b_type, e1, e2, distance_ij, delta, idx, jdx, covalences, nb_covalence, nombre_liaisons, liaisons)
  character(len=*), intent(in) :: b_type, e1, e2
  real(kind=dp), intent(in) :: distance_ij, delta
  integer, intent(in) :: idx, jdx, nb_covalence
  type(covalence_data), intent(in) :: covalences(:)
  integer, intent(inout) :: nombre_liaisons
  type(bond), intent(inout) :: liaisons(:)
  real(kind=dp) :: r1, r2, somme_r, lower, upper

  r1 = get_radius(e1, b_type, covalences, nb_covalence)
  r2 = get_radius(e2, b_type, covalences, nb_covalence)
  somme_r = r1 + r2

  if (somme_r <= 0.0_dp) then
    try_bond_type = .false.
    return
  end if

  lower = somme_r * (1.0_dp - delta)
  upper = somme_r * (1.0_dp + delta)
  if (distance_ij >= lower .and. distance_ij <= upper) then
    nombre_liaisons = nombre_liaisons + 1
    liaisons(nombre_liaisons)%a1 = idx
    liaisons(nombre_liaisons)%a2 = jdx
    ! Conversion du type de liaison en chiffre :
    select case(trim(b_type))
      case('single')
        liaisons(nombre_liaisons)%bond_type = '1'
      case('double')
        liaisons(nombre_liaisons)%bond_type = '2'
      case('triple')
        liaisons(nombre_liaisons)%bond_type = '3'
      case default
        liaisons(nombre_liaisons)%bond_type = b_type
    end select
    try_bond_type = .true.
  else
    try_bond_type = .false.
  end if
end function try_bond_type

  
  real(kind=dp) function get_radius(e, btype, covalences, nb_covalence)
    character(len=*), intent(in) :: e, btype
    type(covalence_data), intent(in) :: covalences(:)
    integer, intent(in) :: nb_covalence
    integer :: k

    get_radius = 0.0_dp
    do k = 1, nb_covalence
      if (trim(covalences(k)%element) == trim(e)) then
        select case(btype)
          case('single')
            get_radius = covalences(k)%rayon_simple
          case('double')
            get_radius = covalences(k)%rayon_double
          case('triple')
            get_radius = covalences(k)%rayon_triple
        end select
        return
      end if
    end do
  end function get_radius

  
  subroutine check_connectivity(n, liaisons, nombre_liaisons, tous_connectes)
    integer, intent(in) :: n, nombre_liaisons
    type(bond), intent(in) :: liaisons(:)
    logical, intent(out) :: tous_connectes
    integer, allocatable :: connectes(:)
    integer :: idx

    allocate(connectes(n))
    connectes = 0
    do idx = 1, nombre_liaisons
      connectes(liaisons(idx)%a1) = connectes(liaisons(idx)%a1) + 1
      connectes(liaisons(idx)%a2) = connectes(liaisons(idx)%a2) + 1
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

end program complete_mol2
