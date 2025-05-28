module mol2_types
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  type :: covalence_data
    character(len=2) :: element
    real(kind=dp)    :: rayon_simple, rayon_double, rayon_triple
  end type covalence_data
  type :: atom_info
    character(len=8) :: elem
    real(kind=dp)    :: x, y, z
  end type atom_info
  type :: bond
    integer           :: a1, a2
    character(len=10) :: type_liaison
  end type bond
end module mol2_types

module mol2_io
  use mol2_types
  implicit none
contains

  subroutine read_covalence(nom_fic, covs, nb_cov)
    character(len=*), intent(in)  :: nom_fic
    type(covalence_data), allocatable, intent(out) :: covs(:)
    integer, intent(out) :: nb_cov
    integer :: u, stat, i, ios
    character(len=256) :: ligne, tligne
    character(len=2)   :: es
    character(len=3)   :: rs
    real(kind=dp)       :: val

    open(newunit=u, file=nom_fic, status='old', action='read', iostat=stat)
    if (stat /= 0) then
      print *, "Erreur ouverture covalence:", trim(nom_fic)
      stop
    end if

    nb_cov = 0
    do
      read(u,'(A)', iostat=stat) ligne
      if (stat<0) exit
      if (stat>0) then
        print *, "Erreur lecture covalence"
        stop
      end if
      tligne = trim(ligne)
      if (len(tligne) < 5 .or. tligne(1:1) == "#") cycle
      nb_cov = nb_cov + 1
    end do
    rewind(u)

    allocate(covs(nb_cov))
    i = 0
    do
      read(u,'(A)', iostat=stat) ligne
      if (stat<0) exit
      tligne = trim(ligne)
      if (len(tligne) < 5 .or. tligne(1:1) == "#") cycle
      i = i + 1
      es = tligne(4:5)
      if (len(tligne) >= 8) then
        rs = tligne(6:8)
        read(rs,*, iostat=ios) val
        covs(i)%rayon_simple = merge(val/100.0_dp, 0.0_dp, ios==0)
      else
        covs(i)%rayon_simple = 0.0_dp
      end if
      if (len(tligne) >= 11) then
        rs = tligne(9:11)
        read(rs,*, iostat=ios) val
        covs(i)%rayon_double = merge(val/100.0_dp, 0.0_dp, ios==0)
      else
        covs(i)%rayon_double =0.0_dp
      end if
      if (len(tligne) >= 14) then
        rs = tligne(12:14)
        read(rs,*, iostat=ios) val
        covs(i)%rayon_triple = merge(val/100.0_dp, 0.0_dp, ios==0)
      else
        covs(i)%rayon_triple = 0.0_dp
      end if
      covs(i)%element = trim(es)
    end do

    close(u)
    nb_cov = min(nb_cov, i)
  end subroutine read_covalence

  subroutine read_mol2(nom_fic, atomes, nb_atomes)
    character(len=*), intent(in) :: nom_fic
    type(atom_info), allocatable, intent(out) :: atomes(:)
    integer, intent(out) :: nb_atomes
    integer :: u, stat, i, idum, ipt
    character(len=256) :: ligne
    character(len=8) :: nm, tp

    open(newunit=u, file=nom_fic, status='old', action='read', iostat=stat)
    if (stat /= 0) then
      print *, "Erreur ouverture MOL2:", trim(nom_fic)
      stop
    end if

    ! lire nombre d'atomes 
    do
      read(u,'(A)', iostat=stat) ligne
      if (stat<0) exit
      if (trim(ligne) == "@<TRIPOS>MOLECULE") then
        read(u,'(A)', iostat=stat) ligne
        read(u,*, iostat=stat) nb_atomes
        exit
      end if
    end do
    if (nb_atomes <= 0) then
      print *, "Nombre d'atomes invalide"
      stop
    end if

    ! attendre section ATOM 
    do
      read(u,'(A)', iostat=stat) ligne
      if (stat<0) then
        print *, "Section ATOM non trouvée"
        stop
      end if
      if (trim(ligne) == "@<TRIPOS>ATOM") exit
    end do

    allocate(atomes(nb_atomes))
    do i = 1, nb_atomes
      read(u,'(A)', iostat=stat) ligne
      if (stat /= 0) then
        print *, "Erreur lecture ATOM"
        stop
      end if
      read(ligne, *, iostat=stat) idum, nm, atomes(i)%x, atomes(i)%y, atomes(i)%z, tp
      if (stat /= 0) then
        print *, "Parsing ATOM"
        stop
      end if
      ipt = index(tp, '.')
      if (ipt > 0) tp = tp(:ipt-1)
      call uppercase(tp)
      atomes(i)%elem = trim(tp)
    end do

    close(u)
  end subroutine read_mol2

  subroutine write_mol2(nom_fic, atomes, nb_atomes, liaisons, nb_liaisons)
    character(len=*), intent(in) :: nom_fic
    type(atom_info), intent(in)  :: atomes(:)
    type(bond),      intent(in)  :: liaisons(:)
    integer, intent(in) :: nb_atomes, nb_liaisons
    integer :: u, i

    open(newunit=u, file=nom_fic, status='replace', action='write')
    write(u,'(A)') "@<TRIPOS>MOLECULE"
    write(u,'(A)') "Ligand_docked"
    write(u,'(I6,1X,I6)') nb_atomes, nb_liaisons  
    write(u,'(A)') "@<TRIPOS>ATOM"
    do i = 1, nb_atomes
      write(u,'(I6,1X,A4,3F10.4,1X,A2,1X,I3,1X,A4,1X,F8.4)') &
           i, trim(atomes(i)%elem), atomes(i)%x, atomes(i)%y, atomes(i)%z,  &
           trim(atomes(i)%elem), 1, "LIG", 0.0000_dp
    end do
    write(u,'(A)') "@<TRIPOS>BOND"
    do i = 1, nb_liaisons
      write(u,'(I6,1X,I4,1X,I4,1X,A1)') &
           i, liaisons(i)%a1, liaisons(i)%a2, trim(liaisons(i)%type_liaison)
    end do
    close(u)
  end subroutine write_mol2

  subroutine uppercase(str)
    character(len=*), intent(inout) :: str
    integer :: i, c
    do i = 1, len(str)
      c = ichar(str(i:i))
      if (c>=ichar('a') .and. c<=ichar('z')) then
        str(i:i) = char(c - 32)
      end if
    end do
  end subroutine uppercase

  pure function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=12) :: s
    write(s,'(I0)') i
  end function itoa

end module mol2_io

module geometry
  use mol2_types
  implicit none
contains

  subroutine calculate_distance(at, i, j, d)
    type(atom_info), intent(in) :: at(:)
    integer, intent(in) :: i, j
    real(kind=dp), intent(out) :: d
    real(kind=dp) :: dx, dy, dz

    dx = at(i)%x - at(j)%x
    dy = at(i)%y - at(j)%y
    dz = at(i)%z - at(j)%z
    d  = sqrt(dx*dx + dy*dy + dz*dz)
  end subroutine calculate_distance

  subroutine build_topology(at, nat, covs, nc, liaisons, nl)
    type(atom_info), intent(in)     :: at(:)
    type(covalence_data), intent(in):: covs(:)
    integer, intent(in)             :: nat, nc
    type(bond), intent(out)         :: liaisons(:)
    integer, intent(out)            :: nl
    integer :: i, j, k
    real(kind=dp) :: d, r1, r2, rcut
    logical :: f1, f2

    nl = 0
    do i = 1, nat-1
      do j = i+1, nat
        f1 = .false.; f2 = .false.
        do k = 1, nc
          if (trim(at(i)%elem) == trim(covs(k)%element)) then
            r1 = covs(k)%rayon_simple; f1 = .true.
          end if
          if (trim(at(j)%elem) == trim(covs(k)%element)) then
            r2 = covs(k)%rayon_simple; f2 = .true.
          end if
          if (f1 .and. f2) exit
        end do
        if (.not.(f1 .and. f2)) cycle
        rcut = (r1 + r2) * 1.35_dp
        call calculate_distance(at, i, j, d)
        if (d <= rcut) then
          nl = nl + 1
          liaisons(nl)%a1 = i
          liaisons(nl)%a2 = j
          liaisons(nl)%type_liaison = '1'
        end if
      end do
    end do
  end subroutine build_topology

  subroutine apply_transform(at, nat, tx, ty, tz, alpha, beta, gamma)
    type(atom_info), intent(inout) :: at(:)
    integer, intent(in)            :: nat
    real(kind=dp), intent(in)      :: tx, ty, tz, alpha, beta, gamma
    real(kind=dp) :: R(3,3), x, y, z, xp, yp, zp
    integer :: i

    ! Construction de la matrice de rotation
    R(1,1) = cos(alpha)*cos(beta)
    R(1,2) = cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)
    R(1,3) = cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)
    R(2,1) = sin(alpha)*cos(beta)
    R(2,2) = sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)
    R(2,3) = sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)
    R(3,1) = -sin(beta)
    R(3,2) = cos(beta)*sin(gamma)
    R(3,3) = cos(beta)*cos(gamma)

    do i = 1, nat
      x = at(i)%x; y = at(i)%y; z = at(i)%z
      xp = R(1,1)*x + R(1,2)*y + R(1,3)*z + tx
      yp = R(2,1)*x + R(2,2)*y + R(2,3)*z + ty
      zp = R(3,1)*x + R(3,2)*y + R(3,3)*z + tz
      at(i)%x = xp; at(i)%y = yp; at(i)%z = zp
    end do
  end subroutine apply_transform

end module geometry

module hbonds
  use mol2_types
  implicit none
contains

  integer function evaluate_hbonds(lig, nlig, site, nsite, lb, nlb, sb, nsb)
    type(atom_info), intent(in) :: lig(:), site(:)
    type(bond),      intent(in) :: lb(:), sb(:)
    integer, intent(in)         :: nlig, nsite, nlb, nsb
    integer :: i, j, k, neighbor
    real(kind=dp) :: d, angle, dx1, dy1, dz1, dx2, dy2, dz2
    real(kind=dp) :: norm1, norm2, dot

    evaluate_hbonds = 0
    do i = 1, nlig
      if (lig(i)%elem == 'H') then
        do j = 1, nsite
          if (site(j)%elem == 'O' .or. site(j)%elem == 'N') then
            dx1 = lig(i)%x - site(j)%x
            dy1 = lig(i)%y - site(j)%y
            dz1 = lig(i)%z - site(j)%z
            d   = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
            if (d >= 2.2_dp .and. d <= 4.0_dp) then
              neighbor = 0
              do k = 1, nsb
                if (sb(k)%a1 == j) then
                  neighbor = sb(k)%a2
                else if (sb(k)%a2 == j) then
                  neighbor = sb(k)%a1
                end if
                if (neighbor /= 0) exit
              end do
              if (neighbor == 0) cycle
              dx2 = site(neighbor)%x - site(j)%x
              dy2 = site(neighbor)%y - site(j)%y
              dz2 = site(neighbor)%z - site(j)%z
              norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
              norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)
              dot   = dx1*dx2 + dy1*dy2 + dz1*dz2
              angle = acos(dot/(norm1*norm2)) * 180.0_dp / acos(-1.0_dp)
              if (angle >= 90.0_dp .and. angle <= 150.0_dp) then
                evaluate_hbonds = evaluate_hbonds + 1
              end if
            end if
          end if
        end do
      end if
    end do
  end function evaluate_hbonds

end module hbonds

module ga_ops
  use mol2_types
  implicit none
contains

  subroutine init_population(pop, npop, ng)
    real(kind=dp), intent(out) :: pop(:,:)
    integer, intent(in)        :: npop, ng
    integer :: i, j
    real(kind=dp) :: r

    call random_seed()
    do i = 1, npop
      do j = 1, 3
        call random_number(r)
        pop(i,j) = r * 2.0_dp * acos(-1.0_dp)
      end do
      do j = 4, ng
        call random_number(r)
        pop(i,j) = (r - 0.5_dp) * 20.0_dp
      end do
    end do
  end subroutine init_population

  subroutine mutate(ind, ng, pm)
    real(kind=dp), intent(inout) :: ind(:)
    integer, intent(in)          :: ng
    real(kind=dp), intent(in)    :: pm
    integer :: j
    real(kind=dp) :: r

    do j = 1, ng
      call random_number(r)
      if (r < pm) then
        call random_number(r)
        if (j <= 3) then
          ind(j) = modulo(ind(j) + (r - 0.5_dp)*0.1745_dp, 2.0_dp*acos(-1.0_dp))
        else
          ind(j) = ind(j) + (r - 0.5_dp)*2.0_dp
        end if
      end if
    end do
  end subroutine mutate

  subroutine tournament_selection(pop, score, npop, sel)
    real(kind=dp), intent(in)  :: pop(:,:), score(:)
    real(kind=dp), intent(out) :: sel(:,:)
    integer, intent(in)        :: npop
    integer :: i, a, b
    real(kind=dp) :: r

    do i = 1, npop
      call random_number(r); a = 1 + int(r*npop)
      call random_number(r); b = 1 + int(r*npop)
      if (score(a) >= score(b)) then
        sel(i,:) = pop(a,:)
      else
        sel(i,:) = pop(b,:)
      end if
    end do
  end subroutine tournament_selection

  subroutine one_point_crossover(p1, p2, c1, c2, ng)
    real(kind=dp), intent(in)  :: p1(:), p2(:)
    real(kind=dp), intent(out) :: c1(:), c2(:)
    integer, intent(in)        :: ng
    integer :: pt
    real(kind=dp) :: r

    call random_number(r)
    pt = 1 + int(r*(ng-1))
    c1(1:pt)    = p1(1:pt)
    c2(1:pt)    = p2(1:pt)
    c1(pt+1:ng) = p2(pt+1:ng)
    c2(pt+1:ng) = p1(pt+1:ng)
  end subroutine one_point_crossover

end module ga_ops

program genetic_docking
  use omp_lib            ! Pour pouvoir ajuster OMP_NUM_THREADS si nécessaire
  use mol2_types
  use mol2_io
  use geometry
  use hbonds
  use ga_ops
  implicit none

  integer, parameter :: pop_size=100, max_gen=500, runs=10
  real(kind=dp), parameter :: p_crossover=0.95_dp, p_mutation=0.05_dp
  real(kind=dp) :: pop(pop_size,6), sel(pop_size,6), newpop(pop_size,6)
  real(kind=dp) :: scores(pop_size), stats(3,(max_gen/10)+1), r_cross
  character(len=256) :: lig_in, lig_out, site_fic, out_dir, run_dir, stats_file
  character(len=12) :: run_str, ind_str
  type(atom_info), allocatable :: lig_atoms(:), site_atoms(:)
  type(covalence_data), allocatable :: covs(:)
  type(bond), allocatable :: lig_bonds(:), site_bonds(:)
  integer :: n_lig, n_site, n_cov, n_lig_b, n_site_b
  integer :: run, gen, i, idx, ret

  call random_seed()

  ! Arguments
  call get_command_argument(1, lig_in)
  call get_command_argument(2, site_fic)
  call get_command_argument(3, out_dir)
  if (len_trim(lig_in)==0 .or. len_trim(site_fic)==0) then
    print *, "Usage: genetic_docking ligand.mol2 site.mol2 [output_dir]"
    stop
  end if
  if (.not.(len_trim(out_dir)>0)) out_dir = 'results2'

  call execute_command_line('mkdir -p '//trim(out_dir), exitstat=ret)
  if (ret /= 0) then
    print *, "Erreur création dossier:", trim(out_dir)
    stop
  end if

  ! Lecture données statiques
  call read_covalence("CoV_radii", covs, n_cov)
  call read_mol2(site_fic, site_atoms, n_site)
  allocate(site_bonds(n_site*(n_site-1)/2))
  call build_topology(site_atoms, n_site, covs, n_cov, site_bonds, n_site_b)

  do run = 1, runs
    write(run_str,'(I2.2)') run
    run_dir = trim(out_dir)//'/run'//trim(run_str)
    call execute_command_line('mkdir -p '//trim(run_dir), exitstat=ret)
    if (ret /= 0) then
      print *, "Erreur création dossier:", trim(run_dir)
      stop
    end if

    ! Initialisation population
    call init_population(pop, pop_size, 6)

    ! Boucle génétique
    do gen = 1, max_gen

      ! Évaluation parallèle de chaque individu
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i, lig_atoms, n_lig, lig_bonds, n_lig_b) &
      !$OMP&       SHARED(pop, scores, lig_in, site_atoms, n_site, covs, n_cov, site_bonds, n_site_b)
      do i = 1, pop_size
        call read_mol2(lig_in, lig_atoms, n_lig)
        allocate(lig_bonds(n_lig*(n_lig-1)/2))
        call build_topology(lig_atoms, n_lig, covs, n_cov, lig_bonds, n_lig_b)
        call apply_transform(lig_atoms, n_lig, pop(i,4), pop(i,5), pop(i,6), &
                             pop(i,1), pop(i,2), pop(i,3))
        scores(i) = evaluate_hbonds(lig_atoms, n_lig, site_atoms, n_site, &
                                    lig_bonds, n_lig_b, site_bonds, n_site_b)
        deallocate(lig_atoms, lig_bonds)
      end do
      !$OMP END PARALLEL DO

      ! Statistiques tous les 10 itérations
      if (mod(gen,10)==0 .or. gen==max_gen) then
        idx = gen/10 + 1
        stats(1,idx) = minval(scores)
        stats(2,idx) = maxval(scores)
        stats(3,idx) = sum(scores)/pop_size
      end if

      ! Sélection
      call tournament_selection(pop, scores, pop_size, sel)

      ! Crossover + mutation
      do i = 1, pop_size, 2
        call random_number(r_cross)
        if (r_cross < p_crossover) then
          call one_point_crossover(sel(i,:), sel(i+1,:), newpop(i,:), newpop(i+1,:), 6)
        else
          newpop(i,:)   = sel(i,:)
          newpop(i+1,:) = sel(i+1,:)
        end if
        call mutate(newpop(i,:),   6, p_mutation)
        call mutate(newpop(i+1,:), 6, p_mutation)
      end do

      pop = newpop

    end do  ! gen

    ! Écriture des conformères finaux (série)
    do i = 1, pop_size
      write(run_str,'(I2.2)') run
      write(ind_str,'(I3.3)') i
      lig_out = trim(run_dir)//'/lig_run'//trim(run_str)//'_ind'//trim(ind_str)//'.mol2'
      call read_mol2(lig_in, lig_atoms, n_lig)
      allocate(lig_bonds(n_lig*(n_lig-1)/2))
      call apply_transform(lig_atoms, n_lig, pop(i,4), pop(i,5), pop(i,6), &
                           pop(i,1), pop(i,2), pop(i,3))
      call build_topology(lig_atoms, n_lig, covs, n_cov, lig_bonds, n_lig_b)
      call write_mol2(lig_out, lig_atoms, n_lig, lig_bonds, n_lig_b)
      deallocate(lig_atoms, lig_bonds)
    end do

    ! CSV des statistiques
    stats_file = trim(out_dir)//'/stats_run'//trim(run_str)//'.csv'
    open(unit=20, file=trim(stats_file), status='replace', action='write')
    write(20,'(A)') 'gen,min,max,mean'
    do idx = 1, size(stats,2)
      write(20,'(I4,",",F8.2,",",F8.2,",",F8.2)') (idx-1)*10, &
           stats(1,idx), stats(2,idx), stats(3,idx)
    end do
    close(20)

  end do  ! run

end program genetic_docking
