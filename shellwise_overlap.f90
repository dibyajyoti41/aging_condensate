! =============================================================================
! shellwise_overlap.f90
!
! Purpose: Computes the pairwise overlap function F(tau) and the pairwise
!          mean-squared change in inter-bead distance <(dr_ij)^2> restricted
!          to beads located within a concentric spherical shell of the droplet.
!
! The shell is defined relative to the droplet centre of mass (CoM):
!   inner radius = cutoff
!   outer radius = cutoff + 50 (= cutoff2)
!
! Both F(tau) and <(dr_ij)^2> are computed for three disjoint shells
! (cutoff = 0, 50, 100) in a single run, and further decomposed into
! intra-chain and inter-chain contributions (as in droplet_overlap.f90).
!
! Key difference from droplet_overlap.f90:
!   A pair (i,j) is included in the average only if BOTH beads i and j lie
!   within the shell at BOTH times t and t+tau. This ensures that the
!   measured dynamics truly reflect beads that reside in the shell throughout
!   the observation window, not beads that pass through transiently.
!
! Input:  input.xyz  -- droplet-only trajectory in XYZ format
! Output: a1b1.dat   -- shellwise results (see write statement for column order)
! =============================================================================

module site_graph_module5
    implicit none
contains
            ! -------------------------------------------------------------------------
    ! Subroutine: generate_site_graph_fortran
    !
    ! Arguments identical to droplet_overlap.f90 with one addition:
    !   cutoff -- inner radius of the shell (distance from droplet CoM).
    !             The outer radius is implicitly cutoff + 50.
    ! -------------------------------------------------------------------------
        subroutine generate_site_graph_fortran(ds1,ds2,a1,a2,a3,a4,eq32,eq62,eq3i,eq6i,eq3o,eq6o,max_nd,cutoff)
    implicit none
    integer :: a1, a2, a3, a4,f1,f2
    integer,parameter:: nd=a3
    double precision,parameter:: L=642.7
    double precision, intent(in):: ds1(a3,a4), ds2(a3,a4)
    double precision, intent(out) :: eq32, eq62
    double precision, intent(out) :: eq3i, eq6i, eq3o, eq6o
    double precision :: rij1, rij2
    double precision :: dx,dy,dz
    integer, intent(in) :: max_nd
    integer :: i, j, c_n1, c_n2, c_n3, c_n4, count_in
    integer :: k, k_max, i1, denom,denomi,denomo,chn_len
    double precision :: sum_overlap1, sum_overlap2, sum_xy1, sum_xy2, diff_val, rc
    double precision :: sintra, sinter, sum_xyi, sum_xyo
    double precision :: xyval1, xyval2
    double precision :: xyvali, xyvalo,r1,r2
    double precision :: dx1,dy1,dz1,dx2,dy2,dz2
    integer,intent(in) :: cutoff
    integer :: cutoff2
    logical:: within_rc1(nd), within_rc2(nd)
    logical :: any_match
!   character(len=1) :: atom_symbol  ! To read the "C"
!    count_in = 11

    
    rc = 6.0d0
    cutoff2=cutoff+50 ! shell thickness is fixed at 50 distance units
    chn_len=93  

!


    denom = 0
    denomi = 0
    denomo = 0
    sum_xy1 = 0.0d0
    sum_xy2 = 0.0d0
    sum_xyi = 0.0d0
    sum_xyo = 0.0d0
    sum_overlap1 = 0.0d0
    sum_overlap2 = 0.0d0
    sintra = 0.0d0
    sinter = 0.0d0

    ! --- Step 1: Compute droplet centre of mass (CoM) for each snapshot ---
    ! The CoM is used as the reference origin for shell membership.
    ! Note: this is a simple arithmetic mean of all bead positions and does
    ! NOT account for PBC wrapping of the droplet across box boundaries.

dx1=0
dy1=0
dz1=0
dx2=0
dy2=0
dz2=0

do i = 1, nd
            dx1 = dx1+ds1(i,1)
            dy1 = dy1+ds1(i,2)
            dz1 = dz1+ds1(i,3)

            dx2 = dx2+ds2(i,1)
            dy2 = dy2+ds2(i,2)
            dz2 = dz2+ds2(i,3)
end do
! Normalise to get CoM coordinates at time t (dx1,dy1,dz1) and t+tau (dx2,dy2,dz2)
dx1=dx1*1.0/nd
dy1=dy1*1.0/nd
dz1=dz1*1.0/nd
dx2=dx2*1.0/nd
dy2=dy2*1.0/nd
dz2=dz2*1.0/nd
    ! --- Step 2: Determine shell membership for every bead at both times ---
    ! Bead i is "within the shell" at time t if its distance from the CoM
    ! satisfies:  cutoff <= r1 < cutoff2
do i = 1, nd
        r1 = sqrt((ds1(i,1)-dx1)**2 + (ds1(i,2)-dy1)**2 + (ds1(i,3)-dz1)**2)
        r2 = sqrt((ds2(i,1)-dx2)**2 + (ds2(i,2)-dy2)**2 + (ds2(i,3)-dz2)**2)
        within_rc1(i) = (r1 .ge. cutoff .and. r1 <= cutoff2)  ! shell membership at t
        within_rc2(i) = (r2 .ge. cutoff .and. r2 <= cutoff2)  ! shell membership at t+tau

end do
    ! --- Step 3: Pairwise overlap restricted to shell-resident pairs ---
    ! A pair (i,j) is counted only if both beads are inside the shell at
    ! both times t and t+tau (persistent shell membership).

do i = 1, nd
    do j = i+2, nd ! j=i+1 excluded to skip directly bonded neighbours
           if ((within_rc1(i) .and. within_rc1(j))) then   ! both in shell at t
           if ((within_rc2(i) .and. within_rc2(j))) then   ! both in shell at t+tau

            ! Distance at time t1 (ds1)
            dx = ds1(i,1) - ds1(j,1)
            dy = ds1(i,2) - ds1(j,2)
            dz = ds1(i,3) - ds1(j,3)

            dx = dx - L * nint(dx*1.0 / L)
            dy = dy - L * nint(dy*1.0 / L)
            dz = dz - L * nint(dz*1.0 / L)

            rij1 = sqrt(dx**2 + dy**2 + dz**2)


             
            ! Distance at time t2 (ds2)
            dx = ds2(i,1) - ds2(j,1)
            dy = ds2(i,2) - ds2(j,2)
            dz = ds2(i,3) - ds2(j,3)

            dx = dx - L * nint(dx*1.0 / L)
            dy = dy - L * nint(dy*1.0 / L)
            dz = dz - L * nint(dz*1.0 / L)

            rij2 = sqrt(dx**2 + dy**2 + dz**2)


            denom = denom + 1
            xyval2=dabs(rij2-rij1)
            sum_xy2=sum_xy2+xyval2**2
            if (rc - xyval2 >= 0.0d0) sum_overlap2 = sum_overlap2+1.0d0

            if (((i-1)/chn_len) .eq. ((j-1)/chn_len))then
            denomi = denomi + 1
             xyvali=dabs(rij2-rij1)
             sum_xyi=sum_xyi+xyvali**2
             if (rc - xyvali >= 0.0d0) sintra = sintra+1.0d0
            else
            denomo = denomo + 1
            xyvalo=dabs(rij2-rij1)
            sum_xyo=sum_xyo+xyvalo**2
            if (rc - xyvalo >= 0.0d0) sinter = sinter+1.0d0
            end if



            end if
            end if
    end do
end do


    if (denom > 0) then
!
        eq32 = sum_xy2/denom
        eq62 = sum_overlap2/denom
    else
!
        eq32 = 0.0d0
        eq62 = 0.0d0

    end if

        if (denomi > 0)then
        eq3i=sum_xyi/denomi
        eq6i=sintra/denomi
        else
        eq3i=0.0
        eq6i=0.0
        end if

        if (denomo > 0)then
        eq3o=sum_xyo/denomo
        eq6o=sinter/denomo
        else
        eq3o=0.0
        eq6o=0.0

        end if


    end subroutine generate_site_graph_fortran
end module site_graph_module5

program traj
    !this code calculates shellwise overlap function and pair-wise mean squared displacement. 
    ! this code also takes only droplet co-ordinates (.xyz format) as input.
    use site_graph_module5
    implicit none
    integer, parameter :: nrows = 83629, ncols=00, nc=3
    integer, parameter :: nlines = ncols * 93
    integer :: i, j, m1, n, iostat_i,jj
    integer :: chn(nrows, ncols),ck1(nrows)
    double precision :: data1(nrows*nlines, nc),r
    double precision :: dsone(nlines, nc), dstwo(nlines, nc)
    double precision :: eq31, eq61, res_val, res_m1, oc, oz,res_eq3,res_eq
    double precision :: oci,ozi,oco,ozo
    double precision :: eq32, eq62,res_val2, res_m12, res_eq32,res_eq2
    double precision :: eq3i,eq6i,res_vali, res_mi, res_eqi
    double precision :: eq3o,eq6o,res_valo, res_mo, res_eqo
    integer :: n1(ncols), n2(ncols),co
    integer:: n11(nrows/2),ck(nrows),ind,j1,nn ! Choose maximum possible c_n3 and nd
    integer, parameter :: max_c_n3 = ncols
    integer, parameter :: max_nd = max_c_n3 * 93   ! Ensure this largeenough for all calls
    character(len=1) :: atom_symbol  ! To read the "C"
    integer :: values(8)
    integer :: coff(5)

    do i=1,3
    coff(i)=50*(i-1) !cutoff 
    end do




    open(unit=3, file='input.xyz', status='old')
    do i = 1, nrows
    ! 1. Read the two header lines per frame
    read(3, *) ! Skips the atom count (5580)
    read(3, *) ! Skips the comment line (Atoms. Timestep: ...)
    do j = 1, nlines
        ind = (i-1)*nlines + j
        read(3, *) atom_symbol, data1(ind, 1), data1(ind, 2), data1(ind, 3)
    end do
    
end do
close(3)


    open(unit=4, file='a1b1.dat', status='replace',action='write',iostat=iostat_i)
    if (iostat_i /= 0) then
        print *, "Error opening output.dat for writing"
        stop
    end if

  ! --- Outer loop over shells ---
          do jj=1,3   
        co=coff(jj) ! inner radius of the current shell

        

        ! --- Loop over lag times tau = m1 ---
        do m1=1,nrows/2
        res_m1 = 0.0
        res_m12 = 0.0
        res_eq3 = 0.0
        res_eq32 = 0.0
        oc = 0.0
        oz = 0.0

        res_mi = 0.0
        res_eqi = 0.0
        ozi = 0.0
        res_mo = 0.0
        res_eqo = 0.0
        ozo = 0.0

        j1=1
        ck=0        
        n11(1:nrows/2)=0

            ! --- Random sampling of reference times (without replacement) ---
            ! For each lag time m1, randomly pick 10000 unique reference frames
            ! from [1, nrows-m1] to use as the reference time t.
            ! Limiting to 10000 samples (rather than nrows/2) reduces compute
            ! cost while providing sufficient statistical averaging. 

!        do while (j1<(nrows/2)-1)
        do while (j1<10000)! 
        call random_number(r)
        ind=int(r*(nrows-m1))+1
        if (ck(ind) .eq. 0 )then
        n11(j1)=ind
        ck(ind)=1
        j1=j1+1
        end if
        end do
        
        ! --- Time-averaging loop over sampled reference times ---
        do n=1,10000
            nn=n11(n)
            dsone = data1(nlines*(nn-1)+1 : nlines*nn, 1:nc)   !trajecory at t
            dstwo = data1(nlines*(nn+m1-1)+1 : nlines*(nn+m1), 1:nc)! trajectory at t+\tau


            call generate_site_graph_fortran(dsone, dstwo, nrows,ncols, nlines, nc, &
                 eq32,eq62,eq3i,eq6i,eq3o,eq6o,max_nd,co)


            res_m12 = res_m12 + eq62
            res_eq32 = res_eq32 + eq32

            res_mi = res_mi + eq6i
            res_eqi = res_eqi + eq3i

            res_mo = res_mo + eq6o
            res_eqo = res_eqo + eq3o




             oc = oc + 1
        end do

        if (oc .gt. 0)then
        res_val2 = res_m12/(oc)
        res_eq2=res_eq32/(oc)

        res_vali = res_mi/(oc)
        res_eqi  = res_eqi/(oc)

        res_valo = res_mo/(oc)
        res_eqo  = res_eqo/(oc)
        else
        print*,'valid trajectory not found'
        end if

! Output columns:
            ! shell_inner_radius, tau(snapshot), F, <|dr_ij|^2>,
            !   F_intra, <|dr_ij|^2>_intra, F_inter, <|dr_ij|^2>_inter
            ! Results are written shell-by-shell (shell1 block, then shell2, shell3).
       write(4,10)co,m1, res_val2,res_eq2, res_vali, res_eqi, res_valo, res_eqo
10     FORMAT(I5,2x,I5,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4)
    end do


    end do

    close(4)


end program traj

