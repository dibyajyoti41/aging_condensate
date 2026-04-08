! =============================================================================
! droplet_overlap.f90
!
! Purpose: Computes the pairwise overlap function F(tau) and the pairwise
!          mean-squared change in inter-bead distance <(dr_ij)^2> for an
!          RNA condensate droplet trajectory (CAG31, 93 beads/chain).
!
! The overlap function for a bead pair (i,j) at lag time tau is defined as:
!   F_ij(tau) = 1  if |r_ij(t+tau) - r_ij(t)| < rc  
!               0  otherwise
! F(tau) is then averaged over all unique non-bonded pairs and over all
! reference times t (time-translation average).
!
! Both F(tau) and <(dr_ij)^2>(tau) are computed separately for:
!   (a) all pairs within the droplet
!   (b) intra-chain pairs only  (i,j on the same chain)
!   (c) inter-chain pairs only  (i,j on different chains)
!
! Input:  input.xyz  -- droplet-only trajectory in XYZ format
!                       (nrows frames, nlines beads per frame)
! Output: a1b1.dat   -- columns: tau, F, <dr^2>, F_intra, <dr^2>_intra,
!                                       F_inter, <dr^2>_inter
! =============================================================================
module site_graph_module5
    implicit none
contains
    ! -------------------------------------------------------------------------
    ! Subroutine: generate_site_graph_fortran
    !
    ! Computes the pairwise overlap function and mean-squared pairwise distance
    ! change between two snapshots ds1 (at time t) and ds2 (at time t+tau).
    !
    ! Arguments:
    !   ds1, ds2  -- bead coordinate arrays for the two snapshots (nd x nc)
    !   a1, a2    -- unused legacy dimension arguments (nrows, ncols)
    !   a3        -- number of beads in the droplet (nd = nlines)
    !   a4        -- number of coordinate dimensions (nc = 3)
    !   eq32      -- output: <|r_ij(t+tau) - r_ij(t)|^2> over all pairs
    !   eq62      -- output: F(tau) -- overlap function over all pairs
    !   eq3i/eq6i -- same quantities restricted to intra-chain pairs
    !   eq3o/eq6o -- same quantities restricted to inter-chain pairs
    ! -------------------------------------------------------------------------        
        subroutine generate_site_graph_fortran(ds1,ds2,a1,a2,a3,a4,eq32,eq62,eq3i,eq6i,eq3o,eq6o,max_nd)
    implicit none
    integer :: a1, a2, a3, a4,f1,f2
    integer,parameter:: nd=a3    ! total number of beads in this snapshot
    double precision,parameter:: L=642.7  ! simulation box side length (LJ units)
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
    double precision :: xyvali, xyvalo
    logical :: any_match

 ! rc: overlap criterion -- a pair (i,j) is counted as "remembered" at lag
    ! tau if the change in their separation |r_ij(t+tau) - r_ij(t)| < rc.   
    ! chn_len: number of beads per polymer chain (CAG31 coarse-grained model).
    rc = 6.0d0  
    chn_len=93

!

 ! Initialise counters and accumulators
    denom = 0  ! total number of valid pairs counted
    denomi = 0 !intra-chain counter
    denomo = 0 !inter-chain counter
    sum_xy1 = 0.0d0
    sum_xy2 = 0.0d0 ! accumulates |dr_ij|^2 over all pairs
    sum_xyi = 0.0d0 ! accumulates |dr_ij|^2 over intra-chain pairs
    sum_xyo = 0.0d0 ! accumulates |dr_ij|^2 over inter-chain pairs
    sum_overlap1 = 0.0d0
    sum_overlap2 = 0.0d0 ! counts pairs satisfying the overlap criterion (all)
    sintra = 0.0d0
    sinter = 0.0d0



    ! Double loop over all unique non-bonded bead pairs (i, j).
    ! j starts at i+2 to exclude the directly bonded neighbour j=i+1

do i = 1, nd
    do j = i+2, nd
            ! --- Pairwise distance at reference time t (snapshot ds1) ---
            dx = ds1(i,1) - ds1(j,1)
            dy = ds1(i,2) - ds1(j,2)
            dz = ds1(i,3) - ds1(j,3)
            ! Apply minimum-image periodic boundary conditions (PBC)
            dx = dx - L * nint(dx*1.0 / L)
            dy = dy - L * nint(dy*1.0 / L)
            dz = dz - L * nint(dz*1.0 / L)

            rij1 = sqrt(dx**2 + dy**2 + dz**2)

            ! --- Pairwise distance at lag time t+tau (snapshot ds2) ---
            dx = ds2(i,1) - ds2(j,1)
            dy = ds2(i,2) - ds2(j,2)
            dz = ds2(i,3) - ds2(j,3)
            ! Apply minimum-image PBC
            dx = dx - L * nint(dx*1.0 / L)
            dy = dy - L * nint(dy*1.0 / L)
            dz = dz - L * nint(dz*1.0 / L)

            rij2 = sqrt(dx**2 + dy**2 + dz**2)

              ! xyval2 = |r_ij(t+tau) - r_ij(t)| : scalar change in pair distance
            denom = denom + 1
            xyval2=dabs(rij2-rij1)
            sum_xy2=sum_xy2+xyval2**2
            ! Overlap condition: pair is valid/'remembered' if distance change < rc
            if (rc - xyval2 >= 0.0d0) sum_overlap2 = sum_overlap2+1.0d0
            ! Classify pair as intra- or inter-chain and accumulate separately
            if (((i-1)/chn_len) .eq. ((j-1)/chn_len))then
             ! --- Intra-chain pair ---        
            denomi = denomi + 1
             xyvali=dabs(rij2-rij1)
             sum_xyi=sum_xyi+xyvali**2
             if (rc - xyvali >= 0.0d0) sintra = sintra+1.0d0 
            else
              ! --- Inter-chain pair ---       
            denomo = denomo + 1
            xyvalo=dabs(rij2-rij1)
            sum_xyo=sum_xyo+xyvalo**2
            if (rc - xyvalo >= 0.0d0) sinter = sinter+1.0d0
            end if




    end do
end do

 ! --- Normalise and return results ---
    ! eq62 = F(tau): fraction of pairs satisfying the overlap criterion
    ! eq32 = <|dr_ij|^2>(tau): mean-squared change in pairwise distance

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
! =============================================================================
! Main program: reads the full droplet trajectory into memory, then computes
! the pairwise overlap function F(tau) for each lag time tau = m1 * dt by
! averaging over randomly sampled, non-repeating reference times t.
!
! System parameters (CAG31 condensate):
!   nrows    -- total number of trajectory frames stored in input.xyz
!   ncols    -- number of chains in the droplet (set at compile time)
!   nlines   -- number of beads per frame = ncols * 93
!   nc = 3   -- spatial dimensions (x, y, z)
! =============================================================================

        use site_graph_module5
    implicit none
    ! --- System dimensions ---
    ! nrows: total frames in the trajectory file.
    ! ncols: number of chains currently in the droplet (set per simulation).
    ! nlines: beads per frame = ncols chains * 93 beads/chain.

    integer, parameter :: nrows = 83629, ncols=00, nc=3
    integer, parameter :: nlines = ncols * 93
    integer :: i, j, m1, n, iostat_i
    integer :: chn(nrows, ncols),ck1(nrows)
    double precision :: data1(nrows*nlines, nc),r
    double precision :: dsone(nlines, nc), dstwo(nlines, nc)
    double precision :: eq31, eq61, res_val, res_m1, oc, oz,res_eq3,res_eq
    double precision :: oci,ozi,oco,ozo
    double precision :: eq32, eq62,res_val2, res_m12, res_eq32,res_eq2
    double precision :: eq3i,eq6i,res_vali, res_mi, res_eqi
    double precision :: eq3o,eq6o,res_valo, res_mo, res_eqo
    integer :: n1(ncols), n2(ncols)
    integer:: n11(nrows/2),ck(nrows),ind,j1,nn ! Choose maximum possible c_n3 and nd
    integer, parameter :: max_c_n3 = ncols
    integer, parameter :: max_nd = max_c_n3 * 93   ! Ensure this largeenough for all calls
        character(len=1) :: atom_symbol  ! To read the "C"
        integer :: values(8)


    ! --- Read entire trajectory into memory ---
    ! Each frame in input.xyz has two header lines followed by nlines bead lines.
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
   ! --- Outer loop over lag times tau = m1 (in units of snapshot interval) ---
    ! m1 runs from 1 to nrows/2 so that both t and t+tau are valid frame indices.
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

        ! --- Random sampling of reference times t (without replacement) ---
        ! n11 holds a set of nrows/2 unique, randomly chosen reference frame
        ! indices in [1, nrows-m1] so that the pair (t, t+tau) is always valid.
        ! ck(ind) acts as a "used" flag to prevent duplicate reference times.

        j1=1
        ck=0        
        n11(1:nrows/2)=0

        do while (j1<nrows/2) ! random frame in [1, nrows-m1]
        call random_number(r) ! accept only if not yet used
        ind=int(r*(nrows-m1))+1
        if (ck(ind) .eq. 0 )then
        n11(j1)=ind
        ck(ind)=1
        j1=j1+1
        end if
        end do
        
        
        ! --- Time-averaging loop over sampled reference times t ---
        do n=1,nrows/2 !this is time averaging loop, 
            nn=n11(n)
            ! Extract bead coordinates for frame nn (time t) and frame nn+m1 (time t+tau)
            dsone = data1(nlines*(nn-1)+1 : nlines*nn, 1:nc)
            dstwo = data1(nlines*(nn+m1-1)+1 : nlines*(nn+m1), 1:nc)


            call generate_site_graph_fortran(dsone, dstwo, nrows,ncols, nlines, nc, &
                 eq32,eq62,eq3i,eq6i,eq3o,eq6o,max_nd)

        ! Accumulate contributions from this reference time
            res_m12 = res_m12 + eq62
            res_eq32 = res_eq32 + eq32

            res_mi = res_mi + eq6i
            res_eqi = res_eqi + eq3i

            res_mo = res_mo + eq6o
            res_eqo = res_eqo + eq3o




             oc = oc + 1
        end do
! --- Normalise by number of reference times ---
        if (oc .gt. 0)then
        res_val2 = res_m12/(oc)
        res_eq2=res_eq32/(oc)

        res_vali = res_mi/(oc)
        res_eqi  = res_eqi/(oc)

        res_valo = res_mo/(oc)
        res_eqo  = res_eqo/(oc)
        else
        print*,'no valid trajectory found'
        end if



        ! Output columns:
        ! tau(snapshot), F(tau), <|dr_ij|^2>(tau),
        !   F_intra, <|dr_ij|^2>_intra,
        !   F_inter, <|dr_ij|^2>_inter                  
        write(4,10) m1, res_val2,res_eq2, res_vali, res_eqi, res_valo, res_eqo 
10     FORMAT(I5,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4,2x,f20.4)
    end do

    close(4)


end program traj
       
