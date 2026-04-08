! =============================================================================
! bpshell4.f90
!
! This Fortran code Computes the pairwise base-pair (BP) energy U_BP for every bead
! pair (i,j) in a single snapshot.
!
!
!   U_BP = U°_bp * exp(U_bp,bond + U_bp,angle + U_bp,dihedral)
!
! where U_bp,bond, U_bp,angle, and U_bp,dihedral encode the distance,
! four bond-angle, and two dihedral constraints required for Watson–Crick
! (WC) geometry in A-form RNA.
!
! A pair is accepted as a base pair if U_BP < -3 k_B T (~-1.77 kcal/mol),
! following the criterion in the Methods section of the paper.
!
! The code computes BPs separately for intra-chain and inter-chain pairs,
! then resolves one-to-one (mutually best-partner) pairings, and writes
! each accepted BP with its shell index relative to the droplet CoM.
!
! Input:  input.lammpstrj  -- LAMMPS custom dump with 6 columns per bead:
!           col 1: global bead index (continuous across chains)
!           col 2: bead type  (encodes nucleotide identity;
!                              G-C pair selected by |type_i - type_j| = 2)
!           col 3: molecule (chain) ID
!           col 4,5,6: x, y, z coordinates (Angstrom)
!
! Output: a1b1.dat  -- one line per accepted BP per frame:
!           frame_index, shell_index, bead_i, bead_j,
!           shell_of_chain_i, shell_of_chain_j, |U_BP(i,j)|
! =============================================================================
module site_graph_module41
    implicit none
contains
    ! Subroutine: generate_site_graph_fortran
    !
    ! For a single snapshot (ds1), computes U_BP(i,j) for all candidate
    ! bead pairs and stores the result in the bp(nd,nd) matrix.
    !
    ! Arguments:
    !   ds1   -- snapshot data (nd x 6): [bead_id, type, chain_id, x, y, z]
    !   a1    -- legacy (nrows)
    !   a2    -- legacy (ncols)
    !   a3    -- nd: number of beads in the droplet (= nlines)
    !   a4    -- number of data columns per bead (= 6)
    !   bp    -- output: bp(i,j) = |U_BP(i,j)| if pair qualifies, else 0
    ! -------------------------------------------------------------------------

        subroutine generate_site_graph_fortran(ds1,a1,a2,a3,a4,bp)
    implicit none
    integer :: a1, a2, a3, a4,f1,f2
    integer,parameter:: nd=a3
    double precision, intent(in):: ds1(a3,a4)
    real, intent(inout):: bp(nd,nd)
    double precision :: rij1(nd,nd)
        integer :: i, j
    double precision :: dx,dy,dz
    double precision :: ubp0,kr,kth,kph,ubpb,ubpa,ubpd,ub
    double precision ::X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4
    double precision ::th1num,th1den1,th1den2,th1,th2,th3,th4,t1c,t2c
    double precision ::NX1,NX2,NY1,NY2,NZ1,NZ2,ph1num,ph1den1,ph1den2,phi1,phi2
    double precision ::ph1c,ph2c,rbp0,arg
    real :: n1(3), n2(3), m1(3), b2(3), dot_m1_b2, sign_phi
    logical :: any_match



     rbp0=13.8   ! equilibrium base-pair separation (Å)
     kr=3.0      ! bond-distance force constant (Å^-2)
     kth=1.5     ! bond-angle force constant (rad^-2)
     kph=0.5     ! dihedral force constant
     t1c=1.8326  ! θ_1 = target angle for θ_{i,j,j-1} and θ_{i-1,i,j} (rad)
     t2c=0.9425  ! θ_2 = target angle for θ_{i,j,j+1} and θ_{i+1,i,j} (rad)
     ph1c=1.8326 ! φ_1 = offset for dihedral φ_{j-1,j,i,i-1} (rad)
     ph2c=1.1345 ! φ_2 = offset for dihedral φ_{j+1,j,i,i+1} (rad)
     ubp0=-5.0   ! U°_bp = -5.0 kcal/mol for G-C WC base pair


        ! --- Nucleotide-type selectivity: only G-C WC base pairs ---
        ! ds1(:,2) encodes the bead type. The SIS model only forms WC base
        ! pairs between G(bead type=3)  and C(bead type=1), not A (bead type=2). If the bead types of i and j differ
        ! by exactly 2, they are a G-C pair. Other combinations (C-C, A-G,
        ! A-A, etc.) are skipped.

   do i = 2, nd-1
    do j = i+1, nd-1
!        if (ds1(i,2) .ne. 2 .and. ds1(j,2) .ne. 2) then
        if (abs(ds1(i,2)-ds1(j,2)) .eq. 2) then
                 ! --- Boundary protection: exclude beads at chain position 1 ---
           if (modulo(int(ds1(i,1)),93) .gt. 1 .and. modulo(int(ds1(j,1)),93) .gt. 1 )then
           !intrachain bp
           if (ds1(i,3) .eq. ds1(j,3))then

              ! Require sequence separation >= 5 along the chain to avoid
              ! trivially close beads         
           if (abs(ds1(i,1)-ds1(j,1)) .ge. 5)then
            dx = ds1(i,4) - ds1(j,4)
            dy = ds1(i,5) - ds1(j,5)
            dz = ds1(i,6) - ds1(j,6)

            rij1(i,j) = sqrt(dx**2 + dy**2 + dz**2)
                ! Distance pre-filter: accept only pairs in [10, 18] Å,
                ! bracketing the equilibrium r_bp,0 = 13.8 Å.
            if (rij1(i,j) .ge. 10 .and. rij1(i,j) .le. 18)then
            ! Bond-distance term: penalises deviation from r_bp,0.
            ubpb=-kr*(rij1(i,j)-rbp0)**2

            ! th1: angle i-j-( j-1 )
            X1=ds1(i,4)
            X2=ds1(j,4)
            X3=ds1(j-1,4)
            Y1=ds1(i,5)
            Y2=ds1(j,5)
            Y3=ds1(j-1,5)
            Z1=ds1(i,6)
            Z2=ds1(j,6)
            Z3=ds1(j-1,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            th1=acos(th1num/(th1den1*th1den2))
            ! th2: angle ( i-1 )-i-j
            X1=ds1(i-1,4)
            X2=ds1(i,4)
            X3=ds1(j,4)
            Y1=ds1(i-1,5)
            Y2=ds1(i,5)
            Y3=ds1(j,5)
            Z1=ds1(i-1,6)
            Z2=ds1(i,6)
            Z3=ds1(j,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            th2=acos(th1num/(th1den1*th1den2))
            ! th3: angle i-j-( j+1 )
            X1=ds1(i,4)
            X2=ds1(j,4)
            X3=ds1(j+1,4)
            Y1=ds1(i,5)
            Y2=ds1(j,5)
            Y3=ds1(j+1,5)
            Z1=ds1(i,6)
            Z2=ds1(j,6)
            Z3=ds1(j+1,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            th3=acos(th1num/(th1den1*th1den2))
            ! th4: angle ( i+1 )-i-j
            X1=ds1(i+1,4)
            X2=ds1(i,4)
            X3=ds1(j,4)
            Y1=ds1(i+1,5)
            Y2=ds1(i,5)
            Y3=ds1(j,5)
            Z1=ds1(i+1,6)
            Z2=ds1(i,6)
            Z3=ds1(j,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            th4=acos(th1num/(th1den1*th1den2))

                  ! U_bp,angle = -k_θ * [(th1-θ1)^2 + (th2-θ1)^2 +
                  !                       (th3-θ2)^2 + (th4-θ2)^2]

            ubpa=-kth*((th1-t1c)**2+(th2-t1c)**2+(th3-t2c)**2+(th4-t2c)**2)

                              ! --- U_bp,dihedral---
                  ! Two dihedral angles:
                  !   phi1: dihedral along the sequence j-1, j, i, i-1
                  !   phi2: dihedral along the sequence j+1, j, i, i+1
                  !
                  ! phi1: plane normals for (j-1,j,i) and (i-1,i,j) planes

            X1=ds1(j-1,4)
            X2=ds1(j,4)
            X3=ds1(i,4)
            X4=ds1(i-1,4)
            Y1=ds1(j-1,5)
            Y2=ds1(j,5)
            Y3=ds1(i,5)
            Y4=ds1(i-1,5)
            Z1=ds1(j-1,6)
            Z2=ds1(j,6)
            Z3=ds1(i,6)
            Z4=ds1(i-1,6)
            ! Normal to plane (j-1,j,i): N1 = (j-1 - j) × (i - j)
            NX1=(Y1-Y2)*(Z3-Z2)-(Z1-Z2)*(Y3-Y2)
            NY1=(Z1-Z2)*(X3-X2)-(X1-X2)*(Z3-Z2)
            NZ1=(X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2)
            ! Normal to plane (i-1,i,j): N2 = (i-1 - i) × (j - i)
            NX2=(Y4-Y3)*(Z3-Z2)-(Z4-Z3)*(Y3-Y2)
            NY2=(Z4-Z3)*(X3-X2)-(X4-X3)*(Z3-Z2)
            NZ2=(X4-X3)*(Y3-Y2)-(Y4-Y3)*(X3-X2)

!!#######################################################
! Signed dihedral via cross product of normals dotted with
                  ! the bond vector b2 = i - j (j to i direction).

n1 = [NX1, NY1, NZ1]
n2 = [NX2, NY2, NZ2]

! Compute cross product of normals
!m1 = cross_product(n1, n2)
m1(1) = n1(2)*n2(3) - n1(3)*n2(2)
m1(2) = n1(3)*n2(1) - n1(1)*n2(3)
m1(3) = n1(1)*n2(2) - n1(2)*n2(1)


! Compute bond vector (from j to i)
b2 = [X3 - X2, Y3 - Y2, Z3 - Z2]

! Determine sign
dot_m1_b2 = dot_product(m1, b2)
sign_phi = sign(1.0d0, dot_m1_b2)

! Compute signed dihedral angle
ph1num = dot_product(n1, n2)
ph1den1 = sqrt(sum(n1**2))
ph1den2 = sqrt(sum(n2**2))
phi1 = sign_phi * acos(ph1num / (ph1den1 * ph1den2))
!!#######################################################
 ! phi2: dihedral along j+1, j, i, i+1

            X1=ds1(j+1,4)
            X2=ds1(j,4)
            X3=ds1(i,4)
            X4=ds1(i+1,4)
            Y1=ds1(j+1,5)
            Y2=ds1(j,5)
            Y3=ds1(i,5)
            Y4=ds1(i+1,5)
            Z1=ds1(j+1,6)
            Z2=ds1(j,6)
            Z3=ds1(i,6)
            Z4=ds1(i+1,6)

            NX1=(Y1-Y2)*(Z3-Z2)-(Z1-Z2)*(Y3-Y2)
            NY1=(Z1-Z2)*(X3-X2)-(X1-X2)*(Z3-Z2)
            NZ1=(X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2)

            NX2=(Y4-Y3)*(Z3-Z2)-(Z4-Z3)*(Y3-Y2)
            NY2=(Z4-Z3)*(X3-X2)-(X4-X3)*(Z3-Z2)
            NZ2=(X4-X3)*(Y3-Y2)-(Y4-Y3)*(X3-X2)
 

!!#######################################################
n1 = [NX1, NY1, NZ1]
n2 = [NX2, NY2, NZ2]

!m1 = cross_product(n1, n2)
m1(1) = n1(2)*n2(3) - n1(3)*n2(2)
m1(2) = n1(3)*n2(1) - n1(1)*n2(3)
m1(3) = n1(1)*n2(2) - n1(2)*n2(1)


! Compute bond vector (from j to i)
b2 = [X3 - X2, Y3 - Y2, Z3 - Z2]

! Determine sign
dot_m1_b2 = dot_product(m1, b2)
sign_phi = sign(1.0d0, dot_m1_b2)

! Compute signed dihedral angle
ph1num = dot_product(n1, n2)
ph1den1 = sqrt(sum(n1**2))
ph1den2 = sqrt(sum(n2**2))
phi2 = sign_phi * acos(ph1num / (ph1den1 * ph1den2))
!!#######################################################
            ! U_bp,dihedral = -k_φ * [(1+cos(phi1+φ1)) + (1+cos(phi2+φ2))]

            ubpd=-kph*((1+cos(phi1+ph1c))+(1+cos(phi2+ph2c)))

             ! U_BP = U°_bp * exp(U_bp,bond + U_bp,angle + U_bp,dihedral)
            ub=ubp0*exp(ubpb+ubpa+ubpd)

                  ! Store |U_BP| if it exceeds the pre-filter threshold (~k_BT).
                  ! The main program applies the full -3 k_B T criterion later.

            if (abs(ub) .ge. 0.59) bp(i,j)=abs(ub)


            

            end if
            end if    
            else
            ! =================================================================
            ! INTER-CHAIN BRANCH: beads i and j on different chains
            ! =================================================================

            dx = ds1(i,4) - ds1(j,4)
            dy = ds1(i,5) - ds1(j,5)
            dz = ds1(i,6) - ds1(j,6)

            rij1(i,j) = sqrt(dx**2 + dy**2 + dz**2)
            if (rij1(i,j) .ge. 10 .and. rij1(i,j) .le. 18)then
            
            ubpb=-kr*(rij1(i,j)-rbp0)**2
            X1=ds1(i,4)
            X2=ds1(j,4)
            X3=ds1(j-1,4)
            Y1=ds1(i,5)
            Y2=ds1(j,5)
            Y3=ds1(j-1,5)
            Z1=ds1(i,6)
            Z2=ds1(j,6)
            Z3=ds1(j-1,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5
             
            arg = th1num/(th1den1*th1den2)
            if (arg > 1.0d0) arg = 1.0d0
            if (arg < -1.0d0) arg = -1.0d0
            th1 = acos(arg)
!            th1=acos(th1num/(th1den1*th1den2))
            
            X1=ds1(i-1,4)
            X2=ds1(i,4)
            X3=ds1(j,4)
            Y1=ds1(i-1,5)
            Y2=ds1(i,5)
            Y3=ds1(j,5)
            Z1=ds1(i-1,6)
            Z2=ds1(i,6)
            Z3=ds1(j,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

           arg = th1num/(th1den1*th1den2)
           if (arg > 1.0d0) arg = 1.0d0
           if (arg < -1.0d0) arg = -1.0d0
           th2 = acos(arg)
!            th2=acos(th1num/(th1den1*th1den2))

            X1=ds1(i,4)
            X2=ds1(j,4)
            X3=ds1(j+1,4)
            Y1=ds1(i,5)
            Y2=ds1(j,5)
            Y3=ds1(j+1,5)
            Z1=ds1(i,6)
            Z2=ds1(j,6)
            Z3=ds1(j+1,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            arg = th1num/(th1den1*th1den2)
            if (arg > 1.0d0) arg = 1.0d0
            if (arg < -1.0d0) arg = -1.0d0
            th3 = acos(arg)
!            th3=acos(th1num/(th1den1*th1den2))

            X1=ds1(i+1,4)
            X2=ds1(i,4)
            X3=ds1(j,4)
            Y1=ds1(i+1,5)
            Y2=ds1(i,5)
            Y3=ds1(j,5)
            Z1=ds1(i+1,6)
            Z2=ds1(i,6)
            Z3=ds1(j,6)
            th1num=(X1-X2)*(X3-X2)+(Y1-Y2)*(Y3-Y2)+(Z1-Z2)*(Z3-Z2)
            th1den1=((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)**0.5
            th1den2=((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)**0.5

            arg = th1num/(th1den1*th1den2)
            if (arg > 1.0d0) arg = 1.0d0
            if (arg < -1.0d0) arg = -1.0d0
            th4 = acos(arg)
!            th4=acos(th1num/(th1den1*th1den2))


            ubpa=-kth*((th1-t1c)**2+(th2-t1c)**2+(th3-t2c)**2+(th4-t2c)**2)

!             print*,'ubpa',real(ubpa),real(th1),real(th2),real(th3),real(th4)
!            j-1,j,i,i-1
            X1=ds1(j-1,4)
            X2=ds1(j,4)
            X3=ds1(i,4)
            X4=ds1(i-1,4)
            Y1=ds1(j-1,5)
            Y2=ds1(j,5)
            Y3=ds1(i,5)
            Y4=ds1(i-1,5)
            Z1=ds1(j-1,6)
            Z2=ds1(j,6)
            Z3=ds1(i,6)
            Z4=ds1(i-1,6)

            NX1=(Y1-Y2)*(Z3-Z2)-(Z1-Z2)*(Y3-Y2)
            NY1=(Z1-Z2)*(X3-X2)-(X1-X2)*(Z3-Z2)
            NZ1=(X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2)

            NX2=(Y4-Y3)*(Z3-Z2)-(Z4-Z3)*(Y3-Y2)
            NY2=(Z4-Z3)*(X3-X2)-(X4-X3)*(Z3-Z2)
            NZ2=(X4-X3)*(Y3-Y2)-(Y4-Y3)*(X3-X2)

            ph1num=(NX1-0)*(NX2-0)+(NY1-0)*(NY2-0)+(NZ1-0)*(NZ2-0)
            ph1den1=((NX1-0)**2+(NY1-0)**2+(NZ1-0)**2)**0.5
            ph1den2=((NX2-0)**2+(NY2-0)**2+(NZ2-0)**2)**0.5

            arg = ph1num/(ph1den1*ph1den2)
            if (arg > 1.0d0) arg = 1.0d0
            if (arg < -1.0d0) arg = -1.0d0
            phi1 = acos(arg)
!            phi1=acos(ph1num/(ph1den1*ph1den2))
!            j+1,j,i,i+1

            X1=ds1(j+1,4)
            X2=ds1(j,4)
            X3=ds1(i,4)
            X4=ds1(i+1,4)
            Y1=ds1(j+1,5)
            Y2=ds1(j,5)
            Y3=ds1(i,5)
            Y4=ds1(i+1,5)
            Z1=ds1(j+1,6)
            Z2=ds1(j,6)
            Z3=ds1(i,6)
            Z4=ds1(i+1,6)

            NX1=(Y1-Y2)*(Z3-Z2)-(Z1-Z2)*(Y3-Y2)
            NY1=(Z1-Z2)*(X3-X2)-(X1-X2)*(Z3-Z2)
            NZ1=(X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2)

            NX2=(Y4-Y3)*(Z3-Z2)-(Z4-Z3)*(Y3-Y2)
            NY2=(Z4-Z3)*(X3-X2)-(X4-X3)*(Z3-Z2)
            NZ2=(X4-X3)*(Y3-Y2)-(Y4-Y3)*(X3-X2)
            
            ph1num=(NX1-0)*(NX2-0)+(NY1-0)*(NY2-0)+(NZ1-0)*(NZ2-0)
            ph1den1=((NX1-0)**2+(NY1-0)**2+(NZ1-0)**2)**0.5
            ph1den2=((NX2-0)**2+(NY2-0)**2+(NZ2-0)**2)**0.5

            arg = ph1num/(ph1den1*ph1den2)
            if (arg > 1.0d0) arg = 1.0d0
            if (arg < -1.0d0) arg = -1.0d0
            phi2 = acos(arg)
!            phi2=acos(ph1num/(ph1den1*ph1den2))

            ubpd=-kph*((1+cos(phi1+ph1c))+(1+cos(phi2+ph2c)))
            ub=ubp0*exp(ubpb+ubpa+ubpd)

            if (abs(ub) .ge. 0.59) bp(i,j)=abs(ub)

            end if
            end if
            end if
        end if
    end do
!    print*,i,'************************************************************************************************'
end do
        

    end subroutine generate_site_graph_fortran
end module site_graph_module41
program traj

    use site_graph_module41
    implicit none
  ! --- System dimensions ---
    ! nrows:   total number of frames in the LAMMPS dump (set per simulation).
    ! ncols:   number of chains in the droplet for this run (must be set > 0).
    ! nc = 6:  columns per bead in the dump (bead_id, type, mol, x, y, z).
    ! nlines:  beads per frame = ncols * 93 (93 beads per CAG31 chain).
    integer, parameter :: nrows = 6590, ncols = 64, nc=6 !ncols number of chains
    integer, parameter :: chn_len=93
    integer, parameter :: nlines = ncols * 93
    integer :: i, j, m1, n, iostat_i,j1,j2
    integer :: chn(nrows, ncols),e1,e2
    integer :: chind(ncols),st1,fn1
    double precision :: data1(nrows*nlines, nc),r
    double precision :: dx,dy,dz,d1,d2,dij
    double precision :: dxch,dych,dzch,chd
    double precision :: dsone(nlines, nc)
    integer :: n1(ncols), n2(ncols),intra,inter,ds
    
    integer:: n11(nrows/2),ck(nrows),ind,nn ! Choose maximum possible c_n3 and nd
    integer, parameter :: max_c_n3 = ncols
    integer, parameter :: max_nd = max_c_n3 * 93   ! Ensure this largeenough for all calls
        ! bp(i,j): stores |U_BP(i,j)| for pairs passing the pre-filter |ub| >= 0.59.
    ! Later filtered to the paper criterion |ub| > 3*0.59 ≈ 1.77 ≈ 3 k_B T.
    real :: bp(max_nd,max_nd),L=642.7
    integer :: best_partner(max_nd)     ! Track best partner for each base
    real    :: best_energy(max_nd)      ! Track highest energy for each base
    logical :: paired(max_nd)
    character(len=20) :: file_name
    integer :: i0,ji,chone,chtwo
    integer :: chain_idx, bead_idx, start_bead
    double precision :: shift_x, shift_y, shift_z
    double precision :: dxc(ncols), dyc(ncols), dzc(ncols)


    ! --- Read entire LAMMPS dump into memory ---
    ! The file has no header lines between frames; every row is a bead record.
    open(unit=3, file='input.lammpstrj', status='old')
    do i = 1, nrows*nlines
        read(3, *) (data1(i, j), j = 1, nc)
    end do
    close(3)

    open(unit=4, file='a1b1.dat', status='unknown')       
 
!=========================================================
!=========================================================


!         bp=0

        
        do n=1,nrows
        bp=0.0
          nn=n
                  ! --- Compute droplet centre of mass (CoM) ---
        ! Simple arithmetic mean of all bead x/y/z in this frame.

            dsone = data1(nlines*(nn-1)+1 : nlines*nn, 1:nc)


dx=sum(dsone(:,4))
dy=sum(dsone(:,5))
dz=sum(dsone(:,6))

dx=dx/max_nd
dy=dy/max_nd
dz=dz/max_nd
        ! --- Assign each chain to a concentric shell around the CoM ---
        ! Shell index = floor(distance_from_CoM / 50 Å).
        ! This matches the shell geometry used in shellwise_overlap.f90
        ! (shell thickness = 50 Å).

do i=1,ncols
st1=(i-1)*chn_len+1
fn1=i*chn_len
dxch=sum(dsone(st1:fn1,4))
dych=sum(dsone(st1:fn1,5))
dzch=sum(dsone(st1:fn1,6))
dxch=dxch/chn_len
dych=dych/chn_len
dzch=dzch/chn_len
chd=((dxch-dx)**2+(dych-dy)**2+(dzch-dz)**2)**0.5
if (int(chd/50.0) .lt. 10)then
chind(i)=int(chd/50.0)
else
chind(i)=10
end if
end do


            call generate_site_graph_fortran(dsone, nrows,ncols, nlines, nc, bp)

                    ! ======================================================================
        ! One-to-one base-pair resolution
        !
        ! The SIS model enforces that each nucleotide forms at most one WC base
        ! pair. Here, for each bead i we record its highest-|U_BP| candidate
        ! partner. We then accept only MUTUAL best-partner pairs (i's best is j
        ! AND j's best is i), marking both as paired to avoid double-counting.
        ! This is consistent with the base-pair criterion in the paper Methods.
        !
        ! Pre-filter: only pairs with |U_BP| > 3 * 0.59 ≈ 1.77 kcal/mol ≈ 3 k_BT
        ! are considered (paper Methods: "If U_BP < -3 k_B T, a base pair is formed").
        ! ======================================================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_partner = 0
best_energy = 0.0
paired = .false.


        ! Accept pair (i,j) only if i's best is j AND j's best is i (mutual),
        ! and neither has been paired already.

do i = 1, max_nd
  do j = i+1, max_nd
    if (bp(i,j) > 3.0*0.59) then

     !===================uncomment when checkin for all base-pairs done=========
      if (bp(i,j) > best_energy(i)) then
        best_energy(i) = bp(i,j)
        best_partner(i) = j
      end if
      if (bp(i,j) > best_energy(j)) then
        best_energy(j) = bp(i,j)
        best_partner(j) = i
      end if
!================================================================



    end if

  end do
end do

 do i = 1, max_nd
  j = best_partner(i)
          ! Accept pair (i,j) only if i's best is j AND j's best is i (mutual),
        ! and neither has been paired already.
  if (j > i .and. best_partner(j) == i .and. .not. paired(i) .and. .not. paired(j)) then
    paired(i) = .true.
    paired(j) = .true.
                ! Shell index of this pair: based on the average distance of
            ! the two bead positions from the droplet CoM.
    d1=((dsone(i,4)-dx)**2+(dsone(i,5)-dy)**2+(dsone(i,6)-dz)**2)**0.5
    d2=((dsone(j,4)-dx)**2+(dsone(j,5)-dy)**2+(dsone(j,6)-dz)**2)**0.5
    dij=(d1+d2)/2
    ds=int(dij/50)
    chone=int((i-1)/chn_len)+1
    chtwo=int((j-1)/chn_len)+1
                ! Output columns:
            ! frame, shell_of_pair, bead_i, bead_j,
            !   shell_of_chain_i, shell_of_chain_j, |U_BP(i,j)|
     write(4, 10) n, ds, i, j, chind(chone),chind(chtwo),bp(i,j)
    end if
end do
!=====================================================================



      end do

10 format(I5, 2x, I5, 2x, I5, 2x, I5, 2x, I5, 2x, I5,2x, F20.7)

close(4)


end program traj


