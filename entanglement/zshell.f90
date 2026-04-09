program calculation
    implicit none
    
    ! --- Variables ---
    integer :: i, j, k, io_status
    integer :: frame_idx
    
    ! --- Histogram Arrays ---
    ! 0:6 corresponds to the shell bins
    double precision :: ent_hist(0:6)   ! Numerator: Total entanglements in this shell
    double precision :: chain_hist(0:6) ! Denominator: Total chains in this shell
    integer :: shell_idx
    
    ! --- Chain Processing Variables ---
    integer :: chain_ent_count          ! Entanglements in current chain
    logical :: reading_chain            ! State flag: Are we inside a chain?
    integer,parameter :: z_nbeads=93                 ! Number of beads in Z1+ chain
    double precision :: z_beads(2000,3) ! Buffer for Z1+ bead coords (assumed max 2000)
    double precision :: chain_cx, chain_cy, chain_cz ! Chain COM
    double precision :: z_box_dummy(3)  ! To skip Z1+ header
    
    ! --- File Reading Variables ---
    double precision :: rx, ry, rz, flag, col5
    double precision :: dist_com
    double precision :: xlo, xhi, ylo, yhi, zlo, zhi
    double precision :: bb(6591*3, 2)   ! Box bounds array
    integer :: b_idx

    ! --- Constants ---
    integer, parameter :: nframe = 6591
    integer, parameter :: rad = 50
    double precision, parameter :: cm = 321.5d0
    double precision, parameter :: box = 643.0d0
    double precision, parameter :: tol = 0.01d0 ! Tolerance for coordinate matching

    ! ==========================================
    ! 1. Read Box Bounds (Existing Logic)
    ! ==========================================
    open(unit=2, file='ab.dat', status='old')
    do i = 1, nframe * 3
        read(2, *, iostat=io_status) bb(i, 1), bb(i, 2)
        if (io_status /= 0) then
            print *, "Error reading boxb6.dat at line", i
            stop
        end if
    end do
    close(2)

    ! ==========================================
    ! 2. Initialize Processing
    ! ==========================================
    ent_hist = 0.0d0
    chain_hist = 0.0d0
    
    open(unit=3, file='cd.dat', status='old')          ! Entanglement data
    open(unit=5, file='gh.dat', status='old') ! Chain Geometry
    
    frame_idx = 0
    reading_chain = .false.
    chain_ent_count = 0

    ! ==========================================
    ! 3. Main Loop (Driven by spc6.dat)
    ! ==========================================
    do
        ! Read line from spc6.dat
        read(3, *, iostat=io_status) rx, ry, rz, flag, col5
        if (io_status /= 0) exit ! End of file

        ! --- CASE A: NEW FRAME HEADER (col5 == 5) ---
        if (int(col5) == 5) then
            frame_idx = frame_idx + 1
            
            ! Sync Z1+: Read the Frame Header (Box dimensions)
            ! This lines up Z1+ with the new frame
!            read(5, *) z_box_dummy(1), z_box_dummy(2), z_box_dummy(3)
            
            if (mod(frame_idx, 100) == 0) print *, "Processing Frame:", frame_idx

        ! --- CASE B: CHAIN BOUNDARY (col5 == 0) ---
        elseif (int(col5) == 0) then
            
            if (.not. reading_chain) then
                ! >>> START OF CHAIN <<<
                reading_chain = .true.
                chain_ent_count = 0
                
                ! 1. Read Chain Info from Z1+
 !               read(5, *) z_nbeads
                
 !               ! 2. Read All Beads for this chain from Z1+
                chain_cx = 0.0d0
                chain_cy = 0.0d0
                chain_cz = 0.0d0
                
                do k = 1, z_nbeads
                    read(5, *) z_beads(k, 1), z_beads(k, 2), z_beads(k, 3)
                    chain_cx = chain_cx + z_beads(k, 1)
                    chain_cy = chain_cy + z_beads(k, 2)
                    chain_cz = chain_cz + z_beads(k, 3)
                end do
                
                ! 3. Calculate Chain COM (Average)
                chain_cx = chain_cx / dble(z_nbeads)
                chain_cy = chain_cy / dble(z_nbeads)
                chain_cz = chain_cz / dble(z_nbeads)

                ! 4. SAFETY CHECK: Match spc6 start bead with Z1+ 1st bead
                if (abs(rx - z_beads(1,1)) > tol) then
                   print *, "MISMATCH Frame", frame_idx, " Start Bead"
                   print *, "spc6:", rx, ry, rz
                   print *, "Z1+ :", z_beads(1,1), z_beads(1,2), z_beads(1,3)
                   stop
                endif

            else
                ! >>> END OF CHAIN <<<
                reading_chain = .false.
                
                ! 1. SAFETY CHECK: Match spc6 end bead with Z1+ last bead
                if (abs(rx - z_beads(z_nbeads,1)) > tol) then
                   print *, "MISMATCH Frame", frame_idx, " End Bead"
                   print *, "spc6:", rx, ry, rz
                   print *, "Z1+ :", z_beads(z_nbeads,1), z_beads(z_nbeads,2), z_beads(z_nbeads,3)
                   stop
                endif
                
                ! 2. Process Chain Logic (Binning by COM)
                
                ! Get Box Bounds for PBC wrapping of the COM
                if (frame_idx > 0 .and. frame_idx <= nframe) then
                    b_idx = (frame_idx - 1) * 3
                    xlo = bb(b_idx + 1, 1); xhi = bb(b_idx + 1, 2)
                    ylo = bb(b_idx + 2, 1); yhi = bb(b_idx + 2, 2)
                    zlo = bb(b_idx + 3, 1); zhi = bb(b_idx + 3, 2)
                    
                    ! Apply PBC to Chain COM relative to dynamic box
                    ! (Matches your original logic for beads)
                    if (chain_cx < xlo) chain_cx = chain_cx + box
                    if (chain_cx > xhi) chain_cx = chain_cx - box
                    
                    if (chain_cy < ylo) chain_cy = chain_cy + box
                    if (chain_cy > yhi) chain_cy = chain_cy - box
                    
                    if (chain_cz < zlo) chain_cz = chain_cz + box
                    if (chain_cz > zhi) chain_cz = chain_cz - box
                    
                    ! Calculate Distance from Droplet COM
                    dist_com = sqrt((chain_cx - cm)**2 + (chain_cy - cm)**2 + (chain_cz - cm)**2)
                    
                    ! Determine Shell Index
                    shell_idx = int(dist_com / rad)
                    if (shell_idx > 6) shell_idx = 6
                    
                    ! 3. Update Histograms
                    chain_hist(shell_idx) = chain_hist(shell_idx) + 1.0d0
                    ent_hist(shell_idx)   = ent_hist(shell_idx) + dble(chain_ent_count)
                end if
            end if

        ! --- CASE C: ENTANGLEMENT (col5 == 1) ---
        elseif (int(col5) == 1) then
            ! Simply increment the counter for the current chain
            chain_ent_count = chain_ent_count + 1
            
        end if
    end do

    close(3)
    close(5)

    ! ==========================================
    ! 4. Write Output (Entanglements per Chain)
    ! ==========================================
    open(unit=4, file='ef.dat', status='replace')
    do i = 0, 6
        if (chain_hist(i) > 0.0d0) then
            write(4, '(I5, 2F15.6)') i, ent_hist(i)/chain_hist(i), chain_hist(i)
        else
            write(4, '(I5, 2F15.6)') i, 0.0d0, 0.0d0
        end if
    end do
    close(4)

end program calculation
