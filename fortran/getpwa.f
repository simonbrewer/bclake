!-------------------------------------------------------------------------------
! getpwa: Finds potential water areas from DEM and linear drainage map
!
! NB: expected draiage map is from r.watershed 
!
! "Provides the "aspect" for each cell measured CCW from East. 
! Multiplying positive values by 45 will give the direction in degrees 
! that the surface runoff will travel from that cell."
! 1: 45:  NE
! 2: 90:  N
! 3: 135: NW
! 4: 180: W
! 5: 225: SW
! 6: 270: S
! 7: 315: SE
! 8: 360: E
! 
!-------------------------------------------------------------------------------

      subroutine getpwa( m, n, dem, ldd, mask,
     >                   drain )
              
      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      integer ldd( m, n ) ! Linear drainage direction (-8 to 8)
      integer mask( m, n ) ! Binary mask (0/1)
 
      ! Internal variables
      integer i,j,k
      integer ii,jj
      integer edge,move
      integer offx(8),offy(8)
      double precision gridelev

      ! Outputs
      integer pwa( m, n ) ! Potential WA (0/1)
      double precision outelv( m, n ) ! Maximum water height in PWA
      integer drain( m, n )

      ! Set up parameters
      offy = (/+1,0,-1,-1,-1,0,+1,+1/) ! Offsets for movement
      offx = (/-1,-1,-1,0,+1,+1,+1,0/)

      !write(*,*) offx

      !-------------------------------------------------------------------------
      ! Main loop through grid cells

      ! Set all PWA to zero
      pwa(:,:) = 0
      ! Reset all outelv to zero
      outelv(:,:) = 0.
      ! Set drainage to zero
      drain(:,:) = 0

      do 10 i=20,20
      !do 10 i=1,m

      do 20 j=20,20
      !do 20 j=1,n

      if (mask(i,j).eq.1) then ! Check if cell is in lake watershed
        write(*,*) "start",i,j,dem(i,j),ldd(i,j)
        ! Tracking counters
        ii = i
        jj = j
        edge = 0
        k = 0
        gridelev = dem(i,j)

        do while (edge.eq.0)
          write(*,*) "current",ii,jj,dem(ii,jj),ldd(ii,jj)
          ! Test elevation
          if (dem(ii,jj).gt.gridelev) then
            pwa(i,j) = 1
          end if

          ! Set drainage output
          drain(ii,jj) = 1

          ! Get direction
          move = ldd(ii,jj)
          ii = ii + offx(move)
          jj = jj + offy(move)
          
          ! Test to see if we've reached the edge          
          if (ii.lt.1.or.ii.gt.m.or.jj.lt.1.or.jj.gt.n) then
            edge = 1
            write(*,*) "edge",ii,jj
          end if
          if (k.gt.1000) then
            exit
          end if

          k = k + 1

        end do
        

      endif ! Mask cell check

20    continue
10    continue

      end
