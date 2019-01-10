!-------------------------------------------------------------------------------
! SWAM: Version of the surface water area model that uses Coe 1998
! description 
!
! Includes code and equations from:
!
! Coe (1997). Simulating continental surface waters: an application to 
! Holocene Northern Africa. J. Clim. 10, 1680-1689
!
! Coe (1998). A linked global model of hydrological processes:
! simulation of modern rivers, lakes and wetlands. JGR, 103, 8885-8899
!
! Ver. 0.1 Takes codes from wca2d.f routine (basic flood model) and
! updates time and linear reservoir code
!-------------------------------------------------------------------------------

      subroutine swam_1t( m, n, dem, ldd, outelv,
     >                    mask, cella, celld, 
     >                    ppt, evap, runoff, baseflow,
     >                    dt, dtu, u, wvl, wse, war )
      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      integer ldd( m, n ) ! Linear drainge direction
      integer mask( m, n ) ! Binary mask (0/1)
      double precision cella( m, n ) ! Cell area (m2)
      double precision celld( m, n ) ! Distance between cells (m)
      double precision ppt( m, n ) ! PPT grid values (mm/d)
      double precision evap( m, n ) ! Evap grid values (mm/d)
      double precision runoff( m, n ) ! Runoff grid values (mm/d)
      double precision baseflow( m, n ) ! Baseflow grid values (mm/d)
      double precision dt ! Current time step (s)
      integer dtu ! Number of steps to run
      double precision u ! Effective velocity (m/s)

      ! Output variables
      double precision wvl( m, n ) ! water volume (m3)
      double precision wse( m, n ) ! water surface elevation (m)
      double precision war( m, n ) ! water area (0-1)

      ! Internal variables
      integer i,j,k ! Counters
      integer ii,jj,kk ! Counters
      double precision fin( m, n ) ! Inflow per time step (m3/s)
      double precision fout( m, n ) ! Outflow per time step (m3/s)
      integer edge,move
      integer offx(8),offy(8) ! Drainage direction codes 
      double precision foutcell,wel,wed,drain

      ! Set up parameters
      ! parameter(g = 9.81) ! (m/s-2)
      ! offx = (/-1,+1,0,0/) ! Offsets for von Neumann neighborhood
      ! offy = (/0,0,-1,+1/)
      ! offx = (/-1,0,+1,-1,+1,-1,0,+1/) ! Offsets for Moore neighborhood
      ! offy = (/-1,-1,-1,0,0,+1,+1,+1/)
      offy = (/+1,0,-1,-1,-1,0,+1,+1/) ! Offsets for r.watershed ldd
      offx = (/-1,-1,-1,0,+1,+1,+1,0/)

      !-------------------------------------------------------------------------
      ! Convert all forcings to m3/s
      do 920 i=1,m
      do 930 j=1,n

      if (mask(i,j).eq.1) then ! Check if cell is in lake watershed
        ! temporary variables to check for units
        runoff(i,j) = (runoff(i,j) * 1e-3) / (60*60*24)
        baseflow(i,j) = (baseflow(i,j) * 1e-3) / (60*60*24)
        ppt(i,j) = (ppt(i,j) * 1e-3) / (60*60*24)
        evap(i,j) = (evap(i,j) * 1e-3) / (60*60*24)

      end if ! Mask loop

930   continue
920   continue
        
      !-------------------------------------------------------------------------
      ! Main loop

      do 10 k=1,dtu

      !-------------------------------------------------------------------------
      ! Start by estimating inflow and outflow
      fin(:,:) = 0.0
      fout(:,:) = 0.0
      drain = 0.0
      do 20 i=1,m

      do 30 j=1,n

      if (mask(i,j).eq.1) then ! Check if cell is in lake watershed
        edge = 0 ! Set edge
        ! Find target downstream cell
        move = ldd(i,j)
        ii = i + offx(move)
        jj = j + offy(move)
 
        ! Test to see if we've left the basin
        if (ii.lt.1.or.ii.gt.m.or.jj.lt.1.or.jj.gt.n) then
          edge = 1
          !write(*,*) "edge",ii,jj
        end if

        ! Test to see if this is a coastal cell
        if (ldd(i,j).lt.0) then
          edge = 1
        end if

        ! Calculate outflow
        ! Get wei and wed
        wel = wse(i,j) + dem(i,j)
        wed = wse(ii,jj) + dem(ii,jj)
        if ((wse(i,j)-dem(i,j)).gt.0.) then
          fout(i,j) = max((wel-wed)*cella(i,j),0.) * u/celld(i,j)
        end if
        if (i.eq.60.and.j.eq.11) then
                write(*,*) k,i,j,edge,wel,wed,fout(i,j)
                write(*,*) wse(i,j),dem(i,j)
        end if

        ! If not an edge, then update inflow and outflow
        if (edge.eq.0) then
          fin(ii,jj) = fin(ii,jj) + fout(i,j)
        else ! If at edge account for water loss in drain reservoir
           drain = drain + fout(i,j)
        end if

      end if ! Mask loop

30    continue
20    continue

      !-------------------------------------------------------------------------
      ! Next update water volume per cell using linear reservoir
      do 120 i=1,m
      do 130 j=1,n

      if (mask(i,j).eq.1) then ! Check if cell is in lake watershed
        ! temporary variables to check for units
        ro_tmp = (runoff(i,j) * 1e-3) / ((60*60*24) * dt) * cella(i,j)
        bf_tmp = (baseflow(i,j) * 1e-3) / ((60*60*24) * dt) * cella(i,j)
        ppt_tmp = (ppt(i,j) * 1e-3) / ((60*60*24) * dt) * cella(i,j)
        evap_tmp = (evap(i,j) * 1e-3) / ((60*60*24) * dt) * cella(i,j)

        dwv = ( (runoff(i,j) + baseflow(i,j)) 
     >          * (1 - war(i,j)) * cella(i,j) ) +
     >        ( (ppt(i,j) - evap(i,j)) 
     >          * war(i,j) * cella(i,j) ) +
     >        ( fin(i,j) - fout(i,j) )

        wvl(i,j) = max(wvl(i,j) + dwv, 0.)
        wse(i,j) = wvl(i,j) / cella(i,j)

        if (wse(i,j).gt.0) then
          war(i,j) = min(wse(i,j), 1.)
        end if
 

      end if ! Mask loop

130   continue
120   continue
        
10    continue

      end
