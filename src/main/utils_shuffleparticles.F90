!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module utils_shuffleparticles
!
! Using only particle positions, this will shuffle those particles so
! that their density profile will match an input density profile.
! Current permitted profiles:
!   sphere-in-box (i.e. uniform density sphere in a uniform density medium)
!   BE sphere with a uniform background from an input array
!   an input array [UNTESTED]
!   using a different particle set
!
! :References: None
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: shuffleparticles

 ! identifiers for exact profile
 integer, parameter :: iuniform = 1
 integer, parameter :: iarray   = 2
 integer, parameter :: ipart    = 3
 ! additional parameters
 integer, parameter :: maxcellcache = 50000
 integer            :: ncall        =     1  ! to track the number of calls to this subroutine for debugging 
 logical            :: firstcall    = .true. ! logical to control debugging print statements
 private

contains
!-------------------------------------------------------------------
!
! Shuffle the particles to match an input profile
! Options and their specific inputs
!  Sphere-in-box:
!     rsphere: radius of the sphere
!     dsphere: density of the sphere
!     dmedium: density of the medium
!  Array with uniform background (designed for a BE sphere in a background medium)
!     ntab: array lengths
!     rtab: array radial positions
!     dtab: array with density positions
!     dmedium: density of the medium
!     dcontrast: density contrast between medium and edge of sphere
!  Array without a background
!     ntab: array lengths
!     rtab: array radial positions
!     dtab: array with density positions
!  Particle distribution:
!     xyzh_parent:  array containing the position of the parent particles
!     pmass_parent: mass of the parent particle
!     n_parent:     number of parent particles
! 
!-------------------------------------------------------------------
subroutine shuffleparticles(iprint,npart,xyzh,pmass,rsphere,dsphere,dmedium,ntab,rtab,dtab,dcontrast, &
                            xyzh_parent,pmass_parent,n_parent,is_setup,prefix)
 use dim,          only:maxneigh,maxp_hard
 use part,         only:vxyzu,divcurlv,divcurlB,Bevol,fxyzu,fext,alphaind
 use part,         only:gradh,rad,radprop,dvdx,rhoh,hrho
 use densityforce, only:densityiterate
 use linklist,     only:ncells,ifirstincell,set_linklist,get_neighbour_list,allocate_linklist
 use kernel,       only:cnormk,wkern,grkern,radkern2
 use boundary,     only:dxbound,dybound,dzbound,cross_boundary
 use kdtree,       only:inodeparts,inoderange
 use domain,       only:isperiodic
 integer,           intent(in)    :: iprint,npart
 real,              intent(in)    :: pmass
 real,              intent(inout) :: xyzh(:,:)
 integer, optional, intent(in)    :: ntab,n_parent
 real,    optional, intent(in)    :: rsphere,dsphere,dmedium,dcontrast,pmass_parent
 real,    optional, intent(in)    :: rtab(:),dtab(:)
 real,    optional, intent(inout) :: xyzh_parent(:,:)
 logical, optional, intent(in)    :: is_setup
 character(len=*) , optional, intent(in)    :: prefix
 integer      :: i,j,k,jm1,ip,icell,ineigh,idebug,ishift,nshiftmax,iprofile,nparterr,nneigh,ncross
 integer,save :: listneigh(maxneigh)
 real         :: stressmax,rmin,rmax,dr,dr1,dedge,dmed
 real         :: xi,yi,zi,hi,hi12,hi14,radi,coefi,rhoi,rhoi1,rij2,qi2,qj2,denom,rhoe,drhoe,err
 real         :: maggradi,maggrade,magshift,rinner,router,gradhi,gradhj
 real         :: dx_shift(3,npart),rij(3),runi(3),grrhoonrhoe(3),grrhoonrhoi(3),signg(3)
 real         :: errmin(3,2),errmax(3,2),errave(3,2),stddev(2,2),rtwelve(12),rnine(9),totalshift(3,npart)
 real         :: twoh21,kernsum,grrhoonrhoe_parent(3,maxp_hard)
 real,save    :: xyzcache(maxcellcache,4)
 real,save    :: xyzcache_parent(maxcellcache,4)
 logical      :: shuffle,at_interface,use_parent_h
 character(len=128) :: prefix0,fmt1,fmt2,fmt3
 !$omp threadprivate(xyzcache,xyzcache_parent,listneigh)

 !--Initialise free parameters
 idebug       =    1 ! 0 = off; 1=print initial & final distribution + errors; 2=print every step; 3=print every step with gradients
 nshiftmax    =  200 ! maximum number of shuffles/iterations
 !--Initialise remaining parameters
 rnine        =    0.
 rtwelve      =    0.
 use_parent_h = .true. ! to prevent compiler warnings

 !--Open debugging files and print the initial particle placements
 if (idebug > 0) then
    if (present(prefix)) then
       prefix0 = trim(prefix)//'_'
    else
       prefix0 = ''
    endif
    if (firstcall) then
       open(unit=332,file=trim(prefix0)//'shuffling_totalshift.dat')
       open(unit=333,file=trim(prefix0)//'shuffling_partposition.dat')
       open(unit=334,file=trim(prefix0)//'shuffling_error_all.dat')
       open(unit=335,file=trim(prefix0)//'shuffling_error_notinterface.dat')
       do j = 334,335
          write(j,"('#',13(1x,'[',i2.2,1x,a11,']',2x))") &
             1,'ishift',&
             2,'min rho err', &
             3,'ave rho err', &
             4,'max rho err', &
             5,'mn dd/d err', &
             6,'av dd/d err', &
             7,'mx dd/d err', &
             8,'min shift/h', &
             9,'ave shift/h', &
            10,'max shift/h', &
            11,'stddev_rho',  &
            12,'stddev_grho', &
            13,'ncall'
       enddo
       firstcall = .false.
    else
       open(unit=332,file=trim(prefix0)//'shuffling_totalshift.dat',        position='append')
       open(unit=333,file=trim(prefix0)//'shuffling_partposition.dat',      position='append')
       open(unit=334,file=trim(prefix0)//'shuffling_error_all.dat',         position='append')
       open(unit=335,file=trim(prefix0)//'shuffling_error_notinterface.dat',position='append')
       do j = 332,335
          write(j,'(a)') ' '
       enddo
    endif
    do i = 1,npart
       write(333,'(I18,1x,16(es18.10,1x),I18)') i,xyzh(1:3,i),rhoh(xyzh(4,i),pmass),rtwelve,ncall
    enddo
 endif

 !--Determine what has been passed in
 if (idebug > 0) then
    write(fmt1,'(a)') '(1x,a)'
    write(fmt2,'(a)') '(1x,a,I8)'
    write(fmt3,'(a)') '(1x,a,Es18.6)'
    if (present(ntab))         write(iprint,fmt2) 'Shuffling: optional variable is present: ntab: ', ntab
    if (present(rsphere))      write(iprint,fmt1) 'Shuffling: optional variable is present: rsphere'
    if (present(dsphere))      write(iprint,fmt1) 'Shuffling: optional variable is present: dsphere'
    if (present(dmedium))      write(iprint,fmt1) 'Shuffling: optional variable is present: dmedium'
    if (present(dcontrast))    write(iprint,fmt1) 'Shuffling: optional variable is present: dcontrast'
    if (present(rtab))         write(iprint,fmt2) 'Shuffling: optional variable is present: rtab: ', size(rtab)
    if (present(dtab))         write(iprint,fmt2) 'Shuffling: optional variable is present: dtab: ', size(rtab)
    if (present(xyzh_parent))  write(iprint,fmt1) 'Shuffling: optional variable is present: xyzh_parent'
    if (present(pmass_parent)) write(iprint,fmt3) 'Shuffling: optional variable is present: pmass_parent: ',pmass_parent
    if (present(n_parent))     write(iprint,fmt2) 'Shuffling: optional variable is present: n_parent: ',n_parent
    if (present(is_setup))     write(iprint,fmt1) 'Shuffling: optional variable is present: is_setup'
 endif

 !--Determine the exact profile
 iprofile = 0
 rmin     = 0.
 rmax     = 0.
 rinner   = 0.
 router   = 0.
 if (present(rsphere) .and. present(dsphere) .and. present(dmedium)) then
    iprofile = iuniform
    dedge    = dsphere
    dmed     = dmedium
    rmin     = rsphere - hrho(dsphere,pmass)
    rmax     = rsphere + hrho(dsphere,pmass)
    rinner   = rsphere - 4.*hrho(dsphere,pmass)
    router   = rsphere + 4.*hrho(dsphere,pmass)
    write(iprint,'(1x,a)') 'Shuffling: Shuffling to a uniform sphere-in-box profile'
 elseif (present(ntab) .and. present(rtab) .and. present(dtab) .and. present(dmedium) .and. present(dcontrast)) then
    iprofile = iarray
    dedge    = dmedium*dcontrast
    dmed     = dmedium
    if (.false.) then
       rmin     = rtab(ntab) - hrho(dedge,pmass)
       rmax     = rtab(ntab) + hrho(dedge,pmass)
       rinner   = rtab(ntab) - 4.*hrho(dedge,pmass)
       router   = rtab(ntab) + 4.*hrho(dedge,pmass)
       ! reset dedge to match the density at rmin
       do j = 1,ntab
          if (rtab(j) < rmin) dedge = dtab(j)
       enddo
    else
       rmin     = rtab(ntab)
       rmax     = rtab(ntab) + 2.*hrho(dedge,pmass)
       rinner   = rtab(ntab) - 2.*hrho(dedge,pmass)
       router   = rtab(ntab) + 6.*hrho(dedge,pmass)
    endif
    write(iprint,'(1x,a)') 'Shuffling: Shuffling to non-uniform sphere-in-box profile (eg a BE sphere in a background medium)'
 elseif (present(ntab) .and. present(rtab) .and. present(dtab)) then
    iprofile = iarray
    rmin     = 2.0*rtab(ntab) ! setting to a large value to be always called
    rmax     = 2.0*rtab(ntab) ! setting to a large value that will never be called
    write(iprint,*) 'Shuffling to a non-uniform sphere'
    write(iprint,*) 'This is untested without a background, so aborting until properly tested'
    return
 elseif (present(xyzh_parent) .and. present(pmass_parent) .and. present(n_parent)) then 
    iprofile = ipart
    if (pmass > pmass_parent) then
       use_parent_h = .true.
    else
       use_parent_h = .false.
    endif
    write(iprint,'(1x,a)') 'Shuffling: Shuffling to a parent array'
#ifdef SPLITTING
    grrhoonrhoe_parent = 0.
!$omp parallel do default (none) &
!$omp shared(xyzh_parent,pmass_parent,grrhoonrhoe_parent,n_parent) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
!$omp private(i,j,xi,yi,zi,hi,hi12,rij,rij2,runi,qi2,qj2,hi14,rhoi1) &
!$omp private(denom,gradhi,gradhj)
    over_pparts: do i=1,n_parent
       xi = xyzh_parent(1,i)
       yi = xyzh_parent(2,i)
       zi = xyzh_parent(3,i)
       hi = xyzh_parent(4,i)
       gradhi = xyzh_parent(5,i)
       hi12  = 1.0/(hi*hi)
       hi14  = hi12*hi12
       rhoi1 = 1.0/rhoh(hi,pmass_parent)
       ! calculate (grad rho) / rho for particle i
       over_pneighbours: do j = 1,n_parent
          if (i==j) cycle over_pneighbours

          rij(1) = xyzh_parent(1,j) - xi
          rij(2) = xyzh_parent(2,j) - yi
          rij(3) = xyzh_parent(3,j) - zi
          gradhj = xyzh_parent(5,j)

#ifdef PERIODIC
          if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
          if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
          if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
#endif
          rij2 = dot_product(rij,rij)
          runi = rij/sqrt(rij2)
          qi2  = rij2*hi12
          qj2  = rij2/xyzh_parent(4,j)**2
          if (qi2 < radkern2) then
             grrhoonrhoe_parent(:,i) = grrhoonrhoe_parent(:,i) - runi*grkern(qi2,sqrt(qi2))*cnormk*gradhi*hi14*rhoi1
          endif
          if (qj2 < radkern2) then
             denom = xyzh_parent(4,j)**4 * rhoh(xyzh_parent(4,j),pmass_parent)
             grrhoonrhoe_parent(:,i) = grrhoonrhoe_parent(:,i) - runi*grkern(qj2,sqrt(qj2))*cnormk*gradhj/denom
          endif
       enddo over_pneighbours
       grrhoonrhoe_parent(:,i) = pmass_parent*grrhoonrhoe_parent(:,i)
    enddo over_pparts
!$omp end parallel do
#endif 
 else
    write(iprint,'(1x,a)') 'Shuffling: Density profile undefined.  Aborting.'
    return
 endif
 dr = rmax - rmin
 if (dr > 0.) then
    dr1 = 1./dr
 else
    dr1 = 0.
 endif

 !--initialise memory for linklist
 if (present(is_setup)) then
    if (is_setup) call allocate_linklist()
 endif

 !--Shuffle particles
 ishift     = 0
 shuffle    = .true.
 totalshift = 0.
 do while (shuffle .and. ishift < nshiftmax)
    ! update densities
    call set_linklist(npart,npart,xyzh,vxyzu)
    call densityiterate(2,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                               fxyzu,fext,alphaind,gradh,rad,radprop,dvdx)

    ! initialise variables for this loop
    dx_shift = 0.
    errmin   = huge(errmin)
    errmax   = 0.
    errave   = 0.
    stddev   = 0.
    nparterr = 0

    ! determine how much to shift by
!$omp parallel do default (none) &
!$omp shared(xyzh,npart,pmass,dx_shift,gradh,rmin,rmax,iprofile,dedge,dmed,dr1,idebug) &
!$omp shared(use_parent_h,n_parent,xyzh_parent,grrhoonrhoe_parent,totalshift) &
!$omp shared(rtab,dtab,ntab,rinner,router,inodeparts,inoderange,ifirstincell,ncells,ncall) &
#ifdef PERIODIC
!$omp shared(dxbound,dybound,dzbound) &
#endif
!$omp private(i,j,k,ip,ineigh,icell,jm1,xi,yi,zi,hi,hi12,rij,rij2,runi,qi2,qj2,hi14,rhoi,rhoi1,radi,coefi) &
!$omp private(grrhoonrhoe,grrhoonrhoi,signg,rhoe,drhoe,denom,nneigh) &
!$omp private(err,maggradi,maggrade,magshift,at_interface) &
!$omp private(twoh21,kernsum) & ! unique to splitting
!$omp reduction(min: errmin) &
!$omp reduction(max: errmax) &
!$omp reduction(+:   errave,stddev,nparterr)
    over_cells: do icell=1,int(ncells)
       k = ifirstincell(icell)

       ! Skip empty cells AND inactive cells
       if (k <= 0) cycle over_cells

       ! Get the neighbour list and fill the cell cache
       call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)

       ! Loop over particles in the cell
       over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
          i  = inodeparts(ip)
          xi = xyzh(1,i)
          yi = xyzh(2,i)
          zi = xyzh(3,i)
          hi = xyzh(4,i)
          hi12  = 1.0/(hi*hi)
          hi14  = hi12*hi12
          rhoi  = rhoh(hi,pmass)
          rhoi1 = 1.0/rhoi
          radi  = sqrt(xi*xi + yi*yi + zi*zi)
          !coefi = 0.0225*hi*hi
          !coefi = 0.075*hi*hi
          coefi = 0.25*hi*hi
          grrhoonrhoi = 0.

          ! calculate (grad rho) / rho for particle i
          over_neighbours: do ineigh = 1,nneigh
             j = abs(listneigh(ineigh))
             if (i==j) cycle over_neighbours

             rij(1) = xyzh(1,j) - xi
             rij(2) = xyzh(2,j) - yi
             rij(3) = xyzh(3,j) - zi
#ifdef PERIODIC
             if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
             if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
             if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
#endif
             rij2 = dot_product(rij,rij)
             if (rij2 < epsilon(rij2)) then
                if (idebug > 0) print*, 'pairing may be occurring between ', i, j
                cycle
             endif
             runi = rij/sqrt(rij2)
             qi2  = rij2*hi12
             qj2  = rij2/xyzh(4,j)**2
             if (qi2 < radkern2) then
                grrhoonrhoi = grrhoonrhoi - runi*grkern(qi2,sqrt(qi2))*cnormk*gradh(1,i)*hi14*rhoi1
             endif
             if (qj2 < radkern2) then
                denom = xyzh(4,j)**4 * rhoh(xyzh(4,j),pmass)
                grrhoonrhoi = grrhoonrhoi - runi*grkern(qj2,sqrt(qj2))*cnormk*gradh(1,j)/denom
             endif
          enddo over_neighbours

          ! calculate the exact (grad rho) / rho
          grrhoonrhoe = 0.0
          drhoe       = 0.0
          rhoe        = 0.0
          if (iprofile==ipart) then
             ! determine the kernel-weighted average grad rho / rho for the parent distribution
             ! need to make this more efficient
             twoh21 = 1./(2.0*xyzh(4,i))**2
             if (xyzh(4,i) < epsilon(xyzh(4,i))) print*, 'fucking i'
             kernsum     = 0.
             over_parents: do j = 1,n_parent
                rij(1) = xyzh_parent(1,j) - xi
                rij(2) = xyzh_parent(2,j) - yi
                rij(3) = xyzh_parent(3,j) - zi
                if (use_parent_h) twoh21 = 1./(2.0*xyzh_parent(4,j))**2
                if (xyzh_parent(4,j) < epsilon(xyzh_parent(4,j))) print*, 'fucking j'
#ifdef PERIODIC
                if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
                if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
                if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
#endif
                rij2 = dot_product(rij,rij)
                qi2  = rij2*twoh21
                if (qi2 < radkern2) then
                   grrhoonrhoe = grrhoonrhoe - grrhoonrhoe_parent(:,j)*wkern(qi2,sqrt(qi2))
                   kernsum     = kernsum + wkern(qi2,sqrt(qi2))
                endif
             enddo over_parents
             if (kernsum > 0.) then
                grrhoonrhoe = grrhoonrhoe/kernsum
             else
                print*, 'Crap... no neighbours!'
             endif
          else
             if (radi < rmin) then
                if (iprofile==iuniform) then
                   rhoe = dedge
                elseif (iprofile==iarray) then
                   j = 1
                   do while (radi > rtab(j) .and. j < ntab)
                      j = j + 1
                   enddo
                   if (j==1) then
                      jm1 = j + 1
                   else
                      jm1 = j - 1
                   endif
                   drhoe = (dtab(j) - dtab(jm1))/(rtab(j) - rtab(jm1))
                   rhoe  = dtab(jm1) +  drhoe*(radi - rtab(jm1))
                endif
             elseif (rmin <= radi .and. radi < rmax) then
                rhoe  = (dedge*(rmax - radi) + dmed*(radi - rmin))*dr1
                drhoe = (dmed - dedge)*dr1
             else  !if (radi > rmax) then
                rhoe = dmed
             endif
             if (rhoe > 0. .and. radi > 0.) then
                grrhoonrhoe(1) = drhoe*xi/(rhoe*radi)
                grrhoonrhoe(2) = drhoe*yi/(rhoe*radi)
                grrhoonrhoe(3) = drhoe*zi/(rhoe*radi)
            endif
          endif
         ! print*, i,grrhoonrhoe,coefi

          ! shift the particles
          grrhoonrhoi  = coefi*pmass*grrhoonrhoi   ! multiply in the correct coefficients for the i'th particle
          grrhoonrhoe  = coefi      *grrhoonrhoe   ! multiply in the correct coefficients for the exact value
          signg        = 1.0                       ! determine whether to add or subtract the terms
          if (iprofile/=ipart) then
             if (xi*grrhoonrhoi(1) > xi*grrhoonrhoe(1)) signg(1) = -1.0
             if (yi*grrhoonrhoi(2) > yi*grrhoonrhoe(2)) signg(2) = -1.0
             if (zi*grrhoonrhoi(3) > zi*grrhoonrhoe(3)) signg(3) = -1.0
          endif
          dx_shift(:,i) = grrhoonrhoi + signg*grrhoonrhoe

          ! limit the particle shift to 0.5h
          magshift  = dot_product(dx_shift(1:3,i),dx_shift(1:3,i))
          if (magshift > 0.25*hi*hi) dx_shift(:,i) = 0.5*hi*dx_shift(:,i)/sqrt(magshift)
          totalshift(:,i) = totalshift(:,i) - dx_shift(:,i)

          ! debugging statements
          if (idebug > 0) then
             if (radi < rinner .or. radi > router .or. iprofile==ipart) then
                at_interface = .false.
                nparterr     = nparterr + 1
             else
                at_interface = .true.
             endif
             ! calculate error in density
             err       = abs(1.0-rhoe*rhoi1)
             errmin(1,1) = min(errmin(1,1),err)
             errmax(1,1) = max(errmax(1,1),err)
             errave(1,1) = errave(1,1) + err
             if (.not. at_interface) then
                errmin(1,2) = min(errmin(1,2),err)
                errmax(1,2) = max(errmax(1,2),err)
                errave(1,2) = errave(1,2) + err
             endif
             ! calculate error in grad rho / rho
             maggrade  = dot_product(grrhoonrhoe,grrhoonrhoe)
             maggradi  = dot_product(grrhoonrhoi,grrhoonrhoi)
             if (maggrade > 0.) then
                err    = abs(1.0-sqrt(maggradi/maggrade))
             else
                err    = sqrt(maggradi)
             endif
             errmin(2,1) = min(errmin(2,1),err)
             errmax(2,1) = max(errmax(2,1),err)
             errave(2,1) = errave(2,1) + err
             if (.not. at_interface) then
                errmin(2,2) = min(errmin(2,2),err)
                errmax(2,2) = max(errmax(2,2),err)
                errave(2,2) = errave(2,2) + err
             endif
             ! calculate the magnitude of the shift
             magshift  = dot_product(dx_shift(1:3,i),dx_shift(1:3,i))
             err       = sqrt(magshift)/hi
             errmin(3,1) = min(errmin(3,1),err)
             errmax(3,1) = max(errmax(3,1),err)
             errave(3,1) = errave(3,1) + err
             if (.not. at_interface) then
                errmin(3,2) = min(errmin(3,2),err)
                errmax(3,2) = max(errmax(3,2),err)
                errave(3,2) = errave(3,2) + err
             endif
             ! calculate the standard deviation compared to the exact values
             stddev(1,1) = stddev(1,1) + (rhoe - rhoi)**2
             stddev(2,1) = stddev(2,1) + (maggrade - maggradi)**2
             if (.not. at_interface) then
                stddev(1,2) = stddev(1,2) + (rhoe - rhoi)**2
                stddev(2,2) = stddev(2,2) + (maggrade - maggradi)**2
             endif
             if (idebug == 3) then
!$omp critical
                write(333,'(I18,1x,16(es18.10,1x),I18)') i,xyzh(1:3,i),rhoh(hi,pmass),&
                dx_shift(:,i),grrhoonrhoi,grrhoonrhoe,magshift,maggradi,maggrade,ncall
!$omp end critical
             endif
          endif
       enddo over_parts
    enddo over_cells
!$omp end parallel do

    ! shift particles
    xyzh(1:3,1:npart) = xyzh(1:3,1:npart) - dx_shift(:,1:npart)

    ! write diagnostics to file
    if (idebug > 0) then
       stddev(:,1) = sqrt(stddev(:,1)/npart)
       stddev(:,2) = sqrt(stddev(:,1)/nparterr)

       write(334,'(I18,1x,11(es18.10,1x),I18)') ishift &
                          ,errmin(1,1),errave(1,1)/npart,errmax(1,1) &
                          ,errmin(2,1),errave(2,1)/npart,errmax(2,1) &
                          ,errmin(3,1),errave(3,1)/npart,errmax(3,1) &
                          ,stddev(1:2,1),ncall
       write(335,'(I18,1x,11(es18.10,1x),I18)') ishift &
                          ,errmin(1,2),errave(1,2)/nparterr,errmax(1,2) &
                          ,errmin(2,2),errave(2,2)/nparterr,errmax(2,2) &
                          ,errmin(3,2),errave(3,2)/nparterr,errmax(3,2) &
                          ,stddev(1:2,2),ncall
       if (idebug==2) then
          do i = 1,npart
             write(333,'(I18,1x,16(es18.10,1x),I18)') i,xyzh(1:3,i),rhoh(xyzh(4,i),pmass),&
             dx_shift(:,i),rnine,ncall
          enddo
       endif
       if (idebug/=1) write(333,*) ' '
    endif

    ! Determine stopping criteria been met [DEVELOP THIS]
    if (.false.) shuffle = .false.

    ! update counter
    ishift = ishift + 1
 enddo

 ! re-adjust particles to ensure none crossed the boundary
#ifdef PERIODIC
 ncross = 0
!$omp parallel do default (none) &
!$omp shared(npart,xyzh,isperiodic) &
!$omp private(i) &
!$omp reduction(+:ncross)
 do i = 1,npart
    call cross_boundary(isperiodic,xyzh(:,i),ncross)
 enddo
!$omp end parallel do
 if (idebug > 0 .and. ncross > 0) write(iprint,'(1x,a,I6)') 'Shuffling: number of particles that cross the boundary: ',ncross
#endif

 ! final debugging print-statements
 if (idebug > 0) then
    do i = 1,npart
       write(332,'(I18,1x,7(es18.10,1x),I18)') i,xyzh(1:4,i),totalshift(1:3,i),ncall
    enddo
    if (idebug == 1) then
       do i = 1,npart
          write(333,'(I18,1x,16(es18.10,1x),I18)') i,xyzh(1:3,i),rhoh(xyzh(4,i),pmass),rtwelve,ncall
       enddo
    endif
    do i = 332,335
       close(i)
    enddo
 endif
 write(iprint,'(1x,2(a,I6))') 'Shuffling: completed with ',ishift,' iterations on call number ',ncall
 ncall = ncall + 1

end subroutine shuffleparticles
!-------------------------------------------------------------------
end module utils_shuffleparticles
