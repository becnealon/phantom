!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module splitmergeutils
!
! None
!
! :References: None
!
! :Owner: Bec Nealon
!
! :Runtime parameters: None
!
! :Dependencies: icosahedron, kernel, part, random
!
 implicit none
! integer, public :: rand_type = 0  ! no randomness: first out of every 13 particles is the parent  (A)
! integer, public :: rand_type = 1  ! 1 of every 13 randomly selected for the parent (B)
! integer, public :: rand_type = 2  ! dice rolled each time (C)
! integer, public :: rand_type = 3  ! Use the tree in the nice code.  gas then splits (D)
! integer, public :: rand_type = 4  ! Use the tree in the hacked code (E)
 integer, public :: rand_type = 5  ! Use the tree in the nice code.  gas & splits simultaneously (F)
 logical, public :: centre_particle = .true.

 contains

!--------------------------------------------------------------------------
!+
!  splits a particle into nchild particles
!+
!--------------------------------------------------------------------------
subroutine split_a_particle(nchild,iparent,xyzh,vxyzu,npartoftype,lattice_type,ires,ichildren)
 use icosahedron, only:pixel2vector,compute_corners,compute_matrices
 use part,        only:copy_particle_all,igas,isplit,set_particle_type,kill_particle,iphase
 use random,      only:ran2
 use physcon,     only:pi
 use vectorutils, only:rotatevec
 integer, intent(in)    :: nchild,iparent,lattice_type,ires,ichildren
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 integer :: j,iseed,ichild,anotherseed
 real    :: dhfac,dx(3),sep,geodesic_R(0:19,3,3), geodesic_v(0:11,3),hchild
 real    :: gamma,beta

 if (lattice_type == 0) then
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 endif
 ! initialise random number generator
 iseed       = -6542
 anotherseed = -2845

 if (centre_particle) then
    dhfac = 1./(nchild+1)**(1./3.)
 else
    dhfac = 1./(nchild)**(1./3.)
 endif
 hchild = xyzh(4,iparent)*dhfac
 sep    = 0.41*xyzh(4,iparent)
 ichild = 0
 beta   = acos(2.0*ran2(iseed)-1.0)
 gamma  = 2.0*pi*ran2(anotherseed)

 do j=0,nchild-1
    ichild = ichild + 1
    ! copy properties
    call copy_particle_all(iparent,ichildren+ichild,.true.)

    ! adjust the position
    if (lattice_type == 0) then
       call pixel2vector(j,ires,geodesic_R,geodesic_v,dx)
    else
       call sample_kernel(iseed,dx)
    endif

    ! rotate the particles on the sphere (consistent for each sphere)
    call rotatevec(dx,(/1.,0.,0./),gamma)
    call rotatevec(dx,(/0.,1.0,0./),beta)
    xyzh(1:3,ichildren+ichild) = xyzh(1:3,iparent) + sep*dx(:)
    xyzh(4,  ichildren+ichild) = xyzh(4,  iparent)*dhfac
    call set_particle_type(ichildren+ichild,isplit)
 enddo
 npartoftype(isplit) = npartoftype(isplit) + nchild

 if (centre_particle) then
    !--fix parent
    xyzh(4,iparent) = xyzh(4,iparent)*dhfac
    call set_particle_type(iparent,isplit)
!    iphase(iparent) = -abs(iphase(iparent))

    !-- tidy up particle types
    npartoftype(igas)   = npartoftype(igas)   - 1
    npartoftype(isplit) = npartoftype(isplit) + 1

 else
    call kill_particle(iparent,npartoftype(:))
 endif

end subroutine split_a_particle

subroutine sample_kernel(iseed,dx)
 use random, only:gauss_random,get_random_pos_on_sphere
 integer, intent(inout) :: iseed
 real :: dx(3),r

 r = 3.
 do while (abs(r) > 2.)
    r = gauss_random(iseed)
 enddo
 dx = get_random_pos_on_sphere(iseed)
 dx = r*dx

end subroutine sample_kernel

!-----------------------------------------------------------------------
!+
! routine that shuffles the particle positions slightly
! (see Vacondio et al. 2013 Equation 11)
! (currently mostly used in testing, not sure if it's important yet)
!+
!-----------------------------------------------------------------------
subroutine shift_particles_Vacondio(npart,xyzh,vxyzu,deltat,beta,shifts)
 use kernel, only:radkern2
 use boundary, only: dxbound,dybound,dzbound
 use part,     only: periodic
 real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)    :: deltat,beta
 integer, intent(in) :: npart
 real, intent(out)   :: shifts(3,npart)
 integer             :: i,j,neighbours
 real                :: rnaught,rij2,dr3,vel2,vmax
 real                :: q2,rij(3),rsum(3)
 logical             :: jhw = .true.

 !if (jhw) print*, 'using minor modifications of jhw'
 vmax = tiny(vmax)
 vel2 = 0.
 shifts = 0.

 do i = 1,npart
    rnaught = 0.
    neighbours = 0
    rsum = 0.
    vel2 = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    if (jhw) then
       vmax = vel2
    else
       if (vel2 > vmax) vmax = vel2
    endif

    over_npart: do j = 1,npart
       if (i == j) cycle over_npart
       rij = xyzh(1:3,j) - xyzh(1:3,i)
       if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
          if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
          if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
          if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
       endif

       rij2 = dot_product(rij,rij)
       q2 = rij2/(xyzh(4,i)*xyzh(4,i))

       if (q2 < radkern2) then
          neighbours = neighbours + 1
          rnaught = rnaught + sqrt(rij2)
          if (jhw) then
             vel2 = dot_product(vxyzu(1:3,j),vxyzu(1:3,j))
             vmax = max(vel2,vmax)
          endif
       endif

       dr3 = 1./(rij2**1.5)
       rsum = rsum + (rij*dr3)
    enddo over_npart

    rnaught = rnaught/neighbours
    shifts(:,i) = beta*rnaught*rnaught*deltat*rsum
    if (jhw) shifts(:,i) = shifts(:,i)*sqrt(vmax)/neighbours
 enddo
 ! Final scaling
 if (.not. jhw) shifts = shifts*sqrt(vmax)
 ! Apply shifts - turned off for the moment
 !xyzh(1:3,1:npart) = xyzh(1:3,1:npart) - shifts

end subroutine shift_particles_Vacondio

!-----------------------------------------------------------------------
!+
! calculates particle shifts according to WVT method from
! Diehls et al 2015
!+
!-----------------------------------------------------------------------
subroutine shift_particles_WVT(npart,xyzh,h0,mu)
 use kernel, only:radkern2,cnormk,kernel_softening
 use part,   only:isdead_or_accreted,iactive,iphase,periodic,rhoh,massoftype,isplit
 use boundary, only: dxbound,dybound,dzbound,xmin,xmax,ymin,zmin,zmax,ymax
 real, intent(inout)    :: xyzh(:,:),h0(:)
 real, intent(in)    :: mu
 integer, intent(in) :: npart
 integer             :: i,j,k,k1,k3,kmax
 real                :: rij2,f,hij,hi12,hi14,f_const,xi,yi,zi,hi,xj,yj,zj,hj,fmax,hcoef,dummy,fcoef
 real                :: q,q2,rij(3)
 real                :: shifts(6,npart),maxshift
 logical             :: really_far,all_r
 logical             :: Diehl = .true.              ! use the Diehl+2015 equation; else use JW's
 logical             :: local_debug = .false.         ! print out a useful fort file to watch the particle shifting
 logical             :: dist_plus_const = .false.    ! r2 = r2 + (0.01*hi)**2; else r2 = r2
 logical             :: local_kernel = .false.       ! calculate forces from the kernel; else f \propto 1/r2
 logical             :: indiv_scalar = .false.       ! scale the max distance on a per-particle basis, not globally


 if (Diehl) then
    kmax  = 1
    all_r = .false.
 else
    kmax  = 2
    all_r = .true.
 endif
 shifts  = 0.
 f_const = 0.25 ! to ensure that f = 0 at 2h
 fmax    = 0.
 hcoef   = 1. ! kinda like mu, but constant
 do k = 1,kmax
    if (k==1) then
       k1 = 1
       k3 = 3
    else
       k1 = 4
       k3 = 6
    endif
    fcoef = 1.
    if (fmax > 1.) fcoef = 1./fmax
!$omp parallel default(none) &
!$omp shared(npart,iphase,xyzh,h0,shifts,k,k1,k3,f_const,dxbound,dybound,dzbound,hcoef) &
!$omp shared(Diehl,dist_plus_const,local_kernel,local_debug,xmin,indiv_scalar,all_r) &
!$omp private(xi,yi,zi,hi,xj,yj,zj,hj,hi12,hi14,rij,rij2,q2,q,hij,f,dummy) &
!$omp firstprivate(fcoef) &
!$omp reduction(max:fmax)
!$omp do
    do i = 1,npart !for each active particle
       ! CHECK ITS ACTIVE
       if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
          hi = h0(i) !xyzh(4,i)
          if (k==1) then
             xi = xyzh(1,i)
             yi = xyzh(2,i)
             zi = xyzh(3,i)
          else
             if (indiv_scalar) then
                fcoef = sqrt(dot_product(shifts(1:3,i),shifts(1:3,i)))
                if (fcoef > 1.) then
                   fcoef = 1./fcoef
                else
                   fcoef = 1.
                endif
             endif
             xi = xyzh(1,i)+shifts(1,i)*0.1*hcoef*hi*fcoef
             yi = xyzh(2,i)+shifts(2,i)*0.1*hcoef*hi*fcoef
             zi = xyzh(3,i)+shifts(3,i)*0.1*hcoef*hi*fcoef
             if (xi < xmin) print*, 'crap',xi
          endif
          hi12 = 1.0/(hi*hi)
          hi14 = hi12*hi12
          over_npart: do j = 1,npart
             if (i == j) cycle over_npart
             hj = h0(j)
             if (k==1) then
                xj = xyzh(1,j)
                yj = xyzh(2,j)
                zj = xyzh(3,j)
             else
                if (indiv_scalar) then
                   fcoef = sqrt(dot_product(shifts(1:3,j),shifts(1:3,j)))
                   if (fcoef > 1.) then
                      fcoef = 1./fcoef
                   else
                      fcoef = 1.
                   endif
                endif
                xj = xyzh(1,j)+shifts(1,j)*0.1*hcoef*hj*fcoef
                yj = xyzh(2,j)+shifts(2,j)*0.1*hcoef*hj*fcoef
                zj = xyzh(3,j)+shifts(3,j)*0.1*hcoef*hj*fcoef
             endif
             rij(1) = xj - xi
             rij(2) = yj - yi
             rij(3) = zj - zi
             if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
                if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
                if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
                if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
             endif
             rij2 = dot_product(rij,rij)
             q2   = rij2*hi12    ! 0 < q < 2 therefore 0 < q2 < 4

             if (q2 < radkern2 .or. all_r) then  ! Might want to play with what h / combination of h we use for q
                if (Diehl) then
                   hij = 0.5*(hi + hj)
                   f   = (hij/(sqrt(rij2) + epsilon(rij2)))**2 - f_const       ! without constant and epsilon,  inf > f > 0.25
                   shifts(1:3,i) = shifts(1:3,i) + hi*f*rij/sqrt(rij2)
                else
                   if (dist_plus_const) rij2 = rij2 + (0.01*hi)**2
                   if (local_kernel) then
                      q = sqrt(q2)
                      call kernel_softening(q2,q,dummy,f)
                      f = -f*cnormk*hi14/(hi*hj)
                   else
                      f = -1./(hi*hj*rij2)  ! The h's are proxy for density which is a proxy for mass
!                      f = -1./rij2  ! [no difference between this and our case above]
                   endif
                   shifts(k1:k3,i) = shifts(k1:k3,i) + f*rij/sqrt(rij2)
                endif
             endif
          enddo over_npart
          fmax = max(fmax,sqrt(dot_product(shifts(k1:k3,i),shifts(k1:k3,i))))
       endif
    enddo
!$omp enddo
!$omp end parallel
    if (local_debug) print*, mu,fmax
 enddo
 ! Final form
 if (Diehl) then
    shifts = shifts*mu
 else
    shifts(1:3,:) = 0.5*(shifts(1:3,:)+shifts(4:6,:))
    if (fmax < 1.0) then
       fcoef = 1.0
    else
       fcoef = 1.0/fmax
    endif
!$omp parallel default(none) &
!$omp shared(npart,h0,shifts,hcoef,fmax,indiv_scalar) &
!$omp firstprivate(fcoef) &
!$omp private(i)
!$omp do
    do i = 1,npart
       if (indiv_scalar) then
          fcoef = sqrt(dot_product(shifts(1:3,i),shifts(1:3,i)))
          if (fcoef > 1.) then
             fcoef = 1./fcoef
          else
             fcoef = 1.
          endif
       endif
       shifts(1:3,i) = -0.1*hcoef*h0(i)*shifts(1:3,i)*fcoef
    enddo
!$omp enddo
!$omp end parallel
 endif

 ! Check, do any of the shifts correspond to more than the box size?
 maxshift = maxval(abs(shifts))
 !if (maxshift > 0.5) print*,'WARNING SHIFTS ARE LARGE',maxshift

 ! Now apply the shifts and update for periodicity
 xyzh(1:3,1:npart) = xyzh(1:3,1:npart) - shifts(1:3,1:npart)
 if (periodic) then
    really_far = .true.
    do while (really_far)
       !print*, 're-checking periodic at mu=',mu
       really_far = .false.
       call shift_for_periodicity(npart,xyzh)
       do i = 1,npart
          if (xyzh(1,i) < xmin) really_far = .true.
          if (xyzh(2,i) < ymin) really_far = .true.
          if (xyzh(3,i) < zmin) really_far = .true.
          if (xyzh(1,i) > xmax) really_far = .true.
          if (xyzh(2,i) > ymax) really_far = .true.
          if (xyzh(3,i) > zmax) really_far = .true.
          if (xyzh(4,i) > 0.5) print*, 'big h at i: ',xyzh(4,i),i
       enddo
    enddo
 endif

 if (local_debug) then
    write(333,*) ' '
    do i = 1,npart
       write(333,*) xyzh(1:4,i),shifts(1:3,i)/xyzh(4,i),rhoh(xyzh(4,i),massoftype(isplit))
    enddo

    if (npart < 1000 .and. mu < 2.d-5) then
       do i = 1,npart
          do j = 1,npart
             write(666,*) i,j,sqrt(dot_product((xyzh(1:3,i)-xyzh(1:3,j)),(xyzh(1:3,i)-xyzh(1:3,j))))
          enddo
          write(666,*) ' '
          write(777,*) i,xyzh(1:3,i)
       enddo
       write(777,*) ' '
    endif
 endif

end subroutine shift_particles_WVT

!-----------------------------------------------------------------------
!+
! merges nchild particles in one parent particle
! parent properties averaged from children
!+
!-----------------------------------------------------------------------
subroutine fancy_merge_into_a_particle(nchild,ichildren,mchild,npart, &
                                       xyzh,vxyzu,npartoftype,iparent)
 use kernel, only:get_kernel,cnormk,radkern
 use part,   only:copy_particle_all,hfact,kill_particle,igas,isplit
 integer, intent(in)    :: nchild,ichildren(nchild),iparent
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)   :: mchild
 integer :: i,j,ichild
 real    :: h1,h31,rij_vec(3)
 real    :: qij,rij,wchild,grkernchild,rho_parent

 !-- copy properties from first child
 call copy_particle_all(ichildren(1),iparent,.true.)

 !-- positions and velocities from centre of mass
 xyzh(1:3,iparent) = 0.
 vxyzu(1:3,iparent) = 0.

 do i=1,nchild
    ichild = ichildren(i)
    xyzh(1:3,iparent) = xyzh(1:3,iparent) + xyzh(1:3,ichild)
    vxyzu(1:3,iparent) = vxyzu(1:3,iparent) + vxyzu(1:3,ichild)
 enddo
 xyzh(1:3,iparent) = xyzh(1:3,iparent)/nchild
 vxyzu(1:3,iparent) = vxyzu(1:3,iparent)/nchild

!-- calculate density at the parent position from the children
!-- procedure described in Vacondio et al. 2013, around Eq 21
 rho_parent = 0.
 do j=1,npart
    h1 = 1./xyzh(4,j)
    h31 = h1**3

    rij_vec = xyzh(1:3,iparent) - xyzh(1:3,j)

    rij = sqrt(dot_product(rij_vec,rij_vec))
    qij = rij*h1

    wchild = 0.
    if (qij < radkern) call get_kernel(qij*qij,qij,wchild,grkernchild)

    rho_parent = rho_parent + (mchild*wchild*cnormk*h31)
 enddo

!-- smoothing length from density
 xyzh(4,iparent) = hfact*(nchild*mchild/rho_parent)**(1./3.)

 ! kill the useless children
 do i=1,nchild
    call kill_particle(ichildren(i),npartoftype(:))
 enddo

 !-- tidy up particle types
 npartoftype(igas) = npartoftype(igas) + 1

end subroutine fancy_merge_into_a_particle

!-----------------------------------------------------------------------
!+
! merges nchild particles in one parent particle
! chooses the first particle to be the de-facto parent
! and discards the rest (super fast)
!+
!-----------------------------------------------------------------------
subroutine fast_merge_into_a_particle(nchild,ichildren,npart, &
           xyzh,vxyzu,npartoftype,iparent)
 use part,   only:copy_particle_all,kill_particle,igas,isplit,set_particle_type
 use random, only:ran2
 integer, intent(in)    :: nchild,ichildren(nchild),iparent
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 integer :: i,ilucky,seed

 if (rand_type==0) then
    ! lucky last child becomes parent
    xyzh(4,iparent) = xyzh(4,ichildren(nchild)) * (nchild)**(1./3.)
 elseif (rand_type==1) then
    ! randomly pick a child to become the parent
    ilucky = int(ran2(seed)*(nchild))+1
    call copy_particle_all(ichildren(ilucky),iparent,.true.)
    xyzh(4,iparent) = xyzh(4,iparent) * (nchild)**(1./3.)
 else
    print*, 'invalid random type',rand_type
    stop
 endif
 call set_particle_type(iparent,igas)

 ! kill the useless children
 do i=1,nchild-1
    call kill_particle(ichildren(i),npartoftype(:))
 enddo

 !-- tidy up particle types
 npartoftype(igas) = npartoftype(igas) + 1
 npartoftype(isplit) = npartoftype(isplit) - 1

end subroutine fast_merge_into_a_particle

!-----------------------------------------------------------------------
!+
! make a new ghost particle from an existing particle
!+
!-----------------------------------------------------------------------
subroutine make_a_ghost(iighost,ireal,npartoftype,npart,nchild_in,xyzh)
  use part, only:copy_particle_all,set_particle_type,ighost
  integer, intent(in)    :: iighost,nchild_in,ireal
  integer, intent(inout) :: npartoftype(:),npart
  real, intent(inout)    :: xyzh(:,:)

  call copy_particle_all(ireal,iighost,.true.)
  npartoftype(ighost) = npartoftype(ighost) + 1
  call set_particle_type(iighost,ighost)
  ! adjust smoothing length
  xyzh(4,iighost) = xyzh(4,iighost) * (nchild_in)**(1./3.)

end subroutine make_a_ghost

!-----------------------------------------------------------------------
!+
! makes split ghost particles from an existing particle
!+
!-----------------------------------------------------------------------
subroutine make_split_ghost(iighost,ireal,npartoftype,npart,nchild,xyzh,vxyzu)
  use part, only:copy_particle_all,set_particle_type,isplitghost,igas,isplit
  integer, intent(in)    :: iighost,nchild,ireal
  integer, intent(inout) :: npartoftype(:),npart
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer :: jj

  call copy_particle_all(ireal,iighost,.true.)
  call split_a_particle(nchild,iighost,xyzh,vxyzu,npartoftype,0,1,iighost)
  do jj = 0,nchild+1
    call set_particle_type(iighost+jj,isplitghost)
  enddo
  npartoftype(isplitghost) = npartoftype(isplitghost) + nchild + 1

  ! because split_a_particle is no longer generic
  npartoftype(igas) = npartoftype(igas) + 1
  npartoftype(isplit) = npartoftype(isplit) - nchild - 1

end subroutine make_split_ghost

!-----------------------------------------------------------------------
!+
! enforce periodicity
! This really should be a routine elsewhere, but it's not
!+
!-----------------------------------------------------------------------
subroutine shift_for_periodicity(npart,xyzh)
 use boundary, only:cross_boundary
 use domain,   only:isperiodic
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 integer                :: i,ncross

 ncross = 0
 !$omp parallel do schedule(static) default(none) &
 !$omp shared(npart,isperiodic,xyzh) &
 !$omp private(i) &
 !$omp reduction(+:ncross)
 do i = 1,npart
    call cross_boundary(isperiodic,xyzh(:,i),ncross)
 enddo
 !$omp end parallel do

end subroutine shift_for_periodicity

!-----------------------------------------------------------------------
!+
! quick and dirty routine to calculate the density from positions
!+
!-----------------------------------------------------------------------
subroutine quick_rho(npart,xyzh_in,rho_ave)
use kernel, only:get_kernel,cnormk,radkern
use part,   only:iamsplit,igas,isplit,massoftype,iphase,periodic
use boundary, only: dxbound,dybound,dzbound
integer, intent(in) :: npart
real,    intent(in) :: xyzh_in(:,:)
real,    intent(out) :: rho_ave
real :: mass,h1,h31,rho_i,rij_vec(3),qij,wchild,rij,grkernchild
integer :: i,j

  rho_ave = 0.
  do i=1,npart
     h1 = 1./xyzh_in(4,i)
     h31 = h1**3
     rho_i = 0.

     do j = 1,npart
       if (iamsplit(iphase(j))) then
         mass = massoftype(isplit)
       else
         mass = massoftype(igas)
       endif

       rij_vec = xyzh_in(1:3,i) - xyzh_in(1:3,j)
       if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
          if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
          if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
          if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
       endif

       rij = sqrt(dot_product(rij_vec,rij_vec))
       qij = rij*h1

       wchild = 0.
       if (qij < radkern) call get_kernel(qij*qij,qij,wchild,grkernchild)

       rho_i = rho_i + (mass*wchild*cnormk*h31)
     enddo
     rho_ave = rho_ave + rho_i
  enddo

  rho_ave = rho_ave/real(npart)

end subroutine quick_rho


!-----------------------------------------------------------------------
end module splitmergeutils
