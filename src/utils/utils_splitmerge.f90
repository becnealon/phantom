!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
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

contains

!--------------------------------------------------------------------------
!+
!  splits a particle into nchild particles
!+
!--------------------------------------------------------------------------
subroutine split_a_particle(nchild,iparent,xyzh,vxyzu, &
           npartoftype,lattice_type,ires,ichildren)
 use icosahedron, only:pixel2vector,compute_corners,compute_matrices
 use part,        only:copy_particle_all,igas,isplit,set_particle_type
 use part,        only:kill_particle
 integer, intent(in)    :: nchild,iparent,lattice_type,ires,ichildren
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 integer :: j,iseed,ichild
 real    :: dhfac,dx(3),sep,geodesic_R(0:19,3,3), geodesic_v(0:11,3)

 if (lattice_type == 0) then
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 else
    ! initialise random number generator
    iseed = -6542
 endif

 sep = xyzh(4,iparent)
 ichild = 0
 dhfac = 1./(nchild)**(1./3.)

 do j=0,nchild-2
    ichild = ichild + 1
    ! copy properties
    call copy_particle_all(iparent,ichildren+ichild)

    ! adjust the position
    if (lattice_type == 0) then
       call pixel2vector(j,ires,geodesic_R,geodesic_v,dx)
    else
       call sample_kernel(iseed,dx)
    endif
    xyzh(1:3,ichildren+ichild) = xyzh(1:3,iparent) + sep*dx(:)
    call set_particle_type(ichildren+ichild,isplit)
    xyzh(4,ichildren+ichild) = xyzh(4,iparent)*dhfac
 enddo

 !--fix parent
 xyzh(4,iparent) = xyzh(4,iparent)*dhfac
 call set_particle_type(iparent,isplit)

 !-- tidy up particle types
 npartoftype(igas) = npartoftype(igas) - 1
 npartoftype(isplit) = npartoftype(isplit) + nchild

 call kill_particle(iparent,npartoftype)

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
subroutine shift_particles(npart,xyzh,vxyzu,deltat,beta,shifts)
 use kernel, only:radkern2
 real, intent(in)    :: xyzh(:,:),vxyzu(:,:)
 real, intent(in)    :: deltat,beta
 integer, intent(in) :: npart
 real, intent(out)   :: shifts(3,npart)
 integer             :: i,j,neighbours
 real                :: rnaught,rij2,dr3,vel2,vmax
 real                :: q2,rij(3),rsum(3)

 vmax = tiny(vmax)
 vel2 = 0.

 do i = 1,npart
    rnaught = 0.
    neighbours = 0
    rsum = 0.
    vel2 = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    if (vel2 > vmax) vmax = vel2

    over_npart: do j = 1,npart
       if (i == j) cycle over_npart
       rij = xyzh(1:3,j) - xyzh(1:3,i)
       rij2 = dot_product(rij,rij)
       q2 = rij2/(xyzh(4,i)*xyzh(4,i))

       if (q2 < radkern2) then
          neighbours = neighbours + 1
          rnaught = rnaught + sqrt(rij2)
       endif

       dr3 = 1./(rij2**1.5)
       rsum = rsum + (rij*dr3)
    enddo over_npart

    rnaught = rnaught/neighbours
    shifts(:,i) = beta*rnaught*rnaught*deltat*rsum
 enddo

 shifts = shifts*sqrt(vel2)

end subroutine shift_particles

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
 call copy_particle_all(ichildren(1),iparent)

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
 integer, intent(in)    :: nchild,ichildren(nchild),iparent
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 integer :: i

 ! lucky last child becomes parent
 call copy_particle_all(ichildren(nchild),iparent)
 xyzh(4,iparent) = xyzh(4,ichildren(nchild)) * (nchild)**(1./3.)
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
subroutine make_a_ghost(iighost,ireal,npartoftype,npart,nchild,xyzh)
  use part, only:copy_particle_all,set_particle_type,ighost
  integer, intent(in)    :: iighost,nchild,ireal
  integer, intent(inout) :: npartoftype(:),npart
  real, intent(inout)    :: xyzh(:,:)

  call copy_particle_all(ireal,iighost)
  npartoftype(ighost) = npartoftype(ighost) + 1
  call set_particle_type(iighost,ighost)

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

  call copy_particle_all(ireal,iighost)
  call split_a_particle(nchild,iighost,xyzh,vxyzu, &
             npartoftype,0,1,iighost)
  do jj = 0,nchild
    call set_particle_type(iighost+jj,isplitghost)
  enddo
  npartoftype(isplitghost) = npartoftype(isplitghost) + nchild
  ! because split_a_particle is no longer generic
  npartoftype(igas) = npartoftype(igas) + 1
  npartoftype(isplit) = npartoftype(isplit) - nchild

end subroutine make_split_ghost

end module splitmergeutils
