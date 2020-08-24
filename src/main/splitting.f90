!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module split
!
! particle splitting module
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: splitutils
!
 implicit none

 public :: init_split

 private

 integer :: nchild = 12

contains

!----------------------------------------------------------------
!+
!  splits particles in the initial conditions
!+
!----------------------------------------------------------------
subroutine init_split(ierr)
 use part, only:xyzh,vxyzu,massoftype,set_particle_type,npartoftype,&
                massoftype,igas,isplit,npart
 use splitmergeutils, only:split_a_particle
 integer, intent(inout) :: ierr
 integer :: ii,jj,ichild
 logical :: split_it

ierr = 0

! split particles that are "inside" the high res boundary
ichild = npart
massoftype(isplit+1) = massoftype(igas)/nchild
do ii = 1,npart
  call check_should_split(xyzh(1:3,ii),split_it)
  if (split_it) then
    print*,'splitting',ii
    npart = npart + nchild - 1
    npartoftype(igas) = npartoftype(igas) - 1
    npartoftype(isplit+1) = npartoftype(isplit+1) + nchild

    call split_a_particle(nchild,ii,xyzh,vxyzu,0,1,ichild)

    ! replace the parent
    call set_particle_type(ii,isplit+1)

    ! set the children
    do jj = 1,nchild-1
      call set_particle_type(ichild+jj,isplit+1)
    enddo

    ichild = ichild + nchild - 1
  endif
enddo

npart = ichild

end subroutine init_split

!----------------------------------------------------------------
!+
!  determines whether the particle should be split or not
! TODO: make fancier, take in parameters from *.in file and
! do either rectangles or circles
!+
!----------------------------------------------------------------
subroutine check_should_split(pos,should_split)
  real, intent(in)     :: pos(3)
  logical, intent(out) :: should_split

  should_split = .false.
  if (pos(2) > 0.5) should_split = .true.

end subroutine check_should_split

end module split
