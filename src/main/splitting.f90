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
! :Dependencies: splitutils,options
!
 use options, only:isplitboundary
 use splitmergeutils
 use part, only:igas,ighost,isplit,isplitghost,iamghost,iamsplit,iphase

 implicit none

 public :: init_split,write_options_splitting,read_options_splitting
 public :: check_split_or_merge,check_ghost_or_splitghost
 public :: update_splitting,merge_particles_wrapper

 private

 integer, parameter, public  :: nchild = 12
 integer, parameter :: nchild_in = nchild + 1
 integer, save      :: children_list(nchild_in) = 0
 ! initialise to 0 to prevent errors when writing to .in file
 real, public       :: splitrad = 0.
 real, public       :: splitbox(4) = 0.
 real, public       :: gzw = 0.

contains

!----------------------------------------------------------------
!+
!  splits particles in the initial conditions
!+
!----------------------------------------------------------------
subroutine init_split(ierr)
 use part, only:xyzh,vxyzu,massoftype,set_particle_type,npartoftype,iamtype,&
                massoftype,npart,copy_particle_all,kill_particle,shuffle_part,isdead_or_accreted
 use splitpart, only: merge_all_particles
 use dim, only:maxp_hard
 integer, intent(inout) :: ierr
 integer :: ii,merge_count,merge_ghost_count,add_npart
 integer(kind=1) :: iphaseii
 logical :: need_to_relax

ierr = 0

! particles must be checked if they are to be split, merged, ghosted or some combination
! 1. delete an existing ghost (in case restarting)
! 2. should it be spilt or merged?
! 3. should ghosts be created?
! (2 and 3 are *separate* steps)

! set masses
if (centre_particle) then
   massoftype(isplit) = massoftype(igas)/nchild_in
else
   massoftype(isplit) = massoftype(igas)/(nchild)
endif
massoftype(ighost) = massoftype(igas)
massoftype(isplitghost) = massoftype(isplit)

! 1. delete existing ghosts
call delete_all_ghosts(xyzh,vxyzu,npartoftype,npart)

! 2. should it be split or merged?

 if (rand_type==5) then
    call update_splitting(npart,xyzh,vxyzu,npartoftype,need_to_relax)
 else

merge_count = 0
add_npart = 0
do ii = 1,npart
  iphaseii = iphase(ii)
  call check_split_or_merge(ii,iphaseii,xyzh,vxyzu,npart,npartoftype,merge_count,add_npart)
enddo
npart = npart + add_npart
if (rand_type==3) then
   call merge_particles_wrapper(npart,xyzh,npartoftype,.false.)
endif
call shuffle_part(npart)

! 3. should a ghost be made?
merge_ghost_count = 0
add_npart = 0
do ii = 1,npart
  call check_ghost_or_splitghost(ii,iphase(ii),xyzh,vxyzu,npart,npartoftype,merge_ghost_count,add_npart)
enddo
npart = npart + add_npart
if (rand_type==3) then
   call merge_particles_wrapper(npart,xyzh,npartoftype,.true.)
endif
if (rand_type==4) then
   call merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu,nchild_in,npart,.true.)
endif

 endif ! end if rand_type==5

print*,'particle splitting is on, you now have',npart,'total particles'
print*,'using randomisation type ',rand_type

end subroutine init_split

!----------------------------------------------------------------
!+
!  Wrapper routine that will command all the splitting & merging
!  and making of ghosts
!+
!----------------------------------------------------------------
subroutine update_splitting(npart,xyzh,vxyzu,npartoftype,need_to_relax)
 use io,   only: fatal
 use part, only:kill_particle,isdead_or_accreted,shuffle_part,iactive,periodic
 use timestep, only:time
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, intent(inout) :: need_to_relax
 integer                :: i,k,add_npart,oldsplits,oldreals
 logical                :: split_it,ghost_it,already_split,rln
 logical                :: make_ghost,send_to_list
 integer, allocatable, dimension(:)   :: ioriginal(:)
 real,    allocatable, dimension(:,:) :: xyzh_split(:,:)
 real :: orbit

 rln = .true.
 need_to_relax = .false.
 oldsplits = npartoftype(isplit)
 oldreals = npartoftype(igas)

 if (rln) then
  ! Hack for testing
  gzw = 0.
  ! Override the splitbox size
  splitbox(:) = 0.

  ! When to do a split or a merge? Can use for either box or disc tests
  orbit = 0.2 !0.25*2.*3.1415926535*5.**1.5

  if (time > 1.*orbit .and. time < 2.*orbit .and. npartoftype(isplit)==0) call split_it_all(npart,xyzh,vxyzu)

  if (time > 2.*orbit .and. time < 3.*orbit .and. npartoftype(igas)==0) call merge_it_all(npart,xyzh,vxyzu)

  if (time > 3.*orbit .and. time < 4.*orbit .and. npartoftype(isplit)==0) call split_it_all(npart,xyzh,vxyzu)

  if (time > 4.*orbit .and. time < 5.*orbit .and. npartoftype(igas)==0) call merge_it_all(npart,xyzh,vxyzu)

  if (time > 5.*orbit .and. time < 6.*orbit .and. npartoftype(isplit)==0) call split_it_all(npart,xyzh,vxyzu)

  if (time > 6.*orbit .and. time < 7.*orbit .and. npartoftype(igas)==0) call merge_it_all(npart,xyzh,vxyzu)

  if (time > 7.*orbit .and. time < 8.*orbit .and. npartoftype(isplit)==0) call split_it_all(npart,xyzh,vxyzu)

  if (time > 8.*orbit .and. time < 9.*orbit .and. npartoftype(igas)==0) call merge_it_all(npart,xyzh,vxyzu)

  if (time > 9.*orbit .and. npartoftype(isplit)==0) call split_it_all(npart,xyzh,vxyzu)

  return !for testing
endif

 if (rand_type /= 5) then
    call fatal('update_splitting','this routine is set up only for rand_type==5')
    stop
 endif
 !
 !--Kill all active ghost particle  &
 !  Convert gas -> splits that have crossed the boundary
 !
 add_npart = 0
 split_em: do i = 1,npart
    if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       if (iamghost(iphase(i))) then
          call kill_particle(i,npartoftype)
          cycle split_em
       endif
       call inside_boundary(xyzh(1:3,i),split_it)
       already_split = iamsplit(iphase(i))

       ! is it a big particle that needs to be split?
       if (split_it .and. .not.already_split) then
          call split_a_particle(nchild,i,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)  ! nchild = 12
          add_npart = add_npart + nchild
       endif
   endif
 enddo split_em
 npart = npart + add_npart
 call shuffle_part(npart)

 ! Now go through and collect anything that should go into set_linklist call
 allocate(xyzh_split(4,npart*4))
 allocate(ioriginal(npart*4))
 k = 0
 add_npart = 0
 make_ghost = .false.
 merge_em: do i = 1,npart
   if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
      send_to_list = .false.
     ! Does the particle need a ghost?
     call inside_ghost_zone(xyzh(1:3,i),ghost_it)
     call inside_boundary(xyzh(1:3,i),split_it)
     already_split = iamsplit(iphase(i))

     ! if it is a merged particle that needs splitghosts, just make these
     if (ghost_it .and. .not.split_it .and. .not.already_split) then
       call make_split_ghost(npart+add_npart+1,i,npartoftype,npart,nchild,xyzh,vxyzu)
       add_npart = add_npart + nchild + 1
     endif

     ! if it is a split particle that needs a merged ghost, send to list
    !if (ghost_it .and. split_it .and. already_split) then
    !   send_to_list = .true.
    ! endif

    ! if it a split that has crossed into the non-split region, send to list
    !if (.not.split_it .and. already_split) then
    !   send_to_list = .true.
    ! endif
    ! for now, just send all splits
    if (already_split) send_to_list = .true.

     ! if particle needs to go to set_linklist call, add it to the list
     if (send_to_list) then
       k = k + 1
       xyzh_split(:,k) = xyzh(:,i)
       ioriginal(k) = i
     endif
   endif
 enddo merge_em

 npart = npart + add_npart

 print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'
 if (k > nchild_in) call merge_particles(npart,k,xyzh,xyzh_split(:,1:k),ioriginal,npartoftype,make_ghost)
 call shuffle_part(npart)

 deallocate(xyzh_split)
 deallocate(ioriginal)

 ! Do we need to shuffle the particles later in evolve? If you never want to, comment out below
 if (npartoftype(igas) /= oldreals .or. npartoftype(isplit) /= oldsplits) need_to_relax = .true.

 if (periodic) call shift_for_periodicity(npart,xyzh)

end subroutine update_splitting

!----------------------------------------------------------------
!+
! Testing routine that splits everything and minimises global error
!+
!----------------------------------------------------------------
subroutine split_it_all(npart,xyzh,vxyzu)
  use icosahedron, only:pixel2vector,compute_corners,compute_matrices
  use random,      only:ran2
  use vectorutils, only:rotatevec
  use physcon,     only:pi
  use splitmergeutils, only:shift_particles_WVT
  use part,        only:set_particle_type,igas,npartoftype,isplit
  use part,        only:iactive,iphase,isdead_or_accreted
  use timestep,    only:dt
  real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
  integer, intent(inout) :: npart
  real :: parents(4,npart)
  real, allocatable, dimension(:,:) :: kids_ref,shifts
  integer :: ii
  integer :: nshifts,nkids,test_ind,add_npart
  real :: err_a,err_b,err_c,err_d,error_p
  real :: mu,test(4),mu_test(4),tol
  logical :: minimum_found
  character(len=40) :: shuffle_type,error_metric_type

  ! Things to choose
  shuffle_type = 'WVT'                ! Vacondio or WVT
  error_metric_type = 'Reyes-Lopez'  ! Reyes-Lopez or FeldmanBonet

  ! Initialise
  nkids = npart*13
  add_npart = 0
  tol = 0.01 !1% between errors

  ! Store the parents separately, evaluate their error
  parents(:,1:npart) = xyzh(:,1:npart)
  if (error_metric_type == 'FeldmanBonet') then
    call calculate_parent_error(npart,parents,error_p) ! only because we're doing all of them
  else
    error_p = 0. !to prevent compiler warnings
  endif

  ! Make *proposed* splits
  allocate(kids_ref(4,nkids),shifts(3,nkids))

  do ii = 1,npart
    if (iactive(iphase(ii)) .and. .not.isdead_or_accreted(xyzh(4,ii))) then
      call split_a_particle(nchild,ii,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)  ! nchild = 12
      add_npart = add_npart + nchild
    endif
  enddo

  kids_ref = xyzh(1:4,1:nkids)

  ! Shuffle kids till error reduced
  ! Here we use the Golden Ratio method for rapid testing,
  ! but in general this should not be used
  ! Procedure is:
  ! 1. Shuffle particles somehow
  ! 2. Evaluate their error somehow
  ! 3. Repeat this till you get a minimum

  minimum_found = .false.

  ! Set bounds for mu or beta (always call it mu for convenience)
  if (shuffle_type == 'WVT') then
    nshifts = 10
    mu_test(1) = 0.0
    call WVT_loop(nshifts,nkids,xyzh,xyzh(4,1:nkids),mu_test(1))
  else
    mu_test(1) = 0.
    call shift_particles_Vacondio(nkids,xyzh(1:4,1:nkids),vxyzu(:,1:nkids),dt,mu_test(1),shifts)
    xyzh(1:3,1:nkids) = kids_ref(1:3,1:nkids) - shifts
  endif

  if (error_metric_type == 'Reyes-Lopez') then
    call calculate_RL_error(npart,parents,nkids,xyzh(1:4,1:nkids),err_a)
  else
    call calculate_global_error(nkids,npart,xyzh,parents,err_a)
    err_a = error_p - err_a
  endif

  ! Reset for next step
  xyzh(1:4,1:nkids) = kids_ref

  if (shuffle_type == 'WVT') then
    mu_test(4) = 0.13 ! this comes from Daniel's PDF, should be limited by Courant condition carefully
    call WVT_loop(nshifts,nkids,xyzh,xyzh(4,1:nkids),mu_test(4))
  else
    mu_test(4) = 10.
    call shift_particles_Vacondio(nkids,xyzh(1:4,1:nkids),vxyzu(:,1:nkids),dt,mu_test(4),shifts)
    xyzh(1:3,1:nkids) = kids_ref(1:3,1:nkids) - shifts
  endif

  if (error_metric_type == 'Reyes-Lopez') then
    call calculate_RL_error(npart,parents,nkids,xyzh(1:4,1:nkids),err_d)
  else
    call calculate_global_error(nkids,npart,xyzh,parents,err_d)
    err_d = error_p - err_d
  endif

  ! Reset for next step
  xyzh(1:4,1:nkids) = kids_ref

  print*,'  mu/beta                   error'
  print*,mu_test(1),err_a
  print*,mu_test(4),err_d

  ! if the error is always zero (this should *not* happen) let me know
  if (err_a < tiny(err_a) .and. err_d < tiny(err_d)) then
    print*,'Something wrong in minimum finding, the error is always zero'
    stop
  endif

   do while (.not.minimum_found)
     mu_test(2) = mu_test(1) + ((mu_test(4) - mu_test(1))/(1.+1.618))
     mu_test(3) = mu_test(1) + mu_test(4) - mu_test(2)

     ! Do each one
     if (shuffle_type == 'WVT') then
       call WVT_loop(nshifts,nkids,xyzh,xyzh(4,1:nkids),mu_test(2))
     else
       call shift_particles_Vacondio(nkids,xyzh(1:4,1:nkids),vxyzu(:,1:nkids),dt,mu_test(2),shifts)
       xyzh(1:3,1:nkids) = kids_ref(1:3,1:nkids) - shifts
     endif

     if (error_metric_type == 'Reyes-Lopez') then
       call calculate_RL_error(npart,parents,nkids,xyzh(1:4,1:nkids),err_b)
     else
       call calculate_global_error(nkids,npart,xyzh,parents,err_b)
       err_b = error_p - err_b
     endif

     ! Reset for next
     xyzh(1:4,1:nkids) = kids_ref

     if (shuffle_type == 'WVT') then
       call WVT_loop(nshifts,nkids,xyzh,xyzh(4,1:nkids),mu_test(3))
     else
       call shift_particles_Vacondio(nkids,xyzh(1:4,1:nkids),vxyzu(:,1:nkids),dt,mu_test(3),shifts)
       xyzh(1:3,1:nkids) = kids_ref(1:3,1:nkids) - shifts
     endif

     if (error_metric_type == 'Reyes-Lopez') then
       call calculate_RL_error(npart,parents,nkids,xyzh(1:4,1:nkids),err_c)
     else
       call calculate_global_error(nkids,npart,xyzh,parents,err_c)
       err_c = error_p - err_c
     endif

     ! Reset for next
     xyzh(1:4,1:nkids) = kids_ref

     ! What we're actually trying to minimise
     test(1) = err_a
     test(2) = err_b
     test(3) = err_c
     test(4) = err_d
     print*,mu_test(2),err_b
     print*,mu_test(3),err_c

     if (abs((test(2) - test(3))/test(3)) < tol) then
       test_ind = minloc(test,dim=1)
       mu = mu_test(test_ind)
       minimum_found = .true.
     endif

     ! if not converged, save for next round
     if (test(1) < test(3) .or. test(2)< test(3)) then
       mu_test(4) = mu_test(3)
       err_d = err_c
     else
       mu_test(1) = mu_test(2)
       err_a = err_b
     endif
   enddo
   print*,'I think I am done, here is the new mu value:',mu,'and error',test(test_ind)

   ! Override parents with the new, shuffled children
   if (shuffle_type == 'WVT') then
     call WVT_loop(nshifts,nkids,xyzh,xyzh(4,1:nkids),mu)
   else
     call shift_particles_Vacondio(nkids,xyzh(1:4,1:nkids),vxyzu(:,1:nkids),dt,mu,shifts)
     xyzh(1:3,1:nkids) = kids_ref(1:3,1:nkids) - shifts
   endif
   npart = npart + add_npart

   ! Tidy up
   deallocate(kids_ref,shifts)

end subroutine split_it_all

!----------------------------------------------------------------
!+
! Testing routine that splits everything and minimises global error
!+
!----------------------------------------------------------------
subroutine merge_it_all(npart,xyzh,vxyzu)
  use part, only:shuffle_part,npartoftype,igas,isplit,iactive,iphase,isdead_or_accreted
  use timestep, only:dt
  integer, intent(inout) :: npart
  real, intent(inout)    :: xyzh(:,:), vxyzu(:,:)
  real, allocatable, dimension(:,:) :: xyzh_split,parents_ref,shifts
  integer, allocatable, dimension(:) :: ioriginal
  integer :: k,i,nkids,nshifts,test_ind
  logical :: make_ghost,minimum_found
  real :: err_a,err_b,err_c,err_d
  real :: err_parent_a,err_parent_d,err_parent_b,err_parent_c
  real :: mu,test(4),mu_test(4),tol
  character(len=40) :: shuffle_type,error_metric_type

  ! Things to choose
  shuffle_type = 'WVT'                ! Vacondio or WVT
  error_metric_type = 'FeldmanBonet'  ! Reyes-Lopez or FeldmanBonet

  make_ghost = .true. ! this is currently overridden for testing
  tol = 0.01 !1% difference in error

  ! Store the splits as a reference
  allocate(xyzh_split(4,npart),ioriginal(npart))
  xyzh_split(:,1:npart) = xyzh(1:4,1:npart)

  ! Find particles to merge and merge them - for now, we simply send every particle for simplicity
  ! but this is not appropriate in any other circumstance
  k = 0
  do i = 1,npart
    if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
      k = k + 1
      xyzh_split(:,k) = xyzh(1:4,i)
      ioriginal(k) = i
    endif
  enddo
  nkids = k

  print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'
  if (k > nchild_in) call merge_particles(npart,k,xyzh,xyzh_split(:,1:k),ioriginal,npartoftype,make_ghost)
  call shuffle_part(npart)

  ! Shuffle till error minimised
  ! Calculate the error of each with respect to the new particle positions
  allocate(parents_ref(1:4,1:npart),shifts(3,npart))
  parents_ref = xyzh(1:4,1:npart)

  ! Shuffle parents (keeping kids still) till error reduced
  ! Here we use the Golden Ratio method for rapid testing,
  ! but in general this should not be used

  minimum_found = .false.

  if (shuffle_type == 'WVT') then
    nshifts = 10
    mu_test(1) = 0.0
    call WVT_loop(nshifts,npart,xyzh,xyzh(4,1:npart),mu_test(1))
  else
    mu_test(1) = 0.
    call shift_particles_Vacondio(npart,xyzh(1:4,1:npart),vxyzu(:,1:npart),dt,mu_test(1),shifts)
    xyzh(1:3,1:npart) = parents_ref(1:3,1:npart) - shifts
  endif

  if (error_metric_type == 'Reyes-Lopez') then
    call calculate_RL_error(npart,xyzh,nkids,xyzh_split(:,1:nkids),err_a)
  else
    call calculate_parent_error(npart,xyzh(:,1:npart),err_parent_a)
    call calculate_global_error(nkids,npart,xyzh_split(:,1:nkids),xyzh,err_a)
    err_a = err_parent_a - err_a
  endif

  ! Reset for next
  xyzh(1:4,1:npart) = parents_ref

  if (shuffle_type == 'WVT') then
    nshifts = 10
    mu_test(4) = 0.13 ! this comes from Daniel's PDF, should be limited by Courant condition carefully
    call WVT_loop(nshifts,npart,xyzh,xyzh(4,1:nkids),mu_test(4))
  else
    mu_test(4) = 10.
    call shift_particles_Vacondio(npart,xyzh(1:4,1:npart),vxyzu(:,1:npart),dt,mu_test(4),shifts)
    xyzh(1:3,1:npart) = parents_ref(1:3,1:npart) - shifts
  endif

  if (error_metric_type == 'Reyes-Lopez') then
    call calculate_RL_error(npart,xyzh,nkids,xyzh_split(:,1:nkids),err_d)
  else
    call calculate_parent_error(npart,xyzh(:,1:npart),err_parent_d)
    call calculate_global_error(nkids,npart,xyzh_split(:,1:nkids),xyzh,err_d)
    err_d = err_parent_d - err_d
  endif

  ! Reset for next
  xyzh(1:4,1:npart) = parents_ref

  print*,'  mu/beta                   error'
  print*,mu_test(1),err_a
  print*,mu_test(4),err_d

  ! if the error is always zero (this should *not* happen) let me know
  if (err_a < tiny(err_a) .and. err_d < tiny(err_d)) then
    print*,'Something wrong in minimum finding, the error is always zero'
    stop
  endif

   do while (.not.minimum_found)
     mu_test(2) = mu_test(1) + ((mu_test(4) - mu_test(1))/(1.+1.618))
     mu_test(3) = mu_test(1) + mu_test(4) - mu_test(2)

     ! Do each one
     if (shuffle_type == 'WVT') then
       call WVT_loop(nshifts,npart,xyzh,xyzh(4,1:nkids),mu_test(2))
     else
       call shift_particles_Vacondio(npart,xyzh(1:4,1:npart),vxyzu(:,1:npart),dt,mu_test(2),shifts)
       xyzh(1:3,1:npart) = parents_ref(1:3,1:npart) - shifts
     endif

     if (error_metric_type == 'Reyes-Lopez') then
       call calculate_RL_error(npart,xyzh,nkids,xyzh_split(:,1:nkids),err_b)
     else
       call calculate_parent_error(npart,xyzh(:,1:npart),err_parent_b)
       call calculate_global_error(nkids,npart,xyzh_split(:,1:nkids),xyzh,err_b)
       err_b = err_parent_b - err_b
     endif

     ! Reset for next
     xyzh(1:4,1:npart) = parents_ref

     if (shuffle_type == 'WVT') then
       call WVT_loop(nshifts,npart,xyzh,xyzh(4,1:nkids),mu_test(3))
     else
       call shift_particles_Vacondio(npart,xyzh(1:4,1:npart),vxyzu(:,1:npart),dt,mu_test(3),shifts)
       xyzh(1:3,1:npart) = parents_ref(1:3,1:npart) - shifts
     endif

     if (error_metric_type == 'Reyes-Lopez') then
       call calculate_RL_error(npart,xyzh,nkids,xyzh_split(:,1:nkids),err_c)
     else
       call calculate_parent_error(npart,xyzh(:,1:npart),err_parent_c)
       call calculate_global_error(nkids,npart,xyzh_split(:,1:nkids),xyzh,err_c)
       err_c = err_parent_c - err_c
     endif

     ! Reset for next
     xyzh(1:4,1:npart) = parents_ref

     ! What we're actually trying to minimise
     test(1) = err_a
     test(2) = err_b
     test(3) = err_c
     test(4) = err_d
     print*,mu_test(2),test(2)
     print*,mu_test(3),test(3)

     if (abs((test(2) - test(3))/test(3)) < tol) then
       test_ind = minloc(test,dim=1)
       mu = mu_test(test_ind)
       minimum_found = .true.
     endif

     ! if not converged, save for next round
     if (test(1) < test(3) .or. test(2)< test(3)) then
       mu_test(4) = mu_test(3)
       err_d = err_c
       err_parent_d = err_parent_c
     else
       mu_test(1) = mu_test(2)
       err_a = err_b
       err_parent_a = err_parent_b
     endif
   enddo
   print*,'I think I am done, here is the new mu value:',mu,'and error',test(test_ind)

  ! Keep final answer
  if (shuffle_type == 'WVT') then
    call WVT_loop(nshifts,npart,xyzh,xyzh(4,1:npart),mu)
  else
    call shift_particles_Vacondio(npart,xyzh(1:4,1:npart),vxyzu(:,1:npart),dt,mu,shifts)
    xyzh(1:3,1:npart) = parents_ref(1:3,1:npart) - shifts
  endif

  ! Tidy up
  deallocate(xyzh_split,ioriginal,parents_ref,shifts)

end subroutine merge_it_all

!----------------------------------------------------------------
!+
! Testing: WVT loop (so we don't copy this 4 times)
!+
!----------------------------------------------------------------
 subroutine WVT_loop(nshifts,n_in,xyzh_in,h0,mu)
   use splitmergeutils, only:shift_particles_WVT
   real, intent(inout) :: xyzh_in(:,:),h0(:)
   real, intent(in)    :: mu
   integer, intent(in) :: n_in,nshifts
   real :: dmu,mu_here
   integer :: ii

   dmu = real(mu/(nshifts+1))
   mu_here = mu

   ! Do the shifts loop
     do ii = 1,nshifts

       ! calculate the shifts according to the WVT method
       call shift_particles_WVT(n_in,xyzh_in,h0,mu_here)

       ! decrease mu for next go around
       mu_here = mu_here - dmu
       !print*,'ran through shifts. mu = ',mu,' for iteration i = ',ii
     enddo

 end subroutine WVT_loop


!----------------------------------------------------------------
!+
! Testing: evaluate error of kids
!+
!----------------------------------------------------------------
 subroutine calculate_global_error(nkids,nparents,kids,parents,error)
   use boundary, only: dxbound,dybound,dzbound
   use part,     only: periodic,massoftype,isplit
   use kernel, only:get_kernel,cnormk,radkern,radkern2
   real, intent(in)    :: kids(:,:),parents(:,:)
   integer, intent(in) :: nparents,nkids
   real, intent(out)   :: error
   integer :: ii,jj
   real    :: x_p(3),rij(3),rij2,wchild,grkernchild,h1,h31,mchild,qij2

   error = 0.
   mchild = massoftype(isplit)

   do ii = 1,nparents
     ! reference parent location
     x_p = parents(1:3,ii)

     ! now check across all children
     do jj = 1,nkids
       h1 = 1./kids(4,jj)
       h31 = h1**3

       rij = x_p - kids(1:3,jj)
       if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
          if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
          if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
          if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
       endif

       rij2 = dot_product(rij,rij)
       qij2 = rij2*h1*h1

       wchild = 0.
       if (qij2 < radkern2) call get_kernel(qij2,sqrt(qij2),wchild,grkernchild)

       ! Add the error
       error = error + (mchild*wchild*cnormk*h31)

     enddo
   enddo

 end subroutine calculate_global_error

 !----------------------------------------------------------------
 !+
 ! Testing: evaluate error of parents
 !+
 !----------------------------------------------------------------
  subroutine calculate_parent_error(npart,parents,error)
    use boundary, only: dxbound,dybound,dzbound
    use part,     only: periodic,massoftype,igas
    use kernel, only:get_kernel,cnormk,radkern,radkern2
    real, intent(in)    :: parents(:,:)
    integer, intent(in) :: npart
    real, intent(out)   :: error
    integer :: ii,jj
    real    :: x_p(3),rij(3),rij2,wchild,grkernchild,h1,h31,mchild,qij2

    error = 0.
    mchild = massoftype(igas)

    do ii = 1,npart
      ! reference parent location
      x_p = parents(1:3,ii)

      ! now check across all parents
      do jj = 1,npart
        h1 = 1./parents(4,jj)
        h31 = h1**3

        rij = x_p - parents(1:3,jj)
        if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
           if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
           if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
           if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
        endif

        rij2 = dot_product(rij,rij)
        qij2 = rij2*h1*h1

        wchild = 0.
        if (qij2 < radkern2) call get_kernel(qij2,sqrt(qij2),wchild,grkernchild)

        ! Add the error
        error = error + (mchild*wchild*cnormk*h31)

      enddo
    enddo

  end subroutine calculate_parent_error

!----------------------------------------------------------------
!+
!  Calculates error in the derivative according to Reyes López et al. 2012
!+
!----------------------------------------------------------------

subroutine calculate_RL_error(nparents,parents,nkids,kids,error)
  use kernel, only:get_kernel,cnormk,radkern2
  use boundary, only: dxbound,dybound,dzbound
  use part,     only: periodic
  integer, intent(in) :: nparents,nkids
  real, intent(in)    :: parents(:,:),kids(:,:)
  real, intent(out)   :: error
  integer :: ii,jj,kk
  real    :: wkid,grkernkid,wparent,grkernparent,error_kids
  real    :: h1,rij_vec(3),rij2,qij2

  ! There are about 10 trillion ways to do this faster, let's do the slow way
  ! for testing

  error = 0.

  ! Use location of parents to evaluate error
  do ii = 1,1

    ! For each parent
    do jj = 1,nparents

      h1 = 1./parents(4,jj)
      rij_vec = parents(1:3,jj) - parents(1:3,ii)

      if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
         if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
         if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
         if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
      endif

      rij2 = dot_product(rij_vec,rij_vec)
      qij2 = rij2*h1*h1

      call get_kernel(qij2,sqrt(qij2),wparent,grkernparent)
      grkernparent = grkernparent*cnormk

      ! Error introduced from the kids (across the whole domain, not just in the family)
      error_kids  = 0.
      do kk = 1,nkids

        h1 = 1./kids(4,kk)
        rij_vec = kids(1:3,kk) - parents(1:3,ii)

        if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
           if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
           if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
           if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
        endif

        rij2 = dot_product(rij_vec,rij_vec)
        qij2 = rij2*h1*h1

        grkernkid = 0.
        if (qij2 < radkern2) call get_kernel(qij2,sqrt(qij2),wkid,grkernkid)

        error_kids = error_kids + grkernkid*cnormk/13.

      enddo

      error = error + (grkernparent - error_kids)**2
    enddo
  enddo

end subroutine calculate_RL_error

!----------------------------------------------------------------
!+
! Testing routine that splits everything and minimises global error
! according to Daniel's research notes (his Eq 17)
!+
!----------------------------------------------------------------
subroutine split_it_all_Price(npart,xyzh,vxyzu)
  use kernel, only:grkern,radkern2,get_kernel,cnormk
  use boundary, only: dxbound,dybound,dzbound
  use part,     only: periodic,massoftype,igas,isplit,iactive,iphase
  use part,     only: npartoftype,isdead_or_accreted,rhoh,gradh
  integer, intent(inout) :: npart
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  real, allocatable      :: shifts(:,:),omega_kids(:),sum_parents(:,:)
  real     :: parents(1:4,1:npart),gradh_parents(npart),mkid,mparent,rij_vec(3),rij2
  real     :: h1a,h1b,qij2a,qij2b,grkerna,grkernb,sum_kids(3)
  real     :: grkernap,grkernbp,omegaa,omegab,omegaap,omegabp,rhoa,rhob,rhoap,rhobp
  real     :: x_a(3),rij_unit(3),maxshift,qij2,sum
  integer  :: nkids,nparents,add_npart,ii,jj,kk


  ! Initialise
  nkids = npart*13
  nparents = npart
  add_npart = 0
  mkid = massoftype(isplit)
  mparent = massoftype(igas)

  ! Store the parents separately for later calculations
  parents(1:4,1:npart) = xyzh(1:4,1:npart)
  gradh_parents(1:npart) = gradh(1,1:npart)

  ! Make *proposed* splits
  allocate(shifts(3,nkids),omega_kids(nkids),sum_parents(3,nkids))

  sum_parents = 0.
  do ii = 1,npart
    if (iactive(iphase(ii)) .and. .not.isdead_or_accreted(xyzh(4,ii))) then

      ! used for the parent sum below
      h1a = 1./parents(4,ii)
      omegaap = 1./gradh_parents(ii)
      rhoap = rhoh(xyzh(4,ii),mparent)

      ! now make the children from the parent
      call split_a_particle(nchild,ii,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)  ! nchild = 12

      ! calculate the contribution to the parents sum for each new child
      ! Here I have *incorrectly* assumed that the parental contribution for each
      ! of its children is the same. This should be changed to interpolate the
      ! parent properties for any off-centre children (probably why the resulting sum is so small!)
      do jj = 1,npart
          if (jj == ii) cycle
          h1b = 1./parents(4,jj)
          omegabp = 1./gradh_parents(jj)
          rhobp = rhoh(xyzh(4,jj),mparent)

          rij_vec = parents(1:3,jj) - xyzh(1:3,ii)

          if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
            if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
            if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
            if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
          endif
          rij2 = dot_product(rij_vec,rij_vec)
          rij_unit = rij_vec/sqrt(rij2)

          qij2a = rij2*h1a*h1a
          grkernap = 0.
          if (qij2a < radkern2) grkernap = grkern(qij2a,sqrt(qij2a))*cnormk

          qij2b = rij2*h1b*h1b
          grkernbp = 0.
          if (qij2b < radkern2) grkernbp = grkern(qij2b,sqrt(qij2b))*cnormk

          sum_parents(:,ii) = sum_parents(:,ii) + (grkernap*rij_unit/(omegaap*rhoap)) &
                                                + (grkernbp*rij_unit/(omegabp*rhobp))
      enddo

      ! save for the other associated children
      do kk = npart+add_npart,npart+add_npart+nchild
        sum_parents(:,kk) = sum_parents(:,ii)
      enddo

      ! add on the new particles
      add_npart = add_npart + nchild
    endif
  enddo

  sum_parents = mparent*sum_parents

  npart = npart + add_npart

  ! I think we need to calculate the omega terms on the proposed kids,
  ! because this has not yet been calculated on them (just inherited by parents, which will not be accurate)
  do ii = 1,nkids
    h1a = xyzh(4,ii)
    rhoa = rhoh(xyzh(4,ii),mkid)
    sum = 0.

    do jj = 1,nkids
      h1b = 1./xyzh(4,jj)
      rij_vec = xyzh(1:3,jj) - xyzh(1:3,ii)

      if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
        if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
        if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
        if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
      endif
      rij2 = dot_product(rij_vec,rij_vec)
      if (ii==jj) then
        qij2 = 0.
      else
        rij_unit = rij_vec/sqrt(rij2)
        qij2 = rij2*h1a*h1a
      endif

      grkerna = 0.
      if (qij2 < radkern2) grkerna = grkern(qij2,sqrt(qij2))*cnormk

      sum = sum + grkerna
    enddo

    omega_kids(ii) = 1. + 3.*mkid*sum/(h1a*rhoa)

  enddo

  ! for each new child work out the proposed shift
  do ii = 1,nkids
    x_a = xyzh(1:3,ii)
    h1a = xyzh(4,ii)
    omegaa = omega_kids(ii)
    rhoa = rhoh(xyzh(4,ii),mkid)

    ! Sum over kids only
    sum_kids = 0.
    do jj = 1,nkids
      h1b = 1./xyzh(4,jj)
      omegab = omega_kids(jj)
      rhob = rhoh(xyzh(4,jj),mkid)

      rij_vec = xyzh(1:3,jj) - xyzh(1:3,ii)

      if (periodic) then  ! Should make this a pre-processor statement since this can be expensive if we're not using periodic BC's
        if (abs(rij_vec(1)) > 0.5*dxbound) rij_vec(1) = rij_vec(1) - dxbound*SIGN(1.0,rij_vec(1))
        if (abs(rij_vec(2)) > 0.5*dybound) rij_vec(2) = rij_vec(2) - dybound*SIGN(1.0,rij_vec(2))
        if (abs(rij_vec(3)) > 0.5*dzbound) rij_vec(3) = rij_vec(3) - dzbound*SIGN(1.0,rij_vec(3))
      endif
      rij2 = dot_product(rij_vec,rij_vec)

      if (ii==jj) then
        qij2a = 0.
        qij2b = 0.
        rij_unit = 0.
      else
        rij_unit = rij_vec/sqrt(rij2)

        qij2a = rij2*h1a*h1a
        grkerna = 0.

        qij2b = rij2*h1b*h1b
        grkernb = 0.
      endif

      if (qij2a < radkern2) grkerna = grkern(qij2a,sqrt(qij2a))*cnormk
      if (qij2b < radkern2) grkernb = grkern(qij2b,sqrt(qij2b))*cnormk

      sum_kids = sum_kids + (grkerna*rij_unit/(omegaa*rhoa)) &
                          + (grkernb*rij_unit/(omegab*rhob))

    enddo
    sum_kids = sum_kids*mkid

    ! Combine the sums
    shifts(:,ii) = -0.5*xyzh(4,ii)**2 * (sum_kids - sum_parents(:,ii))

    !if (modulo(ii,13) == 3) print*,'shifts',shifts(:,ii)
  enddo
!  read*

  maxshift = maxval(abs(shifts))
  if (maxshift > 0.5) print*,'WARNING SHIFTS ARE LARGE',maxshift

  ! Apply the shifts
  xyzh(1:3,1:npart) = xyzh(1:3,1:npart) + shifts(:,1:npart)

! Tidy up
deallocate(shifts,omega_kids,sum_parents)

end subroutine split_it_all_Price

!----------------------------------------------------------------
!+
!  determines whether a particle should be split or merged, then
!  goes ahead and does the appropriate operation
!+
!----------------------------------------------------------------
subroutine check_split_or_merge(ii,iphaseii,xyzh,vxyzu,npart,npartoftype,merge_count,add_npart)
  use part, only:set_particle_type,iamtype,shuffle_part,kill_particle
  use random, only:ran2
  integer, intent(in)    :: ii
  integer(kind=1), intent(in) :: iphaseii
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer, intent(inout) :: npartoftype(:),npart
  integer, intent(inout) :: merge_count,add_npart
  integer :: ilucky, seed
  logical :: already_split,split_it

  call inside_boundary(xyzh(1:3,ii),split_it)
  already_split = iamsplit(iphaseii)
  ! if it should be split but is not, split it
  if (split_it .and. .not.already_split) then

    call split_a_particle(nchild,ii,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)
    add_npart = add_npart + nchild

  ! if it should not be split, but already is, merge it
  else if (.not.split_it .and. already_split .and. rand_type < 3) then
    if (rand_type==2) then
        ilucky = int(ran2(seed)*nchild_in)+1
        if (ilucky==1) then
           ! particle becomes merged
           xyzh(4,ii) = xyzh(4,ii) * nchild_in**(1./3.)
           call set_particle_type(ii,igas)
           npartoftype(igas) = npartoftype(igas) + 1
           npartoftype(isplit) = npartoftype(isplit) - 1
        else
           ! die, particle, die!
           call kill_particle(ii,npartoftype(:))
        endif
    else
       merge_count = merge_count + 1
       children_list(merge_count) = ii
       if (merge_count == nchild_in) then
         call fast_merge_into_a_particle(nchild_in,children_list,npart,xyzh,vxyzu,npartoftype,ii)
         merge_count = 0
       endif
     endif
  endif

end subroutine check_split_or_merge

!----------------------------------------------------------------
!+
!  determines whether a particle should have a ghost, then
!  goes ahead and makes it appropriately - assumes there are no
!  ghosts to start with
!+
!----------------------------------------------------------------
subroutine check_ghost_or_splitghost(ii,iphaseii,xyzh,vxyzu,npart,&
  npartoftype,merge_ghost_count,add_npart)
  use part, only:set_particle_type
  use random, only:ran2
  integer, intent(in)    :: ii
  integer(kind=1), intent(in) :: iphaseii
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer, intent(inout) :: npartoftype(:),npart
  integer, intent(inout) :: add_npart,merge_ghost_count
  integer :: ilucky,seed
  logical :: already_split,ghost_it,already_ghost

  call inside_ghost_zone(xyzh(1:3,ii),ghost_it)
  already_ghost = iamghost(iphaseii)
  if (ghost_it .and. .not.already_ghost) then
    already_split = iamsplit(iphaseii)
    ! this is inside the boundary and should be merged (just use the original parent)
    if (already_split .and. rand_type < 3) then
      if (rand_type==2) then
        ilucky = int(ran2(seed)*nchild_in)+1
        if (ilucky==1) then
           call make_a_ghost(npart+add_npart+1,ii,npartoftype,npart,nchild_in,xyzh)
           add_npart = add_npart + 1
        endif
      else
         merge_ghost_count = merge_ghost_count + 1
         children_list(merge_ghost_count) = ii
         if (merge_ghost_count == nchild_in) then
            if (rand_type==1) then
               ilucky = int(ran2(seed)*nchild_in)+1
               ilucky = children_list(ilucky)
               !print*, ilucky,children_list(ilucky)
            elseif (rand_type==0) then
               ilucky = ii
            else
               print*, 'This should not be possible'
               stop
            endif
            call make_a_ghost(npart+add_npart+1,ilucky,npartoftype,npart,nchild_in,xyzh)
            merge_ghost_count = 0
            add_npart = add_npart + 1
         endif
      endif
    ! this is outside the boundary and should be split
    else if (.not.already_split) then
      call make_split_ghost(npart+add_npart+1,ii,npartoftype,npart,nchild,xyzh,vxyzu)
      add_npart = add_npart + nchild + 1
    endif
  endif
end subroutine check_ghost_or_splitghost
!----------------------------------------------------------------
!+
!  A wrapper programme for merge_particles
!+
!----------------------------------------------------------------
subroutine merge_particles_wrapper(npart,xyzh,npartoftype,make_ghost)
 use part, only:iactive,isdead_or_accreted,set_particle_type,kill_particle
 integer, intent(inout) :: npartoftype(:),npart
 real,    intent(inout) :: xyzh(:,:)
 logical, intent(in)    :: make_ghost
 integer, allocatable, dimension(:)   :: ioriginal(:)
 real,    allocatable, dimension(:,:) :: xyzh_split(:,:)
 integer :: i,k
 logical :: ghost_it,split_it

 allocate(xyzh_split(4,npartoftype(isplit))) ! we can be smarter and choose smaller arrays later
 allocate(ioriginal(npartoftype(isplit)))        ! we can be smarter and choose smaller arrays later

 ! Determine the number of particles
 k = 0
 do i = 1,npart
    if (iphase(i)==isplit .and. iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       if (make_ghost) then
          call inside_ghost_zone(xyzh(1:3,i),ghost_it)
          if (ghost_it) then ! this should be only the split particles that need to have ghosts made of them
             k = k + 1
             xyzh_split(:,k) = xyzh(:,i)
             ioriginal(k) = i
          endif
       else
          call inside_boundary(xyzh(1:3,i),split_it)
          if (.not. split_it) then ! this should be only the split particles that need to be converted into gas particles
             k = k + 1
             xyzh_split(:,k) = xyzh(:,i)
             ioriginal(k) = i
          endif
       endif
    endif
 enddo
 print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'

 if (k > nchild_in) call merge_particles(npart,k,xyzh,xyzh_split(:,1:k),ioriginal,npartoftype,make_ghost)

 deallocate(xyzh_split)
 deallocate(ioriginal)

end subroutine merge_particles_wrapper
!----------------------------------------------------------------
!+
!  This is a subroutine to merge particles, either to create a gas
!  particle or a gas ghost particle
!  if make_ghosts = true, then make ghosts from the split particles
!  if make_ghosts = false, then merge split particles into gas particles
!+
!----------------------------------------------------------------
subroutine merge_particles(npart,ncandidate,xyzh,xyzh_split,ioriginal,npartoftype,make_ghost)
 use linklist,       only:ncells,get_neighbour_list,get_hmaxcell,get_cell_location,set_linklist,ifirstincell
 use kdtree,         only:inodeparts,inoderange
 use part,           only:iactive,isdead_or_accreted,set_particle_type,kill_particle,copy_particle_all
 use mpiforce,       only:cellforce
 use io,             only:fatal
 integer, intent(inout) :: npartoftype(:),npart
 real,    intent(inout) :: xyzh(:,:),xyzh_split(:,:)
 logical, intent(in)    :: make_ghost
 integer, intent(in)    :: ioriginal(:)
 integer, intent(inout) :: ncandidate  ! needs to be inout to satisfy if MPI, but we can't have pre-processor statements in this .f90 file
 integer                :: i,j,k,icell,jmin,n_cell,iave(4),remainder
 real                   :: r2,r2min,hmax,hcell
 type(cellforce)        :: cell
 logical                :: ghost_it,split_it

 ! Make sure only a multiple of nchild_in is ever used here (for even division of particles amongst cells)
 ! (just forget the rest, always a small number and normally gets picked up in the next time-step)
 remainder = modulo(ncandidate,nchild_in)
 ncandidate = ncandidate - remainder

 ! Use the tree to make  the link list, the second last term is not use
 call set_linklist(ncandidate,ncandidate,xyzh_split(:,1:ncandidate),xyzh(:,1:ncandidate),special_split=nchild_in)

 ! Go through all the cells, and create a ghost .or. merge them
 ! in both cases, the new particle will be based upon the centre-most particle in the cell
 iave    = 0
 iave(3) = huge(iave(3))
 hmax    = 0.
 over_cells: do icell=1,int(ncells)
    i = ifirstincell(icell)
    if (i <= 0) cycle over_cells !--skip empty cells AND inactive cells
    n_cell = inoderange(2,icell)-inoderange(1,icell)+1
    iave(1) = iave(1) + 1
    iave(2) = iave(2) + n_cell
    iave(3) = min(iave(3),n_cell)
    iave(4) = max(iave(4),n_cell)
    ! find centre-most particle
    call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
    call get_hmaxcell(icell,hcell)
    !print*, make_ghost, icell, int(ncells), inoderange(2,icell)-inoderange(1,icell)+1,hcell,cell%xpos,cell%xsizei,cell%rcuti
    hmax  = max(hmax,hcell)
    r2min = huge(r2min)
    jmin  = 0
    do j = inoderange(1,icell),inoderange(2,icell)
       k = inodeparts(j)
       r2 = (xyzh_split(1,k) - cell%xpos(1))**2 &
          + (xyzh_split(2,k) - cell%xpos(2))**2 &
          + (xyzh_split(3,k) - cell%xpos(3))**2
       if (r2 < r2min) then
          r2min = r2
          jmin  = j
       endif
    enddo

    ! for each cell generated, if the
    ! centroid is inside ghost region -> make ghost from central particle
    ! centroid is outside split region -> merge splits to make a real
    ! centroid is outside split *and* inside ghost, do both
    k = ioriginal(inodeparts(jmin))
  !  hack to take into account periodic boundaries
    if ((xyzh(1,k) > 0.5)) xyzh(1,k) = xyzh(1,k) - 1.0
    if ((xyzh(1,k) < -0.5)) xyzh(1,k) = xyzh(1,k) + 1.0
    if ((xyzh(2,k) > 0.5)) xyzh(2,k) = xyzh(2,k) - 1.0
    if ((xyzh(2,k) < -0.5)) xyzh(2,k) = xyzh(2,k) + 1.0

    call inside_ghost_zone(xyzh(1:3,k),ghost_it)
    call inside_boundary(xyzh(1:3,k),split_it)

    if (.not.ghost_it .and. .not.split_it) then
      ! this particle should be straight up merged, get rid of the rest
      do j = inoderange(1,icell),inoderange(2,icell)
         k = ioriginal(inodeparts(j))
         if (j==jmin) then
            xyzh(4,k) = xyzh(4,k) * nchild_in**(1./3.)
            call set_particle_type(k,igas)
            npartoftype(igas)   = npartoftype(igas)   + 1
            npartoftype(isplit) = npartoftype(isplit) - 1
         else
            call kill_particle(k,npartoftype(:))
         endif
      enddo
    endif

    if (ghost_it) then
      if (.not.split_it) then
        ! need to make a merged particle and turn the rest into ghosts
        do j = inoderange(1,icell),inoderange(2,icell)
          k = ioriginal(inodeparts(j))
          if (j==jmin) then
              ! copy particle to make the merged
              npart = npart + 1
              call copy_particle_all(k,npart,.true.)
              xyzh(4,npart) = xyzh(4,k) * nchild_in**(1./3.)
              call set_particle_type(npart,igas)
              npartoftype(igas) = npartoftype(igas) + 1
          endif
            ! turn all particles into splitghosts
            call set_particle_type(k,isplitghost)
            npartoftype(isplitghost) = npartoftype(isplitghost) + 1
            npartoftype(isplit) = npartoftype(isplit) - 1
        enddo
      elseif (split_it) then
        ! just need a merged ghost here
        k = ioriginal(inodeparts(jmin))
        npart = npart + 1
        call make_a_ghost(npart,k,npartoftype,npart,nchild_in,xyzh)
      endif
    endif

 enddo over_cells
 print*, 'min, average, and maximum particles per cell: ',iave(3),float(iave(2))/float(iave(1)),iave(4)

end subroutine merge_particles
!----------------------------------------------------------------
!+
!  determines whether the particle is inside boundary and thus
!  should be split
!+
!----------------------------------------------------------------
subroutine inside_boundary(pos,should_split)
  real, intent(in)     :: pos(3)
  logical, intent(out) :: should_split
  real :: rad2

  should_split = .false.
  if (isplitboundary == 0) then !circle
    rad2 = dot_product(pos,pos)
    if (rad2 < splitrad**2) should_split = .true.
  else                         !rectangle box
    if (pos(1) > splitbox(3) .and. pos(1) < splitbox(1) .and. &
        pos(2) > splitbox(4) .and. pos(2) < splitbox(2)) should_split = .true.
  endif

end subroutine inside_boundary

!----------------------------------------------------------------
!+
!  determines whether the particle is inside region where ghosts
!  should exist
!+
!----------------------------------------------------------------
subroutine inside_ghost_zone(pos,should_ghost)
  real, intent(in)     :: pos(3)
  logical, intent(out) :: should_ghost
  real :: rad2
  logical :: in_big_box,in_small_box

  should_ghost = .false.
  in_big_box = .false.
  in_small_box = .false.
  if (isplitboundary == 0) then !circle
    rad2 = dot_product(pos,pos)
    if (rad2 < (splitrad+gzw)**2 .and. &
        rad2 > (splitrad-gzw)**2) should_ghost = .true.
  else                         !rectangle box
    if (pos(1) > (splitbox(3)-gzw) .and. pos(1) < (splitbox(1)+gzw) .and. &
        pos(2) > (splitbox(4)-gzw) .and. pos(2) < (splitbox(2)+gzw)) in_big_box=.true.
    if (pos(1) > (splitbox(3)+gzw) .and. pos(1) < (splitbox(1)-gzw) .and. &
        pos(2) > (splitbox(4)+gzw) .and. pos(2) < (splitbox(2)-gzw)) in_small_box=.true.
    if (in_big_box .and. .not.in_small_box) should_ghost = .true.
  endif

end subroutine inside_ghost_zone

!-----------------------------------------------------------------------
!+
! identifies and kills all ghosts and splitghosts
!+
!-----------------------------------------------------------------------
subroutine delete_all_ghosts(xyzh,vxyzu,npartoftype,npart)
  use part, only:kill_particle,iphase,shuffle_part
  integer, intent(inout) :: npartoftype(:),npart
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer :: ii

  do ii = 1,npart
    if (iamghost(iphase(ii))) then
      call kill_particle(ii,npartoftype)
    endif
  enddo

  call shuffle_part(npart)

end subroutine delete_all_ghosts

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_splitting(iunit,isplitboundary)
 use infile_utils, only:write_inopt,get_optstring
 integer, intent(in) :: iunit,isplitboundary

 write(iunit,"(/,a)") '# options relating to particle splitting'

 !call get_optstring(iboundarytype,boundarytype,string,4)
 call write_inopt(isplitboundary,'boundary type','particles are split inside here, 0: circle, 1: rectangle',iunit)
 call write_inopt(gzw,'boundary width','width of boundary for ghosts in code units',iunit)

 select case(isplitboundary)
 case(0)
    call write_inopt(splitrad,'radius','radius of boundary in code units',iunit)

 case(1)
    call write_inopt(splitbox(3),'left edge','left  edge of boundary in code units',iunit)
    call write_inopt(splitbox(1),'right edge','right edge of boundary in code units',iunit)
    call write_inopt(splitbox(2),'upper edge','upper edge of boundary in code units',iunit)
    call write_inopt(splitbox(4),'lower edge','lower edge of boundary in code units',iunit)
 end select

end subroutine write_options_splitting

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_splitting(name,valstring,imatch,igotall,ierr,isplitboundary)
 use io, only:fatal
 character(len=*), intent(in)    :: name,valstring
 logical,          intent(out)   :: imatch,igotall
 integer,          intent(out)   :: ierr
 integer,          intent(inout) :: isplitboundary
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_splitting'
 logical       :: igotsplitbox,igotsplitrad

 imatch  = .true.
 igotall = .false.
 igotsplitbox = .false.
 igotsplitrad = .false.

 select case(trim(name))
 case('boundary type')
    read(valstring,*,iostat=ierr) isplitboundary
    ngot = ngot + 1
    if (isplitboundary  <  0 .or. isplitboundary > 2) call fatal(label,'silly choice of isplitboundary in input options')
 case('boundary width')
    read(valstring,*,iostat=ierr) gzw
    ngot = ngot + 1
    !if (gzw < tiny(gzw)) call fatal(label,'gzw must be > 0. in input options')
 case default
    imatch = .false.
    select case(isplitboundary)
    case(0)
      call read_options_splitrad(name,valstring,imatch,igotsplitrad,ierr)
      splitbox = 0.
      igotsplitbox = .true.
    case(1)
      call read_options_splitbox(name,valstring,imatch,igotsplitbox,ierr)
      splitrad = 0.
      igotsplitrad = .true.
    end select
  end select

 igotall = (ngot >= 1 .and. igotsplitrad .and. igotsplitbox)

end subroutine read_options_splitting

!-----------------------------------------------------------------------
!+
!  Reads boundary options for a circle from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_splitrad(name,valstring,imatch,igotsplitrad,ierr)
  use io,      only:fatal
  use physcon, only:pi
  character(len=*), intent(in)  :: name,valstring
  logical,          intent(out) :: imatch,igotsplitrad
  integer,          intent(out) :: ierr
  integer, save :: ngot = 0
  character(len=30), parameter :: label = 'read_options_splitrad'

 imatch  = .true.
 igotsplitrad = .false.

 select case(trim(name))
 case('radius')
    read(valstring,*,iostat=ierr) splitrad
    ngot = ngot + 1
    if (splitrad < 0.) call fatal(label,'radius must be > 0. in input options')
 case default
    imatch = .false.
  end select

 igotsplitrad = (ngot >= 1)

end subroutine read_options_splitrad

!-----------------------------------------------------------------------
!+
!  Reads boundary options for a rectangle from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_splitbox(name,valstring,imatch,igotsplitbox,ierr)
  use io,      only:fatal
  use physcon, only:pi
  character(len=*), intent(in)  :: name,valstring
  logical,          intent(out) :: imatch,igotsplitbox
  integer,          intent(out) :: ierr
  integer, save :: ngot = 0
  character(len=30), parameter :: label = 'read_options_splitbox'

 imatch  = .true.
 igotsplitbox = .false.

 select case(trim(name))
 case('left edge')
    read(valstring,*,iostat=ierr) splitbox(3)
    ngot = ngot + 1
 case('right edge')
    read(valstring,*,iostat=ierr) splitbox(1)
    ngot = ngot + 1
 case('upper edge')
    read(valstring,*,iostat=ierr) splitbox(2)
    ngot = ngot + 1
 case('lower edge')
    read(valstring,*,iostat=ierr) splitbox(4)
    ngot = ngot + 1
 case default
    imatch = .false.
  end select

 if (ngot >=4 ) then
   if (splitbox(1) < splitbox(3)) call fatal(label,'right/left boundary wrong way around in input options')
   if (splitbox(2) < splitbox(4)) call fatal(label,'upper/lower boundary wrong way around in input options')
 endif

 igotsplitbox = (ngot >= 1)

end subroutine read_options_splitbox

end module split
