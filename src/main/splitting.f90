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
 public :: inside_boundary,check_split_or_merge,check_ghost_or_splitghost

 private

 integer, parameter, public  :: nchild = 12
 integer, save    :: children_list(nchild) = 0
 real, public     :: splitrad,splitbox(4),gzw

contains

!----------------------------------------------------------------
!+
!  splits particles in the initial conditions
!+
!----------------------------------------------------------------
subroutine init_split(ierr)
 use part, only:xyzh,vxyzu,massoftype,set_particle_type,npartoftype,iamtype,&
                massoftype,npart,copy_particle_all,kill_particle,shuffle_part,isdead_or_accreted
 use dim, only:maxp_hard
 integer, intent(inout) :: ierr
 integer :: ii,merge_count,merge_ghost_count,add_npart
 integer(kind=1) :: iphaseii

ierr = 0

! particles must be checked if they are to be split, merged, ghosted or some combination
! 1. delete an existing ghost (in case restarting)
! 2. should it be spilt or merged?
! 3. should ghosts be created?
! (2 and 3 are *separate* steps)

! set masses
massoftype(isplit) = massoftype(igas)/nchild
massoftype(ighost) = massoftype(igas)
massoftype(isplitghost) = massoftype(isplit)

! 1. delete existing ghosts
call delete_all_ghosts(xyzh,vxyzu,npartoftype,npart)

! 2. should it be split or merged?
merge_count = 0
add_npart = 0
do ii = 1,npart
  iphaseii = iphase(ii)
  call check_split_or_merge(ii,iphaseii,xyzh,vxyzu,npart,npartoftype,merge_count,add_npart)
enddo
npart = npart + add_npart
call shuffle_part(npart)

! 3. should a ghost be made?
merge_ghost_count = 0
add_npart = 0
do ii = 1,npart
  call check_ghost_or_splitghost(ii,iphase(ii),xyzh,vxyzu,npart,npartoftype,merge_ghost_count,add_npart)
enddo
npart = npart + add_npart

print*,'particle splitting is on, you now have',npart,'total particles'

end subroutine init_split

!----------------------------------------------------------------
!+
!  determines whether a particle should be split or merged, then
!  goes ahead and does the appropriate operation
!+
!----------------------------------------------------------------
subroutine check_split_or_merge(ii,iphaseii,xyzh,vxyzu,npart,npartoftype,merge_count,add_npart)
  use part, only:set_particle_type,iamtype,shuffle_part
  integer, intent(in)    :: ii
  integer(kind=1), intent(in) :: iphaseii
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer, intent(inout) :: npartoftype(:),npart
  integer, intent(inout) :: merge_count,add_npart
  logical :: already_split,split_it

  call inside_boundary(xyzh(1:3,ii),split_it)
  already_split = iamsplit(iphaseii)
  ! if it should be split but is not, split it
  if (split_it .and. .not.already_split) then

    call split_a_particle(nchild,ii,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)
    add_npart = add_npart + nchild - 1

  ! if it should not be split, but already is, merge it
  else if (.not.split_it .and. already_split) then
    merge_count = merge_count + 1
    children_list(merge_count) = ii
    if (merge_count == nchild) then
      call fast_merge_into_a_particle(nchild,children_list,npart, &
                                      xyzh,vxyzu,npartoftype,ii)
      merge_count = 0
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
subroutine check_ghost_or_splitghost(ii,iphaseii,xyzh,vxyzu,npart,npartoftype,merge_ghost_count,add_npart)
  use part, only:set_particle_type
  integer, intent(in)    :: ii
  integer(kind=1), intent(in) :: iphaseii
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer, intent(inout) :: npartoftype(:),npart
  integer, intent(inout) :: add_npart,merge_ghost_count
  logical :: already_split,ghost_it,already_ghost

  call inside_ghost_zone(xyzh(1:3,ii),ghost_it)
  already_ghost = iamghost(iphaseii)
  if (ghost_it .and. .not.already_ghost) then
    already_split = iamsplit(iphaseii)
    ! this is inside the boundary and should be merged (just use the original parent)
    if (already_split) then
      merge_ghost_count = merge_ghost_count + 1
      if (merge_ghost_count == nchild) then
        call make_a_ghost(npart+add_npart+1,ii,npartoftype,npart,nchild,xyzh)
        xyzh(4,npart+add_npart+1) = xyzh(4,npart+add_npart+1) * (nchild)**(1./3.)

        merge_ghost_count = 0
        add_npart = add_npart + 1
      endif
    ! this is outside the boundary and should be split
    else if (.not.already_split) then
      call make_split_ghost(npart+add_npart+1,ii,npartoftype,npart,nchild,xyzh,vxyzu)
      add_npart = add_npart + nchild
    endif
  endif
end subroutine check_ghost_or_splitghost
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
! create ghosts and splitghosts
!+
!-----------------------------------------------------------------------
subroutine make_all_ghosts(npart,npartoftype,nchild,xyzh,vxyzu)
  integer, intent(inout) :: npart,npartoftype(:)
  integer, intent(in)    :: nchild
  real, intent(inout)    :: xyzh(:,:),vxyzu(:,:)
  integer :: iighost,merge_count,ii
  integer(kind=1) :: iphaseii
  logical :: ghost_it,split_it,already_ghost,already_split

  iighost = npart
  merge_count = 0
  do ii = 1,npart
    call inside_ghost_zone(xyzh(1:3,ii),ghost_it)
    call inside_boundary(xyzh(1:3,ii),split_it)
    iphaseii = iphase(ii)
    already_ghost = iamghost(iphaseii)
    already_split = iamsplit(iphaseii)
    if (ghost_it .and. .not.already_ghost) then
      if (split_it) then ! this is inside the boundary and should be merged
        merge_count = merge_count + 1
        if (merge_count == nchild) then ! the lucky one that becomes a parent
          iighost = iighost + 1
          call make_a_ghost(iighost,ii,npartoftype,npart,nchild,xyzh)
          xyzh(4,iighost) = xyzh(4,iighost) * (nchild)**(1./3.)
          merge_count = 0
        endif
      else if (.not.already_split) then     ! this is outside the boundary and should be split
        call make_split_ghost(iighost+1,ii,npartoftype,npart,nchild,xyzh,vxyzu)
        iighost = iighost + nchild
      endif
    endif
  enddo

end subroutine

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
    if (gzw < tiny(gzw)) call fatal(label,'gzw must be > 0. in input options')
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
