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
 public :: check_split_or_merge,check_ghost_or_splitghost,compare_reals_and_ghosts
 public :: update_splitting,merge_particles_wrapper

 private

 integer, parameter, public  :: nchild = 12
 integer, save    :: children_list(nchild+1) = 0
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
 use splitpart, only: merge_all_particles
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
if (centre_particle) then
   massoftype(isplit) = massoftype(igas)/(nchild+1)
else
   massoftype(isplit) = massoftype(igas)/(nchild)
endif
massoftype(ighost) = massoftype(igas)
massoftype(isplitghost) = massoftype(isplit)

! 1. delete existing ghosts
call delete_all_ghosts(xyzh,vxyzu,npartoftype,npart)

! 2. should it be split or merged?

 if (rand_type==5) then
    call update_splitting(npart,xyzh,vxyzu,npartoftype)
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
   call merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu,13,npart,.true.)
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
subroutine update_splitting(npart,xyzh,vxyzu,npartoftype)
 use io,   only: fatal
 use part, only:kill_particle,isdead_or_accreted,shuffle_part,iactive
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,k,add_npart
 logical                :: split_it,ghost_it,already_split,already_ghost,make_ghost
 integer, allocatable, dimension(:)   :: iorig(:)
 real,    allocatable, dimension(:,:) :: xyzh_split(:,:)

 if (rand_type /= 5) then
    call fatal('update_splitting','this routine is set up only for rand_type==5')
    stop
 endif
 !
 !--Kill all active ghost particle  &
 !  Convert gas <-> splits that have crossed the boundary
 !
 allocate(xyzh_split(4,npartoftype(isplit)))
 allocate(iorig(npartoftype(isplit)))
 k = 0
 add_npart = 0
 make_ghost = .false.
 split_loop: do i = 1,npart
    if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       if (iamghost(iphase(i))) then
          call kill_particle(i,npartoftype)
          cycle split_loop
       endif
       call inside_boundary(xyzh(1:3,i),split_it)
       already_split = iamsplit(iphase(i))
       if (split_it .and. .not.already_split) then
          call split_a_particle(nchild,i,xyzh,vxyzu,npartoftype,0,1,npart+add_npart)  ! nchild = 12
          add_npart = add_npart + nchild
       elseif (.not.split_it .and. already_split) then ! this should be only the split particles that need to be converted into gas particles
          k = k + 1
          xyzh_split(:,k) = xyzh(:,i)
          iorig(k) = i
       endif
    endif
 enddo split_loop
 npart = npart + add_npart
 print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'
 if (k > 0) call merge_particles(npart,k,xyzh,xyzh_split,iorig,npartoftype,make_ghost)
 call shuffle_part(npart)
 deallocate(xyzh_split)
 deallocate(iorig)
 !
 !--Go through list and create ghost particles as required
 !
 allocate(xyzh_split(4,npartoftype(isplit)))
 allocate(iorig(npartoftype(isplit)))
 k = 0
 add_npart = 0
 make_ghost = .true.
 ghost_loop: do i = 1,npart
    if (iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       call inside_ghost_zone(xyzh(1:3,i),ghost_it)
       already_ghost = iamghost(iphase(i))
       if (ghost_it .and. .not.already_ghost) then
          already_split = iamsplit(iphase(i))
          ! this is inside the boundary and should be merged
          if (already_split) then
             ! make list of particles to make ghost particles
             k = k + 1
             xyzh_split(:,k) = xyzh(:,i)
             iorig(k) = i
          else
             ! make splitghosts
             call make_split_ghost(npart+add_npart+1,i,npartoftype,npart,nchild,xyzh,vxyzu)
             add_npart = add_npart + nchild + 1
          endif
       endif
    endif
 enddo ghost_loop
 npart = npart + add_npart
 print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'
 if (k > 0) call merge_particles(npart,k,xyzh,xyzh_split,iorig,npartoftype,make_ghost)
 deallocate(xyzh_split)
 deallocate(iorig)

end subroutine update_splitting

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
        ilucky = int(ran2(seed)*(nchild+1))+1
        if (ilucky==1) then
           ! particle becomes merged
           xyzh(4,ii) = xyzh(4,ii) * (nchild+1)**(1./3.)
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
       if (merge_count == nchild+1) then
         call fast_merge_into_a_particle(nchild+1,children_list,npart,xyzh,vxyzu,npartoftype,ii)
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
        ilucky = int(ran2(seed)*(nchild+1))+1
        if (ilucky==1) then
           call make_a_ghost(npart+add_npart+1,ii,npartoftype,npart,nchild+1,xyzh)
           add_npart = add_npart + 1
        endif
      else
         merge_ghost_count = merge_ghost_count + 1
         children_list(merge_ghost_count) = ii
         if (merge_ghost_count == nchild+1) then
            if (rand_type==1) then
               ilucky = int(ran2(seed)*(nchild+1))+1
               ilucky = children_list(ilucky)
               !print*, ilucky,children_list(ilucky)
            elseif (rand_type==0) then
               ilucky = ii
            else
               print*, 'This should not be possible'
               stop
            endif
            call make_a_ghost(npart+add_npart+1,ilucky,npartoftype,npart,nchild+1,xyzh)
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
 integer, allocatable, dimension(:)   :: iorig(:)
 real,    allocatable, dimension(:,:) :: xyzh_split(:,:)
 integer :: i,k
 logical :: ghost_it,split_it

 allocate(xyzh_split(4,npartoftype(isplit))) ! we can be smarter and choose smaller arrays later
 allocate(iorig(npartoftype(isplit)))        ! we can be smarter and choose smaller arrays later

 ! Determine the number of particles
 k = 0
 do i = 1,npart
    if (iphase(i)==isplit .and. iactive(iphase(i)) .and. .not.isdead_or_accreted(xyzh(4,i))) then
       if (make_ghost) then
          call inside_ghost_zone(xyzh(1:3,i),ghost_it)
          if (ghost_it) then ! this should be only the split particles that need to have ghosts made of them
             k = k + 1
             xyzh_split(:,k) = xyzh(:,i)
             iorig(k) = i
          endif
       else
          call inside_boundary(xyzh(1:3,i),split_it)
          if (.not. split_it) then ! this should be only the split particles that need to be converted into gas particles
             k = k + 1
             xyzh_split(:,k) = xyzh(:,i)
             iorig(k) = i
          endif
       endif
    endif
 enddo
 print*, 'make_ghost = ',make_ghost,' using ',k,' particles.'

 if (k > 0) call merge_particles(npart,k,xyzh,xyzh_split,iorig,npartoftype,make_ghost)

 deallocate(xyzh_split)
 deallocate(iorig)

end subroutine merge_particles_wrapper
!----------------------------------------------------------------
!+
!  This is a subroutine to merge particles, either to create a gas
!  particle or a gas ghost particle
!  if make_ghosts = true, then make ghosts from the split particles
!  if make_ghosts = false, then merge split particles into gas particles
!+
!----------------------------------------------------------------
subroutine merge_particles(npart,ncandiate,xyzh,xyzh_split,iorig,npartoftype,make_ghost)
 use linklist,       only:ncells,get_neighbour_list,get_hmaxcell,get_cell_location,set_linklist,ifirstincell
 use kdtree,         only:inodeparts,inoderange
 use part,           only:iactive,isdead_or_accreted,set_particle_type,kill_particle
 use mpiforce,       only:cellforce
 integer, intent(inout) :: npartoftype(:),npart
 real,    intent(inout) :: xyzh(:,:),xyzh_split(:,:)
 logical, intent(in)    :: make_ghost
 integer, intent(in)    :: iorig(:)
 integer, intent(inout) :: ncandiate  ! needs to be inout to satisfy if MPI, but we can't have pre-processor statements in this .f90 file
 integer                :: i,j,k,icell,jmin,n_cell,iave(4)
 real                   :: r2,r2min,hmax,hcell
 type(cellforce)        :: cell

 ! Use the tree to make  the link list
 call set_linklist(ncandiate,ncandiate,xyzh_split(:,1:ncandiate),xyzh(:,1:ncandiate)) ! the final term is never used in linklist

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
    ! make a gas ghost from the splits OR merge splits into gas
    if (make_ghost) then
       k = iorig(inodeparts(jmin))
       npart = npart + 1
       call make_a_ghost(npart,k,npartoftype,npart,13,xyzh)
    else
       do j = inoderange(1,icell),inoderange(2,icell)
          k = iorig(inodeparts(j))
          if (j==jmin) then
             xyzh(4,k) = xyzh(4,k) * 13.**(1./3.)
             call set_particle_type(k,igas)
             npartoftype(igas)   = npartoftype(igas)   + 1
             npartoftype(isplit) = npartoftype(isplit) - 1
          else
             call kill_particle(k,npartoftype(:))
          endif
       enddo
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

!-----------------------------------------------------------------------
!+
!  Testing routine to calculate the mean density for the particles in the
!  ghost region, prints out the percentage difference between ghosts and
!  reals in comparable region
!+
!-----------------------------------------------------------------------

subroutine compare_reals_and_ghosts(xyzh,npart)
  use part, only:massoftype,rhoh,iamtype,igas,isplit,ighost,isplitghost
  real, intent(in) :: xyzh(:,:)
  integer, intent(in) :: npart
  real :: m_gas,m_splits,m_ghosts,m_splitghosts
  integer :: ii,itype,gas, splits,ghosts,splitghosts
  logical :: ghost_it

  m_gas = 0.
  gas = 0
  m_splits = 0.
  splits = 0
  m_ghosts = 0.
  ghosts = 0
  m_splitghosts = 0.
  splitghosts = 0

  do ii = 1,npart
    call inside_ghost_zone(xyzh(1:3,ii),ghost_it)
    if (ghost_it) then
      itype = iamtype(iphase(ii))
      if (itype==igas) then
          m_gas = m_gas + rhoh(xyzh(4,ii),massoftype(igas))
          gas = gas + 1
      elseif (itype==isplit) then
          m_splits = m_splits + rhoh(xyzh(4,ii),massoftype(isplit))
          splits = splits + 1
      elseif (itype==ighost) then
          m_ghosts = m_ghosts + rhoh(xyzh(4,ii),massoftype(ighost))
          ghosts = ghosts + 1
      elseif (itype==isplitghost) then
          m_splitghosts = m_splitghosts + rhoh(xyzh(4,ii),massoftype(isplitghost))
          splitghosts = splitghosts + 1
      endif
    endif
  enddo

  m_gas = m_gas/gas
  m_splits = m_splits/splits
  m_ghosts = m_ghosts/ghosts
  m_splitghosts = m_splitghosts/splitghosts

  print*,'-----------------------------------------------'
  print*,'Outside boundary',(m_gas/m_splitghosts)*100.
  print*,'Inside boundary',(m_splits/m_ghosts)*100.
  print*,'-----------------------------------------------'
end subroutine compare_reals_and_ghosts

end module split
