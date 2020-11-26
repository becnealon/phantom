!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module splitpart
!
! None
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: getneighbours, io, part, splitmergeutils, timestep_ind
!

 use splitmergeutils, only:split_a_particle,fancy_merge_into_a_particle,make_a_ghost

 implicit none

contains

subroutine split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                 nchild,lattice_type,ires)
 use io,    only:fatal,error
 use part,  only:igas,copy_particle,shuffle_part,set_particle_type,isplit
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in)    :: lattice_type,ires,nchild
 integer :: ierr,ichild,iparent,j

 ierr = 0

 !--check there is enough memory
 if (size(xyzh(1,:)) < npart*nchild) then
    call error('split_all_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 !--find positions of the new particles
 ichild = npart !to keep track of the kids
 do iparent=1,npart
    ! send in the parent, children return
    ! (the parent acts as the first child, this routine generates nchild-1 new particles
    ! and adjusts the smoothing length on the parent)
    call split_a_particle(nchild,iparent,xyzh,vxyzu,npartoftype,lattice_type,ires,ichild)

    ! undo type, as we don't want to continue simulation with split particles
    do j=1,nchild
      call set_particle_type(ichild+j,igas)
    enddo

    ! for next children
    ichild = ichild + nchild - 1
 enddo

!-- new npartoftype
npartoftype(igas) = npartoftype(isplit)
npartoftype(isplit) = 0

 !-- new npart
 npart = npart * nchild

 call shuffle_part(npart)

 !--new masses
 massoftype(:) = massoftype(:)/nchild

 if (ierr /= 0) call fatal('splitpart','could not split particles')

end subroutine split_all_particles

subroutine merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                                nchild,nactive_here,fancy_merging)
 use part,          only:igas,isplit,iphase,kill_particle,delete_dead_or_accreted_particles
 use part,          only:isdead_or_accreted,copy_particle
 use timestep_ind,  only:nactive
 use io,            only:fatal,error
 use getneighbours, only:generate_neighbour_lists,neighb,neighcount,neighmax
 integer, intent(inout) :: npart,npartoftype(:)
 integer, intent(in)    :: nchild
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, optional, intent(in) :: fancy_merging
 integer, optional, intent(in) :: nactive_here
 real,    allocatable, dimension(:,:) :: xyzh_split(:,:)
 integer :: ierr,nparent,remainder, i,k,npart0
 integer :: on_list(2*npart), children_list(nchild)
 integer :: neighbours(neighmax),neigh_count
 integer :: ichild,iparent,child_found,m
 real    :: rik(3),rik2,rikmax,r2_neighbours(neighmax)
 logical :: merge_stochastically,found
 logical :: ghostify = .true.

 ierr = 0

 !-- how should the particles be merged? stochastic is fastest but not *quite*
 !   as accurate as averaging properties over the children
 merge_stochastically = .true.
 if (present(fancy_merging)) then
    merge_stochastically = .false.
    print*,'Doing fancy merging, this takes a bit longer'
 endif

 !-- how many active particles? If called from moddump, it must be provided
 if (present(nactive_here)) nactive = nactive_here

 if (nactive < npart) then
    call delete_dead_or_accreted_particles(npart,npartoftype)
    print*,'Discarding inactive particles'
 endif

 if (ghostify) then
    allocate(xyzh_split(4,npart))
    xyzh_split(:,1:npart) = xyzh(:,1:npart)
    nactive = 0
    do i = 1,npart
       if (iphase(i)/=isplit) then
          xyzh_split(4,i) = 0. ! so we only find neighbours of split particles
       else
          nactive = nactive + 1
       endif
    enddo
 endif
 print*, 'nactive,nchild=',nactive,nchild
 !-- check number of parent particles
 nparent = floor(real(nactive)/real(nchild))
 remainder = mod(nactive,nchild)
 if (remainder/nactive > 0.01) then
    ! discarding a couple of particles is ok, just don't want it to be too many
    call error('merge_particles','need to merge evenly, make sure npart(active)/nchild ~ integer')
    ierr = 1
    return
 elseif (remainder > 0) then
    ! we just forget about these later on
    print*,' ignoring',remainder,'particles to get an even split'
 endif

 !--check there is enough memory
 if (size(xyzh(1,:)) < nparent + npart) then
    call error('merge_particles','not enough memory, increase MAXP and recompile')
    ierr = 1
    return
 endif

 iparent = 0
 ichild  = 0
 on_list = -1
 i       = 0
 npart0  = npart
 children_list = 0
 neighbours    = 0
 neigh_count   = 0

 if (merge_stochastically) then
    !-- quick stochastic merging
    do i = 1,npart,nchild
       iparent = iparent + 1
       call copy_particle(i,npart+iparent)
       xyzh(4,npart+iparent) = xyzh(4,i) * (nchild)**(1./3.)
    enddo
 else
    !-- slower merging that averages properties of children to find new parent
    !-- find all the neighbours first
    if (ghostify) then
       call generate_neighbour_lists(xyzh_split,vxyzu,npart,'dummy',write_neighbour_list=.false.)
    else
       call generate_neighbour_lists(xyzh,vxyzu,npart,'dummy',write_neighbour_list=.false.)
    endif
    print*, 'Done generating neighbours'
    print*, 'nparent = ',nparent
    !-- now find the children
!!!    over_parent: do while (iparent < nparent) ! until all the children are found  !if (.not.ghostify)
    i = 0
    over_parent: do while (i < npart0) ! if (ghostify)
       i = i + 1
       ! already on the list
       if (xyzh_split(4,i) < epsilon(xyzh_split(4,i))) cycle over_parent! if (ghostify)
       if (on_list(i) > 0) cycle over_parent


       ! first child
       iparent = iparent + 1
       ichild = 1
       on_list(i) = i
       children_list(ichild) = i
       neigh_count = neighcount(i)
       neighbours(:) = neighb(i,:)

       ! calculate distance to neighbours
       do k = 1, neigh_count
          m = neighbours(k)
          rik = xyzh(1:3,m) - xyzh(1:3,i)
          r2_neighbours(k) = dot_product(rik,rik)
       enddo

       finding_children: do while (ichild < nchild)
          found = .false.
          child_found = -1
          rikmax = huge(rikmax)

          !first check over the neighbours
          if (.not.found) then
             over_neighbours: do k=1, neigh_count
                m = neighbours(k)
                if (on_list(m) > 0) cycle over_neighbours
                if (r2_neighbours(k) < rikmax) then
                   rikmax = r2_neighbours(k)
                   child_found = m
                endif
             enddo over_neighbours
             if (child_found > 0) found = .true.
          endif

          ! otherwise, check everyone else
          if (.not.found) then
             over_all: do k = 1,npart0
                if (on_list(k) > 0) cycle over_all
                rik = xyzh(1:3,k) - xyzh(1:3,i)
                rik2 = dot_product(rik,rik)
                if (rik2 < rikmax) then
                   rikmax = rik2
                   child_found = k
                endif
             enddo over_all
             if (child_found > 0) found = .true.
          endif
          ! error if no children found for the parent particle
          if (.not.found) then
             call error('mergepart','no child found for parent particle')
             ierr = 1
          endif

          ! save the child to the list
          ichild = ichild + 1
          children_list(ichild) = child_found
          on_list(child_found)  = child_found
       enddo finding_children


       if (ghostify) then
          ! copy central particle into a ghost
          npart = npart + 1
          call make_a_ghost(npart,i,npartoftype,npart,13,xyzh)
       else
          ! send in children, parent returns
          ! parents temporarily stored after all the children
          call fancy_merge_into_a_particle(nchild,children_list, massoftype(igas), &
                                  npart,xyzh,vxyzu,npartoftype,npart+iparent)
       endif
    enddo over_parent
 endif

 if (ghostify) then
    deallocate(xyzh_split)
 else
    !-- move the new parents
    do i = 1,nparent
       call copy_particle(npart+i,i)
    enddo

    !-- kill all the useless children
    do i=nparent+1,npart
       call kill_particle(i,npartoftype)
    enddo

    !--tidy up
    call delete_dead_or_accreted_particles(npart,npartoftype)

    !--update npartoftype
    npartoftype(igas) = nparent
    npart = nparent
    massoftype(:) = massoftype(:) * nchild
 endif
 print*, 'DONE merge_all_particles'
 if (ierr /= 0) call fatal('moddump','could not merge particles')

! I think this method is perfect by fluke.  To split, we convert a gas into a split, and add 12 splits at the end of the list.
! In this case, since we go through i = 1,npart, we naturally run through all the former parents first, and remove their children.
! So due to the placement of the children on the list, the above routine perfectly undoes what the merging did.  This will not happen
! so smoothly in a live calculation.

end subroutine merge_all_particles

end module splitpart
