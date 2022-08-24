!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module relax_particles
  !
  ! Relaxing particles to a reference distribution (described by a separate set of particles)
  !
  ! :References: None
  !
  ! :Owner: Rebecca Nealon
  !
  ! :Runtime parameters: None
  !
  ! :Dependencies: None
  !
  implicit none

contains

  ! Subroutine to relax the new set of particles to a reference particle distribution
  ! Doesn't care if there's a mix of phases or not

  subroutine relaxparticles(npart,xyzh,npart_ref,xyzh_ref,force_ref,pmass_ref,n_toshuffle,to_shuffle)
    use deriv,       only:get_derivs_global
    use part,        only:massoftype,igas,fxyzu,isplit,iphase,iamghost
    integer,           intent(in)    :: npart,npart_ref,n_toshuffle
    real,              intent(in)    :: force_ref(3,npart_ref),xyzh_ref(5,npart_ref)
    real,              intent(in)    :: pmass_ref(npart_ref)
    integer,           intent(in)    :: to_shuffle(n_toshuffle)
    real,              intent(inout) :: xyzh(:,:)
    real,  allocatable :: a_ref(:,:)
    real :: ke,maxshift,pmassi
    logical :: converged,write_output,a_ghost
    integer :: ishift,ii,iref,igoal,nshifts

    if (n_toshuffle == 0) return ! it shouldn't have been called if this was the case

    print*,'Relaxing',n_toshuffle,' particles the heavenly way from',npart_ref,'references.'

    ! Initialise for the loop
    converged = .false.
    ishift = 0
    nshifts = 200

    ! If output required
    write_output = .true.
    if (write_output) then
      iref = 12
      open(unit=iref,file='heavenly_shuffling.dat',status='replace',position='append')
    endif

    ! a_ref stores the accelerations at the locations of the new particles as interpolated from the old ones
    allocate(a_ref(3,npart))

    ! This gets fxyz of the new particles at their new locations
    call get_derivs_global()

    do while (.not.converged)

      ! These are the accelerations at the locations of the new particles, interpolated from the parents
      call get_reference_accelerations(npart,xyzh,a_ref,npart_ref,xyzh_ref,force_ref,pmass_ref,n_toshuffle,to_shuffle)

      ! Shift the particles by minimising the difference between the acceleration at the new particles and
      ! the interpolated values (i.e. what they do have minus what they should have)
      call shift_particles(npart,xyzh,a_ref,n_toshuffle,to_shuffle,ke,maxshift)

      ! Todo: cut-off criteria
      if (ishift >= nshifts) converged = .true.
      ishift = ishift + 1

    enddo

    ! Tidy up
    deallocate(a_ref)
    if (write_output) close(iref)

  end subroutine relaxparticles

  !----------------------------------------------------------------
  !+
  ! Calculates the shift for each particle, using the reference
  ! and interpolated accelerations
  ! (based off the routine in relax_star)
  !+
  !----------------------------------------------------------------

  subroutine shift_particles(npart,xyzh,a_ref,n_toshuffle,to_shuffle,ke,maxshift)
    use dim,      only:periodic
    use part,     only: rhoh,vxyzu,fxyzu,massoftype,isplit,igas,iamtype,iphase
    use eos,      only: get_spsound
    use deriv,    only:get_derivs_global
    use options,  only:ieos
    use boundary, only:cross_boundary
    use domain,   only:isperiodic
    integer, intent(in)     :: npart,n_toshuffle
    real,    intent(inout)  :: xyzh(4,npart)
    real,    intent(in)     :: a_ref(3,npart)
    integer, intent(in)     :: to_shuffle(n_toshuffle)
    real,    intent(out)    :: ke,maxshift
    real :: hi,rhoi,cs,dti,dx(3),vi(3),err
    integer :: nlargeshift,i,iamtypei,ncross,j

    ke = 0.
    nlargeshift = 0
    ncross = 0
    maxshift = tiny(maxshift)
    !$omp parallel do schedule(guided) default(none) &
    !$omp shared(npart,xyzh,vxyzu,fxyzu,massoftype,ieos,iphase,a_ref,maxshift) &
    !$omp shared(isperiodic,ncross,to_shuffle,n_toshuffle) &
    !$omp private(i,dx,dti,cs,rhoi,hi,vi,iamtypei,err) &
    !$omp reduction(+:nlargeshift,ke)
    do j=1,n_toshuffle
      if (to_shuffle(j) == 0) cycle
      i = to_shuffle(j)
      iamtypei = iamtype(iphase(i))
      hi = xyzh(4,i)
      rhoi = rhoh(hi,massoftype(iamtypei))
      cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
      dti = 0.3*hi/cs   ! h/cs

      dx  = 0.5*dti**2*(fxyzu(1:3,i) - a_ref(1:3,i))
      if (sqrt(dot_product(dx,dx)) > maxshift) maxshift = sqrt(dot_product(dx,dx))
      if (dot_product(dx,dx) > hi**2) then

        dx = dx / sqrt(dot_product(dx,dx)) * hi  ! Avoid large shift in particle position !check with what James has done
        nlargeshift = nlargeshift + 1
      endif
      xyzh(1:3,i) = xyzh(1:3,i) + dx(:)
      if (periodic) call cross_boundary(isperiodic,xyzh(:,i),ncross)
      vi(1:3) = dx(:)/dti ! fake velocities, so we can estimate the magnitude of the shift
      ke = ke + dot_product(vi,vi)

      ! Output for testing purposes
      err = sqrt(dot_product(dx,dx))/hi
      if (err > maxshift) maxshift = err

    enddo
    !$omp end parallel do
    if (nlargeshift > 0) print*,'Warning: Restricted dx for ', nlargeshift, 'particles'

    !
    ! get forces on particles
    !
    call get_derivs_global()

  end subroutine shift_particles

  !----------------------------------------------------------------
  !+
  ! Interpolates the accelerations at the locations of the new particles
  ! from the old set of particles (the reference particles)
  !+
  !----------------------------------------------------------------

  subroutine get_reference_accelerations(npart,xyzh,a_ref,npart_ref,xyzh_ref,&
    force_ref,pmass_ref,n_toshuffle,to_shuffle)
    use dim,          only:periodic
    use kernel,       only:wkern,grkern,radkern2,cnormk
    use boundary,     only:dxbound,dybound,dzbound,xmax,xmin,ymax,ymin,zmax,zmin
    use part,         only:massoftype,rhoh,iamtype,iphase
    integer, intent(in) :: npart,npart_ref,n_toshuffle
    real,    intent(in) :: xyzh(4,npart),force_ref(3,npart_ref),xyzh_ref(5,npart_ref),pmass_ref(npart_ref)
    integer, intent(in) :: to_shuffle(n_toshuffle)
    real,    intent(out) :: a_ref(3,npart)
    real :: xi,yi,zi,hi,rij(3),h21,qj2,rij2,rhoj,h31,mass_ref,pmassi
    integer :: i,j,iamtypej,k

    a_ref(:,:) = 0.

    !$omp parallel do schedule(guided) default (none) &
    !$omp shared(xyzh,xyzh_ref,npart,npart_ref,force_ref,a_ref,pmass_ref,to_shuffle,n_toshuffle,massoftype,iphase) &
    !$omp shared(mass_ref,dxbound,dybound,dzbound) &
    !$omp private(i,j,xi,yi,zi,rij,h21,h31,rhoj,rij2,qj2,pmassi)

    ! Over the new set of particles that are to be shuffled
    over_new: do k = 1,n_toshuffle
      if (to_shuffle(k) == 0) cycle over_new
      i = to_shuffle(k)
      xi = xyzh(1,i)
      yi = xyzh(2,i)
      zi = xyzh(3,i)
      pmassi = massoftype(iphase(i))

      ! Over the reference set of particles to which we are matching the accelerations
      over_reference: do j = 1,npart_ref  ! later this should only be over active particles
        rij(1) = xyzh_ref(1,j) - xi
        rij(2) = xyzh_ref(2,j) - yi
        rij(3) = xyzh_ref(3,j) - zi
        mass_ref = pmass_ref(j)

        if (periodic) then
          if (abs(rij(1)) > 0.5*dxbound) rij(1) = rij(1) - dxbound*SIGN(1.0,rij(1))
          if (abs(rij(2)) > 0.5*dybound) rij(2) = rij(2) - dybound*SIGN(1.0,rij(2))
          if (abs(rij(3)) > 0.5*dzbound) rij(3) = rij(3) - dzbound*SIGN(1.0,rij(3))
        endif

        h21 = 1./(xyzh_ref(4,j))**2
        h31 = 1./(xyzh_ref(4,j))**3
        rhoj = rhoh(xyzh_ref(4,j),mass_ref)

        rij2 = dot_product(rij,rij)
        qj2  = rij2*h21
        if (qj2 < radkern2) then
          ! Interpolate acceleration at the location of the new particle
          a_ref(:,i) = a_ref(:,i) + force_ref(:,j)*wkern(qj2,sqrt(qj2))*cnormk*h31/rhoj
        endif

      enddo over_reference
    enddo over_new

    !$omp end parallel do

  end subroutine get_reference_accelerations

end module relax_particles
