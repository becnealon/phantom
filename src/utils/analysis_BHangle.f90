!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to crudely estimate the change in angular momentum of BH
! when using a fixed potential
!
! :References: Gerosa et al. 2020
!
! :Owner: Bec Nealon and Enrico Ragusa
!
! :Runtime parameters: None
!
! :Dependencies: part
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'BHangle'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only:massoftype,igas,xyzmh_ptmass,vxyz_ptmass
 use vectorutils,  only:cross_product3D,rotatevec
 use units,        only:utime,udist
 use physcon,      only:pi,years,pc
 use extern_lensethirring, only:sin_spinangle,cos_spinangle,blackhole_spin
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: particlemass,time
 integer, parameter :: iwrite = 24
 integer    :: ii
 logical    :: iexist
 real, save :: t_old = -1.,ctheta_pred
 real       :: sum_dJdt(3),pmass,fixed_BH(3),Li(3),r,dJdti(3)
 real       :: xi(3),vi(3),dJdt_mag,dJdt_unit(3),L_star(3),dcthetadt,dt,ctheta
 real       :: L_star_unit(3),R_star(3),R_star_mag,R_star3,temp(3)
 real       :: gamma,t_b,R_b,t_inspiral,f,L_tot(3),fixed_BH_unit(3),vec_mag
 real       :: J(3),rotate_about_z,rotate_about_y,L_unit(3),J_unit(3)
 character(len=40) :: filename

 ! initialise
 if (t_old < 0.0) then
   t_old = time
   ctheta_pred = cos_spinangle
 endif
 dt = time - t_old

 ! calculate dJ/dt, equation 37 of Gerosa et al. 2020
 ! NB: G=c=1 because iexternalforce=9,11
 fixed_BH = (/blackhole_spin*sin_spinangle,0.,blackhole_spin*cos_spinangle/)
 fixed_BH_unit = fixed_BH/sqrt(dot_product(fixed_BH,fixed_BH))
 sum_dJdt = 0.
 pmass = massoftype(igas)
 L_tot = 0.

 do ii=1,npart
   xi = xyzh(1:3,ii)
   vi = vxyzu(1:3,ii)
   r = sqrt(dot_product(xi,xi))
   call cross_product3D(xi,vi,Li)
   Li = pmass*Li
   call cross_product3D(fixed_BH,Li,dJdti)

   sum_dJdt = sum_dJdt + dJdti/(r**2)
   L_tot = L_tot + Li
 enddo
 sum_dJdt = -sum_dJdt*4.*pi

 ! calculate dctheta/dt, equation 38
 dJdt_mag = sqrt(dot_product(sum_dJdt,sum_dJdt))
 dJdt_unit = sum_dJdt/dJdt_mag
 call cross_product3D(xyzmh_ptmass(1:3,1),vxyz_ptmass(1:3,1),L_star)
 L_star = xyzmh_ptmass(4,1)*L_star
 vec_mag = sqrt(dot_product(L_star,L_star))
 L_star_unit = L_star/vec_mag
 dcthetadt = dot_product(dJdt_unit,L_star_unit)

 ! actual angle between disc and BH, cos(theta)
 ! (equivalent answer achieved with star rather than disc)
 vec_mag = sqrt(dot_product(L_tot,L_tot))
 L_unit = L_tot/vec_mag
 ctheta = dot_product(fixed_BH_unit,L_unit)
 ! predicted angle, from Eq 38
 ctheta_pred = ctheta_pred + dt*dcthetadt
 if (abs(ctheta_pred) > 1.0) ctheta_pred=1.0

 ! calculate R_star^-3
 ! this is a proxy for kappa, which is R_star^-3*(bunch of constants set by IC)
 R_star = xyzmh_ptmass(1:3,1)
 R_star_mag = sqrt(dot_product(R_star,R_star))
 R_star3 = R_star_mag**(-3)

 ! calculate new R_star (forcing inspiral), equations 47 and 48
 gamma = 1.5
 t_b = 1.0E6 * years/utime
 R_b = 0.05 * pc/udist
 f = 0.1
 t_inspiral = t_b/f * (R_star_mag/R_b)**gamma
 R_star = R_star - (dt*R_star/t_inspiral)
 xyzmh_ptmass(1:3,1) = R_star

 ! rotate the disc and star *as though* the BH has moved for next step
 ! rotation is done in 2 stages: we rotate around the z axis
 ! to be coincident with the x axis, then rotate about the y axis
 ! we look for the rotation required to move the new J
 ! back to it's original fixed position

J = fixed_BH + sum_dJdt*dt !where J "should" be
temp = (/J(1),J(2),0./)
vec_mag = sqrt(dot_product(temp,temp))
! rotation in the x-y plane, lines J up with x-axis
rotate_about_z = -acos(dot_product((/1.,0.,0./),temp/vec_mag))*sign(1.0,temp(2))
call rotatevec(J,(/0.,0.,1.0/),rotate_about_z)

vec_mag = sqrt(dot_product(J,J))
J_unit = J/vec_mag
! rotation in z-x plane, lines J up with fixed_BH position
rotate_about_y = -acos(dot_product((fixed_BH_unit),J_unit))
! this is required if the BH is retrograde
if (J_unit(1) > 0. .and. cos_spinangle < 0.) rotate_about_y = -rotate_about_y
call rotatevec(J,(/0.,1.,0./),rotate_about_y)

! now we have angles, rotate the star
call rotate_it(xyzmh_ptmass(1:3,1),vxyz_ptmass(1:3,1),rotate_about_z,rotate_about_y)

! rotate the gas particles
 do ii = 1,npart
   xi = xyzh(1:3,ii)
   vi = vxyzu(1:3,ii)

   call rotate_it(xi,vi,rotate_about_z,rotate_about_y)

   xyzh(1:3,ii) = xi
   vxyzu(1:3,ii) = vi
 enddo

 ! update values for next time-step
 t_old = time

 ! output to file
 filename='BH_angle.out'
 inquire(file=filename,exist=iexist)
 if (.not.iexist .or. time < tiny(time)) then
   open(unit=iwrite,file=filename,status="replace")
   write(iwrite,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
              1,'time', &
              2,'R_*^-^3', &
              3,'|dJ/dt|', &
              4,'theta pred', &
              5,'theta'
 else
   open(unit=iwrite,file=filename,status="old",position="append")
 endif
 write(iwrite,'(5(es18.10,1X))') time,R_star3,dJdt_mag,&
                acos(ctheta_pred)*180./pi,acos(ctheta)*180./pi
 close(unit=iwrite)

 write(*,*) "End of analysis BHangle"

 return

end subroutine do_analysis

  ! Rotate position and velocity vectors about z and y axis

subroutine rotate_it(pos,vel,z_angle,y_angle)
  use vectorutils,  only:rotatevec
    real, intent(inout) :: pos(3),vel(3)
    real, intent(in)    :: z_angle,y_angle

  call rotatevec(pos,(/0.,0.,1.0/),z_angle)
  call rotatevec(pos,(/0.,1.0,0./),y_angle)
  call rotatevec(vel,(/0.,0.,1.0/),z_angle)
  call rotatevec(vel,(/0.,1.0,0./),y_angle)

end subroutine rotate_it

end module analysis
