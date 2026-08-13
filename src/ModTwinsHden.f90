!*******************************************************************************
!                         
!                           ModTwinsHden.f90
!
! This module contains subrutines that calculates exospheric hydrogen density
!  from TWINS data
!
! Originally created by Gonzalo Cucho-Padin, Code 674.
! Modified to a separate module by S.-B. Kang, Code 673.
!
!*******************************************************************************

  !module twins_mod
  module ModTwinsHden

  public :: & 
    NumHour,&
    twins_init,&
    twins_geocorona

  private ! except

  integer,parameter :: nparam = 33 ! number of parameters per hour
  integer,parameter :: ntmatrix = 9 ! number of elements in the transformation matix
  integer,parameter :: nhours = 72 ! number of hours in the file
  integer :: NumHour ! current hour
  real,dimension(:,:) :: ReadParam(0:71,1:nparam)
  real,dimension(:,:) :: ReadTMatrix(0:71,1:ntmatrix)
  real,dimension(:) :: param2beUsed(1:nparam)
  real,dimension(:) :: Tmatrix2beUsed(1:ntmatrix)

 contains
!-----------------------------------------------------------------------------
  subroutine twins_init(storm)
!-----------------------------------------------------------------------------
  implicit none
  integer :: i
  character*13   :: storm
  character*1000 :: StringHeader
  logical :: FileExist

  !Get twins parameter from a file
  inquire(file=trim(storm)//'_HDENS.dat',exist=FileExist)
  if (.not.FileExist) then
     write(*,*) '*** ERROR **************************'
     write(*,*) 'File, '//trim(storm)//'_HDENS.dat does NOT exist'
     write(*,*) '************************************'
  endif

  open(300, FILE=trim(storm)//'_HDENS.dat',status='old')
  read(300,*) StringHeader ! Reading the header (dummy read)

  !Reading parameters
  read(300,*) ReadParam(0,1:nparam) !only 1 set of parameters (static case)
  !write(*,*) 'let see params',ReadParam(0,1:nparam)

  ! Get values for transformation matrix
  read(300,*) StringHeader ! Readung the header (dummy read)

  ! Reading transformation matrix values
  do i=0,nhours-1
        read(300,*) ReadTMatrix(i,1:ntmatrix)
  end do
  !write(*,*) 'let see matrixvals',ReadTMatrix(0,1:ntmatrix)
  !Initial value for NumHour
  NumHour=0

  !Initial values for param2beUsed
  param2beUsed(1:nparam) = ReadParam(NumHour,1:nparam)

  !Initial values for Tmatrxi2beUsed
  Tmatrix2beUsed(1:ntmatrix)=ReadTMatrix(NumHour,1:ntmatrix)

  !Closing file
  close(300)

  end subroutine twins_init

!-----------------------------------------------------------------------------
  subroutine twins_geocorona(x,y,z,hden)
!-----------------------------------------------------------------------------
! Geocorona model based on  twins data for 2015
  real :: x,y,z,hden,rad,radxy,cos_phi,cos_theta,phi_h,theta,part1,part2
  real :: cos2theta,sin2theta,cos3theta,sin3theta
  real :: Pi=3.14159265358979  ! constants

! Parameters for static quiet time H density using TWINS data fitted to SPH 3rd order
  param2beUsed(1:nparam) = ReadParam(0,1:nparam)  !params won't change in this version
  Tmatrix2beUsed(1:ntmatrix) = ReadTMatrix(NumHour,1:ntmatrix)
 
! Converting SM values to GSE
  ! fixed by SBK array index starting from 1
  x_gse = x*Tmatrix2beUsed(1) + y*Tmatrix2beUsed(4) + z*Tmatrix2beUsed(7)
  y_gse = x*Tmatrix2beUsed(2) + y*Tmatrix2beUsed(5) + z*Tmatrix2beUsed(8)
  z_gse = x*Tmatrix2beUsed(3) + y*Tmatrix2beUsed(6) + z*Tmatrix2beUsed(9)
  ! previously written by GCP
  !x_gse = x*Tmatrix2beUsed(0) + y*Tmatrix2beUsed(3) + z*Tmatrix2beUsed(6)
  !y_gse = x*Tmatrix2beUsed(1) + y*Tmatrix2beUsed(4) + z*Tmatrix2beUsed(7)
  !z_gse = x*Tmatrix2beUsed(2) + y*Tmatrix2beUsed(5) + z*Tmatrix2beUsed(8)

! Geometric calculations
  rad = sqrt(x_gse*x_gse+y_gse*y_gse+z_gse*z_gse)
  radxy = sqrt(x_gse*x_gse+y_gse*y_gse)
  cos_phi = x_gse / radxy
  cos_theta = z_gse / rad
  phi_h = acos(cos_phi)
  theta = acos(cos_theta)
  cos2theta = (cos(theta))**2
  sin2theta = (sin(theta))**2
  cos3theta = (cos(theta))**3
  sin3theta = (sin(theta))**3

! Formula for SPH 3rd Order
  part1 = param2beUsed(1) * (rad**(-param2beUsed(2)))*(param2beUsed(3)**(1/rad)) * sqrt(4.*Pi)
  part2 = sqrt(1./(4. * Pi)) +  &
 (param2beUsed(4) + param2beUsed(5)*log(rad))*sqrt(3./(4. * Pi)) * cos(theta) + &
 (param2beUsed(6) + param2beUsed(7)*log(rad))*cos(phi_h)*(-sqrt(3./(8. * Pi)))*sin(theta) + &
 (param2beUsed(8) + param2beUsed(9)*log(rad))*sin(phi_h)*(-sqrt(3./(8. * Pi)))*sin(theta) + &
 (param2beUsed(10) + param2beUsed(11)*log(rad))*sqrt(5./(4.*Pi)) * (1.5 * cos2theta - 0.5) + &
 (param2beUsed(12) + param2beUsed(13)*log(rad))*cos(phi_h)*(-sqrt(15./(8.*Pi))*sin(theta)*cos(theta)) + &
 (param2beUsed(14) + param2beUsed(15)*log(rad))*sin(phi_h)*(-sqrt(15./(8.*Pi))*sin(theta)*cos(theta)) + &
 (param2beUsed(16) + param2beUsed(17)*log(rad))*cos(2.*phi_h)*(0.25 *sqrt(15./(2.*Pi))*sin2theta) + &
 (param2beUsed(18) + param2beUsed(19)*log(rad))*sin(2.*phi_h)*(0.25 *sqrt(15./(2.*Pi))*sin2theta) + &
 (param2beUsed(20) + param2beUsed(21)*log(rad))*sqrt(7./(4. * Pi))*((2.5*cos3theta) - 1.5*cos(theta)) + &
 (param2beUsed(22) + param2beUsed(23)*log(rad))*cos(phi_h)*(-0.25)*sqrt(21./(4.*Pi))*sin(theta)*(5.*cos2theta-1) + &
 (param2beUsed(24) + param2beUsed(25)*log(rad))*sin(phi_h)*(-0.25)*sqrt(21./(4.*Pi))*sin(theta)*(5.*cos2theta-1) + &
 (param2beUsed(26) + param2beUsed(27)*log(rad))*cos(2.*phi_h)*(0.25)*sqrt(105./(2.*Pi))*sin2theta*cos(theta) + &
 (param2beUsed(28) + param2beUsed(29)*log(rad))*sin(2.*phi_h)*(0.25)*sqrt(105./(2.*Pi))*sin2theta*cos(theta) + &
 (param2beUsed(30) + param2beUsed(31)*log(rad))*cos(3.*phi_h)*(-0.25)*sqrt(35./(4.*Pi))*sin3theta + &
 (param2beUsed(32) + param2beUsed(33)*log(rad))*sin(3.*phi_h)*(-0.25)*sqrt(35./(4.*Pi))*sin3theta
 
 hden = part1*part2

 if (hden.lt.0) hden=0.

  end subroutine twins_geocorona
!-----------------------------------------------------------------------------

  end module ModTwinsHden
