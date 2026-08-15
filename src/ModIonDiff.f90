!*******************************************************************************
!                           ModIonDiff.f90
!
! Module contains Ion diffusion in momentum space, which includes diffuse due to
!  EMIC waves.
!
!*******************************************************************************
  module ModIonDiff
    implicit none

    public :: &
      ! subroutines
      init_ion_diff,&
      diffuse_UV_ions,& 
      unit_test_diffuse_ions,& 
      ! variables
      iEmicDiff,&
      EmicHPowerInPs,&  ! H band
      EmicHePowerInPs,& ! He band
      EmicOPowerInPs,&  ! O band
      EmicHPowerOutPs,&  ! H band
      EmicHePowerOutPs,& ! He band
      EmicOPowerOutPs    ! O band

    private   ! except

    !\
    ! variables from cimi
    integer :: ir,ip,iw,ik,ns,ijs
    integer,allocatable :: js(:)
    real,allocatable :: &
               xk(:),lnp(:,:),dlnp(:),gridp(:,:),ekev(:,:)
    real :: dt,dlnK,xme

    integer nKp
    real,allocatable :: tKp(:),zkpa(:)

    ! constants
    real :: pi=3.14159265358979
    real    :: re_m=6.3712e6	            ! earth's average radius (m)

    ! EMIC wave difusion
    integer :: iEmicDiff=1

    integer :: ipEmic,&   ! size of fpe/fce array for EMIC wave diffusion   
               iePA       ! size of energy arrary for "
    real,allocatable,dimension(:) :: &
               EmicOmpe,& ! fpe/fce array
               BEmic,&    ! EMIC wave band array (H+,He+,O+)
               ePA        ! pitch-angle array
    real,allocatable,dimension(:,:,:,:) :: &
         EmicHDaa,EmicHDpp,EmicHDap,&    ! H-band EMIC waves' Daa, Dap, Dpp
         EmicHeDaa,EmicHeDpp,EmicHeDap,& ! He-band
         EmicODaa,EmicODpp,EmicODap      ! O-band
                   !dimension: ns, fpe/fce, iw, ik  
    ! EMIC waves amplitude (nT) 
    real EmicPower0
    real,allocatable,dimension(:,:) :: &
    !  inside plasmasphere
         EmicHPowerInPs,&  ! H band
         EmicHePowerInPs,& ! He band
         EmicOPowerInPs ,& ! O band
    !  outside plasmasphere
         EmicHPowerOutPs ,& ! H band
         EmicHePowerOutPs,& ! He band
         EmicOPowerOutPs    ! O band

    ! CRRES EMIC wave data
  integer,parameter :: nCRRES=14
  real,dimension(nCRRES,3) :: EMIC_CRRES_H,EMIC_CRRES_He
  real,dimension(nCRRES) :: L_CRRES
  
  data EMIC_CRRES_H(:,1) /0.,0.,0.,0.    ,0.    ,0.     ,5.e-5  ,6.25e-5,1.25e-5,0.     ,0.     ,0.,0.,0./
  data EMIC_CRRES_H(:,2) /0.,0.,0.,2.5e-4,4.5e-4,4.75e-4,7.75e-4,1.e-3  ,1.3e-3 ,1.18e-3,4.25e-4,0.,0.,0./
  data EMIC_CRRES_H(:,3) /0.,0.,0.,3.e-3 ,8.5e-3,7.e-3  ,7.e-3  ,2.05e-2,2.5e-2 ,2.e-2  ,1.25e-2,0.,0.,0./
  data EMIC_CRRES_He(:,1) /0.,0.,0.,2.e-4 ,3.e-5 ,2.5e-4 ,2.5e-4,5.e-4,0.001 ,5.e-4 ,5.e-4,0.,0.,0./ 
  data EMIC_CRRES_He(:,2) /0.,0.,0.,3.5e-3,5.e-3 ,0.01   ,0.01  ,0.02 ,2.5e-2,2.5e-2,0.005,0.,0.,0./ 
  data EMIC_CRRES_He(:,3) /0.,0.,0.,5.e-4 ,7.e-2 ,0.1    ,0.1   ,0.1  ,8.e-2 ,3.e-2 ,0.001,0.,0.,0./
  data L_CRRES /1.,2.,3., 3.5 ,4. ,4.5 ,5. ,5.5 ,6. ,6.5 ,7. ,8.,9.,10./

  contains
!***********************************************************************
!                        init_ion_diff     
!
!  Routine receives parameters required to calculate pitch-angle
!   diffusion coefficient (Daa) due to field line curvature 
!   (curvature scattering) from cimi.
!***********************************************************************
  subroutine init_ion_diff(ns_,ir_,ip_,iw_,ik_,ijs_,js_,&
                           xk_,dlnK_,lnp_,dlnp_,gridp_,ekev_,&
                           nKp_,tKp_,zkpa_,xme_,dt_)
  integer,intent(in) :: ns_,ir_,ip_,iw_,ik_,ijs_,nKp_
  integer,intent(in) :: js_(ns_)
  real,intent(in) :: xk_(0:ik_+1),dlnK_,lnp_(ns_,0:iw_+1),dlnp_(ns_),gridp_(ns_,0:iw_+1),&
                     ekev_(ns_,0:iw_+1),tKp_(nKp_),zkpa_(nKp_),xme_,dt_

  ns=ns_
  ir=ir_
  ip=ip_
  iw=iw_
  ik=ik_
  ijs=ijs_
  dt=dt_

  dlnK=dlnK_

  nKp=nKp_

  xme=xme_    ! dipole moment, (T m^3)

  if (.not.allocated(js)) allocate(js(ns))
  js(1:ns)=js_(1:ns)

  if (.not.allocated(xk)) allocate(xk(0:ik+1))
  xk(0:ik+1)=xk_(0:ik+1)

  if (.not.allocated(lnp)) allocate(lnp(ns,0:iw+1))
  lnp(1:ns,0:iw+1)=lnp_(1:ns,0:iw+1)

  if (.not.allocated(dlnp)) allocate(dlnp(ns))
  dlnp(1:ns)=dlnp_(1:ns)

  if (.not.allocated(gridp)) allocate(gridp(ns,0:iw+1))
  gridp(1:ns,0:iw+1)=gridp_(1:ns,0:iw+1)

  if (.not.allocated(ekev)) allocate(ekev(ns,0:iw+1))
  ekev(1:ns,0:iw+1)=ekev_(1:ns,0:iw+1)

  if (.not.allocated(tKp)) allocate(tKp(nKp))
  tKp(1:nKp)=tKp_(1:nKp)

  if (.not.allocated(zkpa)) allocate(zkpa(nKp))
  zkpa(1:nKp)=zkpa_(1:nKp)

  if (.not.allocated(EmicHPowerInPs)) &
     allocate(EmicHPowerInPs(ir,ip))

  if (.not.allocated(EmicHePowerInPs)) &
     allocate(EmicHePowerInPs(ir,ip))

  if (.not.allocated(EmicOPowerInPs)) &
     allocate(EmicOPowerInPs(ir,ip))

  if (.not.allocated(EmicHPowerOutPs)) &
     allocate(EmicHPowerOutPs(ir,ip))

  if (.not.allocated(EmicHePowerOutPs)) &
     allocate(EmicHePowerOutPs(ir,ip))

  if (.not.allocated(EmicOPowerOutPs)) &
     allocate(EmicOPowerOutPs(ir,ip))

  ! read EMIC diffusion coef.
  call read_emic_diff_coef

  end subroutine init_ion_diff

!*****************************************************************************
!                           read_emic_diff_coef
!
! Routine reads EMIC wave diffusion coeff.
!*****************************************************************************
 subroutine read_emic_diff_coef

 integer i,j,k,m,n,iband
 real pa,daa,dap,dpp,EEo,daadpp,DD0,pkeVlog(0:iw+1)
 real,allocatable,dimension(:,:) :: wDpp,wDaa,wDap
 real,allocatable,dimension(:) :: Daa1,Dap1,Dpp1
 ! variables for EMIC
 integer :: ipEmic1,iwe,ipa1,iHband=0,iHeband=0,iReadPA=0
 real  :: LEmic1
 real,allocatable,dimension(:) :: &
     EmicKeV,EmicKeVlog,EmicOmpe1,ePA1
 logical :: DoReadEmic=.false.
 logical :: DoAllocateEmic=.false.
 character header*80,fhead*33

 ! Read H and He band EMIC for proton Daa, Dap, Dpp at L = 4.
   ipEmic=0
   EmicPower0=1.
   allocate(BEmic(3)) ! H,He,and O band 
   if (iEmicDiff>0) then
      do iband=1,2   ! 1=H-band, 2=He-band
         if (iband==1) write(*,*) ' ... reading H band EMIC waves Daa, Dap, Dpp'
         if (iband==2) write(*,*) ' ... reading He band EMIC waves Daa, Dap, Dpp'
         do n=1,ijs 
            pkeVlog(:)=log10(ekev(n,:))     ! electron energy grid
            if (iband==1.and.js(n)==1) then
               open(unit=48,file='D_EMIC_H_p.dat',status='old')
               write(*,*) ' ... reading D_EMIC_H_p.dat'
               DoReadEmic=.true.
            else if (iband==1.and.js(n)==2) then
               open(unit=48,file='D_EMIC_H_o.dat',status='old')
               write(*,*) ' ... reading D_EMIC_H_o.dat'
               DoReadEmic=.true.
            else if (iband==2.and.js(n)==1) then
               open(unit=48,file='D_EMIC_He_p.dat',status='old')
               write(*,*) ' ... reading D_EMIC_He_p.dat'
               DoReadEmic=.true.
            else if (iband==2.and.js(n)==2) then
               open(unit=48,file='D_EMIC_He_o.dat',status='old')
               write(*,*) ' ... reading D_EMIC_He_o.dat'
               DoReadEmic=.true.
            endif 

            if (DoReadEmic) then
               read(48,*) header
               read(48,*) LEmic1
               read(48,*) ipEmic1,iwe,ipa1  ! fpe/fce, energy, and a0 bins
               ! allocate Emic Daa, Dap, Dpp
               if (iband==1) then
                  if (.not.allocated(EmicHeDaa)) &
                     allocate (EmicHeDaa(ijs,ipEmic1,0:iw+1,ipa1),&
                               EmicHeDap(ijs,ipEmic1,0:iw+1,ipa1), &
                               EmicHeDpp(ijs,ipEmic1,0:iw+1,ipa1))
                  if (.not.allocated(EmicHDaa)) &
                     allocate (EmicHDaa(ijs,ipEmic1,0:iw+1,ipa1),&
                               EmicHDap(ijs,ipEmic1,0:iw+1,ipa1), &
                               EmicHDpp(ijs,ipEmic1,0:iw+1,ipa1))
                     BEmic(1)=xme/(re_m*LEmic1)**3
                     BEmic(2)=xme/(re_m*LEmic1)**3
               endif
               allocate (ePA1(ipa1),EmicOmpe1(ipEmic1),EmicKeV(iwe),EmicKeVlog(iwe))
               allocate (wDaa(iwe,ipa1),wDap(iwe,ipa1),wDpp(iwe,ipa1))
               allocate (Daa1(iwe),Dap1(iwe),Dpp1(iwe))
               read(48,*) EmicOmpe1
               read(48,*) EmicKeV
               EmickeVlog(:)=log10(EmicKeV(:))
               do j=1,ipEmic1
                  do k=1,iwe
                     read(48,*) header
                     read(48,*) header
                     do m=1,ipa1
                        read(48,*) pa,Daa,Dap,Dpp
                        wDaa(k,m)=Daa
                        wDap(k,m)=Dap
                        wDpp(k,m)=Dpp
                        if (iReadPA==0) ePA1(m)=pa
                     enddo ! end of m
                     iReadPA=1
                  enddo    ! end of k, iwe
                  ! map coefficients to ekev grid
                  do m=1,ipa1
                     Daa1(:)=wDaa(:,m)
                     Dap1(:)=wDap(:,m)
                     Dpp1(:)=wDpp(:,m)
                     do k=0,iw+1
                        if      (iband==1) then
                           EmicHDaa(n,j,k,m)=0.
                           EmicHDap(n,j,k,m)=0.
                           EmicHDpp(n,j,k,m)=0.
                        else if (iband==2) then
                           EmicHeDaa(n,j,k,m)=0.
                           EmicHeDap(n,j,k,m)=0.
                           EmicHeDpp(n,j,k,m)=0.
                        endif
                        if(pkeVlog(k).ge.EmicKeVlog(1).and.&
                           pkeVlog(k).le.EmickeVlog(iwe))then
                           call lintp(EmicKevlog,Daa1,iwe,pkeVlog(k),DD0)
                           EmicHeDaa(n,j,k,m)=DD0
                           if (iband==1) EmicHDaa(n,j,k,m)=DD0
                           if (iband==2) EmicHeDaa(n,j,k,m)=DD0
                           call lintp(EmicKevlog,Dpp1,iwe,pkeVlog(k),DD0)
                           EmicHeDpp(n,j,k,m)=DD0
                           if (iband==1) EmicHDpp(n,j,k,m)=DD0
                           if (iband==2) EmicHeDpp(n,j,k,m)=DD0
                           if (iband==1) DaaDpp=sqrt(EmicHDaa(n,j,k,m)*EmicHDpp(n,j,k,m))
                           if (iband==2) DaaDpp=sqrt(EmicHeDaa(n,j,k,m)*EmicHeDpp(n,j,k,m))
                           call lintp(EmicKevlog,Dap1,iwe,pkeVlog(k),DD0)
                           if (DaaDpp.lt.abs(DD0)) DD0=DD0*DaaDpp/abs(DD0)
                           EmicHeDap(n,j,k,m)=DD0
                           if (iband==1) EmicHDap(n,j,k,m)=DD0
                           if (iband==2) EmicHeDap(n,j,k,m)=DD0
                        endif
                     enddo ! end if k,iw
                  enddo    ! end of mapping coefficients to ekev grid
               enddo       ! end of j
               if (ipEmic==0.and.iband==1) then
                  ipEmic=ipEmic1
                  iePA=ipa1
                  allocate(EmicOmpe(ipEmic))
                  allocate(ePA(iepa))
                  EmicOmpe(:)=EmicOmpe1(:)
                  ePA(:)=ePA1(:)
               else
                  if (ipEmic/=ipEmic1) write(*,*) ' **WARNING: ipEmic /= ipEmic1'      
                  if (iepa/=ipa1) write(*,*) ' **WARNING: iepa /= ipa1'      
               endif       
               close(48)
               deallocate (EmicOmpe1,ePA1,EmicKeV,EmicKeVlog)
               deallocate (wDaa,wDap,wDpp,Daa1,Dap1,Dpp1)
               DoReadEmic=.false.
            endif ! end of if (DoReadEmic) then
         enddo    ! end of do n=1,ijs 
      enddo       ! end of do iband=1,2   ! 1=H-band, 2=He-band
   endif          ! end of if (iEmicDiff>0) then

  end subroutine read_emic_diff_coef


!****************************************************************************
!                             diffuse_UV_ions  
!
!  This subroutine solves diffusion in velocity space 
!   in (U,V)=(lnp,lnK) coordinates.
!****************************************************************************
  recursive subroutine diffuse_UV_ions(n,i,j,Gjac0,t,ro1,bo1,y0,tya0,rppa1,ompe1,f2)

  integer,intent(in) :: n,i,j  
  real,   intent(in) :: t,ro1,bo1,rppa1,ompe1
  real,intent(in),contiguous :: Gjac0(:,:),y0(:),tya0(:)
  real,intent(inout),contiguous :: f2(:,:)   ! psd * jacobian ~ # of particle

  integer m,k,ipEmic1,iexit,mpa
  ! note, contiguous arrays start from 1 and Gjac,y1,tya1 start from 0
  real :: Gjac(0:iw+1,0:ik+1),y1(0:ik+1),tya1(0:ik+1)
  real EmicHFact,&   ! H band EMIC wave amplitude factor
       EmicHeFact,&   ! He band EMIC wave amplitude factor
       EmicOFact,&   ! O band EMIC wave amplitude factor
       cosa,pp3, &
       f2d(iw,ik),BLbo(3)
  real :: bwH,bwHe,bwO
  real rbo,ao,DDm,DDp,y_2,kak,dKda,dtU2,dtV2
  real Daa,Dap,Dpp,DUU(0:iw+1,0:ik+1),DVV(0:iw+1,0:ik+1),DUV(0:iw+1,0:ik+1)
  real :: wDaa(0:iw+1,iePA),wDap(0:iw+1,iePA),wDpp(0:iw+1,iePA)
  real :: wPA(iePA),DD(0:max(iw,ik)+1),um(max(iw,ik)),up(max(iw,ik))
  real :: a1d(max(iw,ik)),b1d(max(iw,ik)),c1d(max(iw,ik)),f0(0:max(iw,ik)+1),fr(max(iw,ik)),fnew(max(iw,ik))

  real :: iDaa=1,&
          iDpp=1

 ! xlam=0.5              ! implicitness in solving diffusion equation
 ! alam=1.-xlam          ! xlam=1: fully implicit; xlam=0.5: Crank-Nicolson
 mpa=iePA
 wPA(1:mpa)=ePA(1:mpa)

 ! Gjac,y1,tya1
 Gjac(0:iw+1,0:ik+1)=Gjac0(1:iw+2,1:ik+2)
 y1(0:ik+1)=y0(1:ik+2)
 tya1(0:ik+1)=tya0(1:ik+2)

 dtU2=dt/dlnp(n)/dlnp(n)
 dtV2=dt/dlnK/dlnK
 pp3=(gridp(n,1)/gridp(n,0))**3
 if (iEmicDiff/=0) then
    ! initialize Bw arrays
    EmicHpowerInPs(i,j)=0.
    EmicHepowerInPs(i,j)=0.
    EmicHpowerOutPs(i,j)=0.
    EmicHepowerOutPs(i,j)=0.
    BwH =0.
    BwHe=0.
    BwO =0.

    call EmicPower(t,ro1,bwH,bwHe,bwO)
    BLbo(:)=BEmic(:)/bo1

    iexit=1
    ! outside plasmasphere
    if (ro1.gt.rppa1) then
       ! no EMIC wave
       if (BwH==0..and.BwHe==0.) goto 9999
       EmicHpowerOutPs(i,j)=BwH
       EmicHepowerOutPs(i,j)=BwHe
       EmicOpowerOutPs(i,j)=BwO
       call locate1(EmicOmpe,ipEmic,ompe1,ipEmic1)
       if (ipEmic1==0) goto 9999                   ! beyond fpe/fce range
       EmicHFact=BwH/EmicPower0*BLbo(1)
       EmicHeFact=BwHe/EmicPower0*BLbo(2)
       EmicOFact=BwO/EmicPower0*BLbo(3)
       if (ipEmic1.lt.ipEmic) then
          if ((EmicOmpe(ipEmic1+1)-Ompe1)<(Ompe1-EmicOmpe(ipEmic1))) &
             ipEmic1=ipEmic1+1
       else
          ipEmic1=ipEmic
       endif
       iexit=0
       do k=0,iw+1
          wDaa(k,:)=EmicHDaa(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODaa(n,ipEmic1,k,:)*EmicOFact
          wDap(k,:)=EmicHDap(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDap(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDap(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODap(n,ipEmic1,k,:)*EmicOFact
          wDpp(k,:)=EmicHDpp(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODpp(n,ipEmic1,k,:)*EmicOFact
       enddo
    endif
    ! insdie plasmasphere
    if (ro1.le.rppa1) then
       ! no EMIC wave
       if (BwH==0..and.BwHe==0.) goto 9999
       EmicHpowerInPs(i,j)=BwH
       EmicHepowerInPs(i,j)=BwHe
       EmicOpowerInPs(i,j)=BwO
       call locate1(EmicOmpe,ipEmic,ompe1,ipEmic1)
       if (ipEmic1==0) goto 9999                   ! beyond fpe/fce range
       EmicHFact=BwH/EmicPower0*BLbo(1)
       EmicHeFact=BwHe/EmicPower0*BLbo(2)
       EmicOFact=BwO/EmicPower0*BLbo(3)
       if (ipEmic1.lt.ipEmic) then
          if ((EmicOmpe(ipEmic1+1)-Ompe1)<(Ompe1-EmicOmpe(ipEmic1))) &
             ipEmic1=ipEmic1+1
       endif
       iexit=0
       do k=0,iw+1
          wDaa(k,:)=EmicHDaa(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODaa(n,ipEmic1,k,:)*EmicOFact
          wDap(k,:)=EmicHDap(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDap(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDap(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODap(n,ipEmic1,k,:)*EmicOFact
          wDpp(k,:)=EmicHDpp(n,ipEmic1,k,:)*EmicHFact&
                   +EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact
                   ! when O band Emic waves are included
                   !+EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact&
                   !+EmicODpp(n,ipEmic1,k,:)*EmicOFact
       enddo
    endif

    if (iexit.eq.1) goto 9999

    ! Setup DUU, DVV, DUV
    rbo=ro1*re_m*sqrt(bo1)
    do m=0,ik+1
       y_2=y1(m)**2
       cosa=sqrt(1.-y_2)
       dKda=2.*rbo*cosa*tya1(m)/y_2
       kak=dKda/xk(m)
       ao=asin(y1(m))*180./pi

       if (iEmicDiff/=0) then
          if (ao.lt.wPA(1)) ao=wPA(1)
          if (ao.gt.wPA(mpa)) ao=wPA(mpa)
       endif
       do k=0,iw+1
          if (iEmicDiff/=0.) then
             call lintp(wPA,wDaa(k,:),mpa,ao,Daa)
             call lintp(wPA,wDap(k,:),mpa,ao,Dap)
             call lintp(wPA,wDpp(k,:),mpa,ao,Dpp)
          else
             Daa=0.
             Dap=0.
             Dpp=0.
          endif
          DUU(k,m)=Dpp
          DUV(k,m)=kak*Dap
          DVV(k,m)=kak*kak*Daa
       enddo
    enddo


    ! Setup 2D psd in (V,U) and Gjac
    do m=1,ik
       do k=1,iw
          f2d(k,m)=f2(k,m)/Gjac(k,m)     ! psd in cimi grid
       enddo
    enddo

    ! do diffusion in U(lnp)
    if (iDpp.eq.1) then
       do m=1,ik
          do k=0,iw+1
             DD(k)=Gjac(k,m)*DUU(k,m)
          enddo
          do k=1,iw
             DDm=0.5*(DD(k)+DD(k-1))
             DDp=0.5*(DD(k)+DD(k+1))
             um(k)=dtU2*DDm/Gjac(k,m)
             up(k)=dtU2*DDp/Gjac(k,m)
          enddo
          do k=1,iw
             a1d(k)=-um(k)
             b1d(k)=1.+(um(k)+up(k))
             c1d(k)=-up(k)
          enddo
          ! start diffusion in U 
          f0(1:iw)=f2d(1:iw,m)
          f0(0)=f0(1)             ! f0 is psd, psd(0)=psd(1) for e-
          f0(iw+1)=f0(iw)/pp3     ! f0 is psd, f2(iw+1)=f2(iw) for e-
          fr(1:iw)=f0(1:iw)
          fr(1)=fr(1)+um(1)*f0(0)
          fr(iw)=fr(iw)+up(iw)*f0(iw+1)
          call tridiagonal(a1d,b1d,c1d,fr,iw,fnew)
          f2d(:,m)=fnew(:)
       enddo
    endif   ! end of if (iDpp.eq.1)

    ! do diffusion in V(lnK)
    if (iDaa.eq.1) then
       do k=1,iw
          do m=0,ik+1
             DD(m)=Gjac(k,m)*DVV(k,m)
          enddo
          do m=1,ik
             DDm=0.5*(DD(m)+DD(m-1))
             DDp=0.5*(DD(m)+DD(m+1))
             um(m)=dtV2*DDm/Gjac(k,m)
             up(m)=dtV2*DDp/Gjac(k,m)
          enddo
          do m=1,ik
             a1d(m)=-um(m)
             b1d(m)=1.+(um(m)+up(m))
             c1d(m)=-up(m)
          enddo
          ! start diffusion in V
          f0(1:ik)=f2d(k,1:ik)
          fr(1:ik)=f0(1:ik)
          ! boundary condition
          b1d(1)=a1d(1)+b1d(1)    ! flat distribution near 90 deg pitch angle
          b1d(ik)=b1d(ik)+c1d(ik) ! flat distribution near zero deg pitch angle 
          call tridiagonal(a1d,b1d,c1d,fr,ik,fnew)
          f2d(k,:)=fnew(:)
       enddo
    endif    ! end of if (iDaa.eq.1)

    ! Update f2
    do m=1,ik
       do k=1,iw
          f2(k,m)=f2d(k,m)*Gjac(k,m)
       enddo
    enddo
9999   continue

 endif

end subroutine diffuse_UV_ions 
 
 
!*****************************************************************************'
!                          EmicPower  
!
! This subroutine calculates the EMIC intensity 
!  base on the CRRES data (Kersten et al. 2014)
!
! EMIC wave amplitude is categorized by Kp index
!
! L = 3.5 - 7, every 0.5
!*****************************************************************************
  recursive subroutine EmicPower(t,ro1,BwH,BwHe,BwO)

  real,intent(in) :: t,ro1
  real,intent(out) :: BwH,BwHe,BwO   ! wave amplitude in nT
  real Kp
  integer iKp,jKp,i,j
  
  !call locate1(tKp(1:nKp),nKp,t,iKp)
  call lintp(tKp(1:nKp),zkpa(1:nKp),nKp,t,Kp)

  if      (Kp<2. ) then
     jKp=1
  else if (Kp>=4.) then
     jKp=3
  else!if (zkpa(iKp)>=2..and.zkpa(iKp)<4.) then
     jKp=2
  endif

  ! interpolate EMIC wave intensity
  call lintp(L_CRRES,EMIC_CRRES_H(:,jKp) ,nCRRES,ro1,BwH )
  call lintp(L_CRRES,EMIC_CRRES_He(:,jKp),nCRRES,ro1,BwHe)
  BwO=0.

  end subroutine EmicPower


! ************************************************************************
!                          tridiagonal
!
!  This subroutine solves inversion of tridiagonal matrix
!   through Thomas's tridiagonal matrix algorithm.
!
!  This subroutine is copied from cimi.f90
!
!  b1*fi-1 + b2*fi + b3*fi+1 = b4,  i= 1,2,3,..,n
!
!  Inputs : b1,b2,b3,b4,n
!  Outputs : f
! ************************************************************************
     recursive subroutine tridiagonal(b1,b2,b3,b4,n,f)
     implicit none

     integer,intent(in) :: n
     real,intent(in),contiguous :: b1(:),b2(:),b3(:),b4(:)
     real,intent(out) :: f(n)
     real c0,c1(n),c2(n)
     integer i

     c1(1)=b3(1)/b2(1)
     c2(1)=b4(1)/b2(1)
     do i=2,n
        c0=b2(i)-b1(i)*c1(i-1)
        c1(i)=b3(i)/c0
        c2(i)=(b4(i)-b1(i)*c2(i-1))/c0
     enddo

     f(n)=c2(n)
     do i=n-1,1,-1
        f(i)=c2(i)-c1(i)*f(i+1)
     enddo

     end subroutine tridiagonal

!-----------------------------------------------------------------------
      recursive subroutine lintp(xx,yy,n,x,y)

!  This subroutine is copied from cimi.f90
!-----------------------------------------------------------------------
!  Routine does 1-D interpolation.  xx must be increasing or decreasing
!  monotonically. If x is beyound xx, it will be forced inside xx range.

      
      integer,intent(in) :: n
      real,intent(in) :: x
      !real,intent(in),contiguous :: xx(:),yy(:)
      real,intent(in) :: xx(n),yy(n)
      real,intent(out) :: y
      integer i,j,jl,ju,jm
      real x1,minxx,maxxx,d

!  Make sure xx is increasing or decreasing monotonically
      do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            write(*,*) ' lintp: xx is not increasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            write(*,*) ' lintp: xx is not decreasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
      enddo

!  Make sure x is inside xx range
      minxx=minval(xx)
      maxxx=maxval(xx)
      x1=x
      if (x1.lt.minxx) x1=minxx 
      if (x1.gt.maxxx) x1=maxxx 

!    initialize lower and upper values
!
      jl=1
      ju=n
!
!    if not dne compute a midpoint
!
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
!
!    now replace lower or upper limit
!
        if((xx(n).gt.xx(1)).eqv.(x1.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
!
!    try again
!
      go to 10
      endif
!
!    this is j
!
      j=jl      ! if x.le.xx(1) then j=1
!                 if x.gt.xx(j).and.x.le.xx(j+1) then j=j
!                 if x.gt.xx(n) then j=n-1
      d=xx(j+1)-xx(j)
      y=(yy(j)*(xx(j+1)-x1)+yy(j+1)*(x1-xx(j)))/d

      end subroutine lintp


!--------------------------------------------------------------------------
      recursive subroutine locate1(xx,n,x,j)

!  This subroutine is copied from cimi.f90
!--------------------------------------------------------------------------
!  Routine return a value of j such that x is between xx(j) and xx(j+1).
!  xx must be increasing or decreasing monotonically. If not, the locate will
!  stop at the turning point.
!  If xx is increasing:
!     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
!     If x < xx(1), j=0  and if x > xx(n), j=n
!  If xx is decreasing:
!     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
!     If x > xx(1), j=0  and if x < xx(n), j=n

      integer,intent(in) :: n
      real,intent(in) :: x
      real,intent(in),contiguous :: xx(:)
      integer,intent(out) :: j

      integer nn,i,jl,ju,jm

!  Make sure xx is increasing or decreasing monotonically
      nn=n
      monoCheck: do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
      enddo monoCheck
      if (nn.ne.n) then
         write(*,*)'locate1: xx is not increasing or decreasing monotonically'
         write(*,*)'n,x ',n,x
         write(*,*)'xx ',xx
         stop
      endif

      jl=0
      ju=nn+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(nn).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl

      end subroutine locate1


!****************************************************************************
!                        unit_test_diffuse_ions  
!
!  This soubroutine test ion diffusion calculation
!****************************************************************************
 subroutine unit_test_diffuse_ions
 integer,parameter :: ns_=1,ir_=1,ip_=1,iw_=48,ik_=40,nKp_=4,&
                      ijs_=1

 integer :: js_(ns_)

 real :: tKp_(nKp_),zkpa_(nKp_),dt_=3600. ,&
         xk_(0:ik_+1),dlnK_,lnp_(ns_,0:iw_+1),&
         dlnp_(ns_),gridp_(ns_,0:iw_+1),&
         ekev_(ns_,0:iw_+1),&
         xme_=7.72e15,&                 ! dipole moment in T m3
         GridEmin=0.1, GridEmax=500.,&  ! min and max energy in keV
         Eo=938272.,&                   ! proton rest mass in keV
         Gridpmin,Gridpmax,dlnGridp,&   ! min, max, and log momentum grids
         GridKmin=58.,dGridK=0.267,&    ! T^0.5 m
         lnGridK,&                      ! log (T^0.5 m)
         cCspeed=2.99792458e8,&         ! speed of light, m/s              
         cCharge=1.60217663e-19,&       ! elementary charge, Coulomb
         cKeVtoEv=1000.,&               ! keV to eV         
         bm(0:ik_+1),&                  ! B at the mirror points corresponding to xk
         y(0:ik_+1),&                   ! sin(a0) corresponding to xk
         tya(0:ik_+1)                   ! T(a0) corresponding to xk
         

  real :: ro1=4.,bo1,rppa1=6., ompe1=20.,&
          t=2.,&
          Gjac(0:iw_+1,0:ik_+1),&   ! Jacobian 
          f2(iw_,ik_),&             ! psd * Jacobian, # particle in dpdK
          E0=10.,&  ! 10 keV
          psd0=1.e1 ! unitless


  integer :: k,m,iKp,n=1,i=1,j=1,&
             iEner=27,&           ! ~ 10keV     
             iiter,niter=9     

 !
 js_(1:ns_)=1

 ! momentum grids
 Gridpmin=sqrt(GridEmin*(GridEmin+2.*Eo))*cKeVtoEv*cCharge/cCspeed
 Gridpmax=sqrt(GridEmax*(GridEmax+2.*Eo))*cKeVtoEv*cCharge/cCspeed
 dlnGridP=log(GridPmax/GridPmin)/float(iw_)
 do k=0,iw_+1
    lnp_(n,k)=log(gridpmin)+(float(k)-0.5)*dlnGridp
    gridp_(n,k)=exp(lnp_(n,k))
    ekev_(n,k)=sqrt(Eo**2 + (gridp_(n,k)*cCspeed/cCharge/cKeVtoEv)**2)-Eo
 enddo        
! TEST
!write(*,*) 'gridp '
!write(*,'(1p10E9.2)') gridp_(n,:)

 ! Grids in K (T^0.5 m)
 do m=0,ik_+1
    lnGridK=log(GridKmin)+float(m-1)*dGridK 
    xk_(m)=exp(lnGridK) 
 enddo

 ! dipole bm
 bm(0:ik_+1)=&
 (/1.88496E-06, 1.88668E-06, 1.88892E-06, 1.89185E-06, 1.89568E-06,& 
   1.90067E-06, 1.90720E-06, 1.91572E-06, 1.92685E-06, 1.94139E-06,&
   1.95912E-06, 1.98203E-06, 2.01196E-06, 2.05104E-06, 2.10209E-06,&
   2.17001E-06, 2.25939E-06, 2.37613E-06, 2.53211E-06, 2.73828E-06,&
   3.01398E-06, 3.38653E-06, 3.89192E-06, 4.56922E-06, 5.52919E-06,& 
   6.86238E-06, 8.72844E-06, 1.14093E-05, 1.54000E-05, 2.12938E-05,&
   2.99599E-05, 4.33165E-05, 6.35242E-05, 9.51649E-05, 1.50478E-04,& 
   2.28428E-04, 3.67746E-04, 5.82127E-04, 9.38879E-04, 1.54061E-03,&
   2.48867E-03, 4.19828E-03/)

  ! sin of equatorial pitch-angle
  y(0:ik_+1)=&
 (/0.99851,0.99806,0.99746,0.99669,0.99568,0.99438,0.99267,0.99046,0.98760,0.98389,&
   0.97943,0.97375,0.96648,0.95723,0.94554,0.93062,0.91203,0.88934,0.86152,0.82845,&
   0.78965,0.74495,0.69490,0.64133,0.58301,0.52332,0.46402,0.40586,0.34934,0.29708,&
   0.25046,0.20829,0.17200,0.14053,0.11176,0.09070,0.07149,0.05682,0.04474,0.03493,&
   0.02748,0.02116/)

  ! helical path (gyro-bounce path)
  tya(0:ik_+1)=& 
 (/6.53811E-01,6.53811E-01,6.53811E-01,6.53811E-01,6.53811E-01,&
   6.53811E-01,6.53811E-01,6.53811E-01,6.53811E-01,6.53811E-01,&
   6.54843E-01,6.56388E-01,6.58407E-01,6.61044E-01,6.64487E-01,&
   6.71778E-01,6.82816E-01,6.97232E-01,7.13606E-01,7.33270E-01,&
   7.57429E-01,7.86025E-01,8.15947E-01,8.48223E-01,8.80728E-01,&
   9.16570E-01,9.54753E-01,9.92423E-01,1.02782E+00,1.06469E+00,&
   1.10399E+00,1.14404E+00,1.16811E+00,1.19352E+00,1.21302E+00,&
   1.23595E+00,1.25446E+00,1.27410E+00,1.28901E+00,1.30294E+00,&
   1.31452E+00,1.32450E+00/)

  ! Jacobian
  do m=0,ik_+1
     do k=0,iw_+1
        Gjac(k,m)=xk_(m)/bm(m)**1.5*gridp_(n,m)**3
     enddo
  enddo
! TEST
!write(*,'(1p10E9.2)') Gjac(:,ik)


  ! time in sec
  do iKp=1,nKp_
     tKp_(iKp)=float(iKp)
     zkpa_(iKp)=5.
  enddo

  bo1=bm(1)*y(1)**2
  dlnp_(1:n)=dlnGridp
  dlnK_=dGridK

! TEST
!  write(*,'(1p10E9.2)') ekev_(n,1:iw_)

  call init_ion_diff(ns_,ir_,ip_,iw_,ik_,ijs_,js_,&
                     xk_,dlnK_,lnp_,dlnp_,gridp_,ekev_,&
                     nKp_,tKp_,zkpa_,xme_,dt_)

  ! phase space density
  do k=1,iw_
     do m=1,ik_
        f2(k,m)=psd0*ekev_(n,k)/E0*exp(-ekev_(n,k)/E0)*y(m)*Gjac(k,m)
     enddo
  enddo

  write(*,'(a,f9.2)') 'energy (keV) :',ekev_(n,iEner)
  write(*,*) 'psd(iEner,:)'
  write(*,*) 'psd before'
  write(*,'(1p10E9.2)') f2(iEner,1:ik_)/Gjac(iEner,1:ik_)


! TEST
!EmicHDaa(n,:,:,:)=1.e-3
!EmicHeDaa(n,:,:,:)=1.e-3
!EmicHDpp(n,:,:,:)=0.
!EmicHeDpp(n,:,:,:)=0.

  ! calculate ion diffusion in UV space
  call diffuse_UV_ions(n,i,j,Gjac,t,ro1,bo1,y,tya,rppa1,ompe1,f2)

  write(*,'(a,f9.1)') 'psd after, t= ',dt 
  write(*,'(1p10E9.2)') f2(iEner,1:ik_)/Gjac(iEner,1:ik_)

  ! calculate ion diffusion in UV space
  do iiter=1,niter
     call diffuse_UV_ions(n,i,j,Gjac,t,ro1,bo1,y,tya,rppa1,ompe1,f2)
  enddo

  write(*,'(a,f9.1)') 'psd after, t= ',dt*10
  write(*,'(1p10E9.2)') f2(iEner,1:ik_)/Gjac(iEner,1:ik_)


  write(*,'(a,10f8.4)')'H EMIC power (inside/outside plasmasphere): ',&
                        EmicHPowerInPs(i,j),EmicHPowerOutPs(i,j)
  write(*,'(a,10f8.4)')'He EMIC power(inside/outside plasmasphere): ',&
                        EmicHePowerInPs(i,j),EmicHePowerOutPs(i,j)
  write(*,'(a)')' Daa H (rad^2/s): '
  write(*,'(1p10E8.1)') EmicHDaa(n,ipEmic,:,ik_)
  write(*,'(a)')' Daa He (rad^2/s): '
  write(*,'(1p10E8.1)') EmicHeDaa(n,ipEmic,:,ik_)
  write(*,'(a)')' pitch-angle (deg): '
  write(*,'(10f8.2)') asin(y)*180./acos(-1.)

  end subroutine unit_test_diffuse_ions

          
  end module           


