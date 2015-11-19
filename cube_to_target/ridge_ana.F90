#define SPEEDUPLOOP
module ridge_ana

use rotation 
USE reconstruct
use shr_kind_mod, only: r8 => shr_kind_r8

IMPLICIT NONE
private


public find_ridges
public remapridge2target
public remapridge2tiles
public paintridgeoncube

public anglx_target,aniso_target,mxdis_target,hwdth_target
public mxvrx_target,mxvry_target,bsvar_target,wghts_target,riseq_target
public ang22_target,anixy_target,clngt_target,cwght_target,count_target
public nsubr,grid_length_scale,fallq_target

  REAL, allocatable  :: MXVRX(:,:,:),MXDIS(:,:,:),MNSLP(:,:,:),ANGLX(:,:,:),ANISO(:,:,:),XS(:),YS(:)
  REAL, allocatable  :: XSPK(:,:,:),YSPK(:,:,:),MXDS0(:,:,:),MXDS1(:,:,:),SFT0(:,:,:),SFT1(:,:,:)
  REAL, allocatable  :: PKHTS(:,:,:),VLDPS(:,:,:),RWPKS(:,:,:),RWVLS(:,:,:),ANGLL(:,:,:)
  REAL, allocatable  :: BSVAR(:,:,:),HWDTH(:,:,:),NPKS(:,:,:),NVLS(:,:,:),MXVRY(:,:,:)
  REAL, allocatable  :: RISEQ(:,:,:),FALLQ(:,:,:),ANIXY(:,:,:),MXDSP(:,:,:),CLNGTH(:,:,:),MXDS2(:,:,:)

    REAL(KIND=dbl_kind), allocatable ::  ALP0(:,:),BET0(:,:),LAT0(:,:,:),LON0(:,:,:)
    REAL(KIND=dbl_kind), allocatable ::  ALP1(:,:),BET1(:,:),LAT1(:,:,:),LON1(:,:,:)

  real(r8), allocatable, dimension(:,:) :: anglx_target,aniso_target,mxdis_target,hwdth_target
  real(r8), allocatable, dimension(:,:) :: mxvrx_target,mxvry_target,bsvar_target,wghts_target 
  real(r8), allocatable, dimension(:,:) :: ang22_target,anixy_target,clngt_target,cwght_target
  real(r8), allocatable, dimension(:,:) :: count_target,riseq_target,fallq_target
  !!,rwpks_target

    INTEGER (KIND=int_kind),allocatable :: UQRID(:,:,:) 

  real(r8), allocatable, dimension(:,:) :: anglx_tiles,aniso_tiles,mxdis_tiles,hwdth_tiles
  real(r8), allocatable, dimension(:,:) :: mxvrx_tiles,mxvry_tiles,bsvar_tiles,wghts_tiles 
  real(r8), allocatable, dimension(:,:) :: ang22_tiles,agnpk_tiles,bgnpk_tiles
  integer,  allocatable, dimension(:,:) :: uqrid_tiles
  integer,  allocatable, dimension(:)   :: numbr_tiles,error_tiles



    REAL, allocatable ::  wt1p(:,:)

    REAL, allocatable :: agnom(:),bgnom(:)

    REAL(r8) :: grid_length_scale

    REAL (KIND=dbl_kind), PARAMETER :: pi        = 3.14159265358979323846264338327
    REAL (KIND=dbl_kind), PARAMETER :: earth_radius        = 6371.0

    integer, parameter             ::   NANG=16
    integer, parameter             ::   NSUBR = NANG
    integer, parameter             ::   NTILES = 250

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_ridges ( terr_dev, terr_raw, ncube, nhalo, nsb, nsw )
!------------------------------------------------
!  INPUTS.
!      NSW = size of window used for ridge analysis
!      




    REAL (KIND=dbl_kind), &
            DIMENSION(ncube,ncube,6), INTENT(IN) :: terr_dev
    REAL (KIND=dbl_kind), &
            DIMENSION(ncube,ncube,6), INTENT(IN) :: terr_raw

       INTEGER (KIND=int_kind), INTENT(IN) :: ncube, nhalo, nsb, nsw
       INTEGER (KIND=int_kind) :: i,j,np,ncube_halo,ipanel,N,norx,nory,ip



    REAL (KIND=dbl_kind),                                            &
         DIMENSION(1-nhalo:ncube+nhalo, 1-nhalo:ncube+nhalo, 6) :: terr_halo
    REAL (KIND=dbl_kind),                                            &
         DIMENSION(1-nhalo:ncube+nhalo, 1-nhalo:ncube+nhalo, 6) :: terr_dev_halo
    REAL ,                                            &
         DIMENSION(1-nhalo:ncube+nhalo, 1-nhalo:ncube+nhalo, 6) :: terr_halo_r4
    REAL ,                                            &
         DIMENSION(1-nhalo:ncube+nhalo, 1-nhalo:ncube+nhalo, 6) :: terr_dev_halo_r4

    REAL  ,                                                          &
         DIMENSION(1-nhalo:ncube+nhalo )                          :: xv,yv,alph,beta


    REAL(KIND=dbl_kind)  :: lon_r8, lat_r8, cosll, dx, dy, dcube2, ampfsm,dbet,dalp,diss,diss00

    CHARACTER(len=1024) :: ofile$,ve$

!---------------------------------------------------------------------------------------
!  B E G I N   C A L C U L A T I O N S
!---------------------------------------------------------------------------------------

    DO np = 1, 6
     CALL CubedSphereFillHalo_Linear_extended(terr_raw, terr_halo(:,:,np), np, ncube+1,nhalo)  
     CALL CubedSphereFillHalo_Linear_extended(terr_dev, terr_dev_halo(:,:,np), np, ncube+1,nhalo)  
    END DO


   
    DO i=1-nhalo,ncube+nhalo
       xv(i)=1.*i   !  xv,yv are 'SW' corners
       yv(i)=1.*i   
       !xv(i)=1.*i - 0.5   !  xv,yv are cell centers
       !yv(i)=1.*i - 0.5
    END DO
    DO i=1-nhalo,ncube+nhalo
       alph(i) =  ( xv(i) - 0.5 - ncube/2 )*(pi/2.)/ncube
       beta(i) =  ( yv(i) - 0.5 - ncube/2 )*(pi/2.)/ncube
       !alph(i) =  ( xv(i) - ncube/2 )*(pi/2.)/ncube
       !beta(i) =  ( yv(i) - ncube/2 )*(pi/2.)/ncube
    END DO

    grid_length_scale = ( alph(1)-alph(0) )*earth_radius


    write(*,*) "ALPHA(1:NCUBE) ",alph(1)*180./pi,alph(ncube)*180./pi

    write(*,*) " GRID scale ",grid_length_scale

    allocate( agnom(ncube) )
    allocate( bgnom(ncube) )
    agnom(1:ncube) = alph(1:ncube)
    bgnom(1:ncube) = beta(1:ncube)

      terr_dev_halo_r4 = terr_dev_halo
      terr_halo_r4     = terr_halo


      ncube_halo = size( terr_halo_r4, 1)




! Allocate arrays to hold results on reduced Cubed-sphere grid
!--------------------------------------------------------------

  N = ncube_halo
  allocate( xs(N/nsb) )
  allocate( ys(N/nsb) )
  ys=-999.
  xs=-999.

  allocate( bsvar( N/nsb , N/nsb , 6) )
  bsvar=0. !-999.

  allocate( mxvrx( N/nsb , N/nsb , 6) )
   mxvrx=0.
  allocate( mxvry( N/nsb , N/nsb , 6) )
   mxvry=0.
  allocate( mxdis( N/nsb , N/nsb , 6) )
   mxdis=0.
  allocate( mxds2( N/nsb , N/nsb , 6) )
   mxds2=0.
  allocate( anglx( N/nsb , N/nsb , 6) )
   anglx=-999.
  allocate( angll( N/nsb , N/nsb , 6) )
   angll=-999.
  allocate( aniso( N/nsb , N/nsb , 6) )
   aniso=0.
  allocate( anixy( N/nsb , N/nsb , 6) )
   anixy=0.
  allocate( mnslp( N/nsb , N/nsb , 6) )
   mnslp=0.
  allocate( xspk( N/nsb , N/nsb , 6) )
   xspk=0.
  allocate( yspk( N/nsb , N/nsb , 6) )
   yspk=0.
  allocate( mxdsp( N/nsb , N/nsb , 6) )
   mxdsp=0.
  allocate( mxds0( N/nsb , N/nsb , 6) )
   mxds0=0.
  allocate( mxds1( N/nsb , N/nsb , 6) )
   mxds1=0.
  allocate( hwdth( N/nsb , N/nsb , 6) )
   hwdth=0.
  allocate( npks( N/nsb , N/nsb , 6) )
   npks=0.
  allocate( nvls( N/nsb , N/nsb , 6) )
   nvls=0.
  allocate( sft0( N/nsb , N/nsb , 6) )
   sft0=0.
  allocate( sft1( N/nsb , N/nsb , 6) )
   sft1=0.
  allocate( vldps( N/nsb , N/nsb , 6) )
   vldps=0.
  allocate( pkhts( N/nsb , N/nsb , 6) )
   pkhts=0.
  allocate( rwvls( N/nsb , N/nsb , 6) )
   rwvls=0.
  allocate( rwpks( N/nsb , N/nsb , 6) )
   rwpks=0.
  allocate( uqrid( N/nsb , N/nsb , 6) )
   uqrid=-1
  allocate( riseq( N/nsb , N/nsb , 6) )
   riseq=-1.
  allocate( fallq( N/nsb , N/nsb , 6) )
   fallq=1.
  allocate( clngth( N/nsb , N/nsb , 6) )
   clngth=1.

       norx= N/nsb !size(anglx,1)
       nory= N/nsb !size(anglx,1)


          allocate( ALP0(norx,nory))
            ALP0  = -9999.d+0
          allocate( BET0(norx,nory))
            BET0  = -9999.d+0
          allocate( LON0(norx,nory,6))
              LON0  = -9999.d+0
          allocate( LAT0(norx,nory,6))
              LAT0  = -9999.d+0
          allocate( ALP1(norx,nory))
            ALP1  = -9999.d+0
          allocate( BET1(norx,nory))
            BET1  = -9999.d+0
          allocate( LON1(norx,nory,6))
             LON1  = -9999.d+0
          allocate( LAT1(norx,nory,6))
              LAT1  = -9999.d+0




    dcube2 = 0.5 * ncube 

#ifdef DEBUGRIDGE
    do ipanel=4,4
#else
    do ipanel=1,6   !1, 6
#endif
       call ANISO_ANA( terr_dev_halo_r4(:,:,ipanel) , terr_halo_r4(:,:,ipanel) , XV , YV , ncube_halo ,NSB,NSW, ipanel )


       
       do j=1,nory
       do i=1,norx
          uqrid( i,j, ipanel ) = (ipanel-1) * norx * nory + (j-1)*norx + i
       end do
       end do

       do j=1,nory
       do i=1,norx
          alp0(i,j) = (PI/4.)*(xs(i)-dcube2)/dcube2
          bet0(i,j) = (PI/4.)*(ys(j)-dcube2)/dcube2
       end do
       end do
       do j=1,nory
       do i=1,norx
          alp1(i,j) =  alp0(i,j) + 0.01*SIN( ANGLX(i,j,ipanel)*PI/180. )   
          bet1(i,j) =  bet0(i,j) + 0.01*COS( ANGLX(i,j,ipanel)*PI/180. )   
       end do
       end do

       do j=1,nory
       do i=1,norx
          call CubedSphereRLLFromABP(alp0(i,j), bet0(i,j) , ipanel, lon_r8, lat_r8)
          lon0(i,j,ipanel) = lon_r8
          lat0(i,j,ipanel) = lat_r8
          call CubedSphereRLLFromABP(alp1(i,j), bet1(i,j) , ipanel, lon_r8, lat_r8)
          lon1(i,j,ipanel) = lon_r8
          lat1(i,j,ipanel) = lat_r8
       end do
       end do


      do j=1,nory
       do i=1,norx
          dx = COS( lat0(i,j,ipanel))*(lon1(i,j,ipanel)-lon0(i,j,ipanel))
          dy = (lat1(i,j,ipanel)-lat0(i,j,ipanel))
          if ( dx < 0.0 ) dy  = -1.*dy  ! 
          COSLL  = dy /sqrt( dx**2 + dy**2 )
          ANGLL(i,j,ipanel) = ACOS( COSLL )*180./PI
       end do
       end do



      write(*,*) " Finished with Panel = ",   ipanel

    end do

    hwdth = hwdth ! * grid_length_scale



#if 1

    ofile$ = "TEST.dat"
    OPEN (unit = 31, file= trim(ofile$) ,form="UNFORMATTED" )

    write(31) ncube_halo , grid_length_scale
    write(31) xv,yv
    write(31) terr_halo_r4
    write(31) terr_dev_halo_r4
   


write(31) size(anglx,1),size(anglx,2)
write(31) xs,ys
write(31) mxvrx
write(31) bsvar
write(31) mxdis
write(31) anglx
write(31) aniso
write(31) mnslp
write(31) angll
write(31) xspk
write(31) yspk
write(31) mxds0
write(31) mxds1
write(31) sft0
write(31) sft1
write(31) hwdth
write(31) npks
write(31) mxvry

write(31) nvls
write(31) pkhts
write(31) vldps
write(31) rwpks
write(31) rwvls

write(31) lon0
write(31) lon1
write(31) lat0
write(31) lat1

write(31) uqrid
write(31) riseq
write(31) fallq
write(31) clngth

write(31) mxds2

   CLOSE(31)
!!!STOP
#endif

 

end subroutine find_ridges
!----------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ANISO_ANA( AA,AARAW,X,Y,N,NSB,NSW,IP)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER,            intent(IN)  ::  N,NSB,NSW,IP
  REAL,               intent(IN)  ::  AA(N,N),X(N),Y(N),AARAW(N,N)


  real, allocatable :: SUBA(:,:), SUBARW(:,:),SUBX(:),SUBY(:),RT(:,:),RTX(:),XRT(:),RTXMN(:),RTXSLP(:),rtx_dt(:)
  real, allocatable :: PKLC(:), RTY(:),RTRW(:,:),RTRWX(:),DERTX(:),DERTY(:),CUSP(:),face(:)
  real, allocatable :: silux(:), sillx(:), siluy(:), silly(:)

  logical,allocatable :: lhgts(:),lflats(:),lsides(:)
  logical :: Keep_Cuestas=.true.
  

  real :: THETRAD,PI,swt,ang,rotmn,rotvar,mnt,var,xmn,xvr,basmn,basvar,mn2,var2
  integer :: i,j,l,m,n2,mini,maxi,minj,maxj,ns0,ns1,iorn(1),jj, &
            ipkh(1),ivld(1),ift0(1),ift1(1),i2,ii,ipksv(1)

  real :: vvaa(NANG),qual(NANG),dex(NANG),beta(NANG),alph,xpkh(NANG),ang00
  real :: dex0(nang),dex1(nang),xft0(NANG),xft1(NANG),HWDX(NANG),xvld(NANG)
        
  real :: NPKX(NANG),NVLX(NANG),vva2(NANG),pkht(NANG),vldp(NANG),rwpk(NANG)
  real :: rwvl(NANG),RISEX(NANG),FALLX(NANG),LNGTH(NANG) 
  real :: dex_dt(NANG)

  ns0=nsw/2+1
  ns1=ns0+nsw+1
 

  PI = 2*ACOS(0.0)

! Allocate work arrays for ridge analysis
!-----------------------------------------
  allocate( suba( 2*nsw+1, 2*nsw+1 ) )
  allocate( subarw( 2*nsw+1, 2*nsw+1 ) )
  allocate( rt( 2*nsw+1, 2*nsw+1 ) )
  allocate( rtrw( 2*nsw+1, 2*nsw+1 ) )
  allocate( subx( 2*nsw+1  ) )
  allocate( suby( 2*nsw+1  ) )
  allocate( rtx( nsw+1  ) )
  allocate( rtx_dt( nsw+1  ) )
  allocate( rty( nsw+1  ) )
  allocate( cusp( nsw+1  ) )
  allocate( dertx( nsw+1-1  ) )
  allocate( derty( nsw+1-1  ) )
  allocate( rtrwx( nsw+1  ) )
  allocate( rtxslp( nsw+1-1  ) )
  allocate( pklc( 2:nsw+1-1  ) )
  allocate( rtxmn( nsw+1  ) )
  allocate( xrt( nsw+1  ) )
  allocate( face( nsw+1  ) )

  allocate( silux( nsw+1  ) )
  allocate( sillx( nsw+1  ) )
  allocate( siluy( nsw+1  ) )
  allocate( silly( nsw+1  ) )

  allocate( Lhgts( nsw+1  ) )
  allocate( Lflats( nsw+1  ) )
  allocate( Lsides( nsw+1  ) )


  n2=n/2
  mini=nsw/nsb+1
  maxi=N/nsb-nsw/nsb   !+1
  minj=nsw/nsb+1
  maxj=N/nsb-nsw/nsb   !+1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check ranges for orientation analysis
  !***************************************************
     write(*,*) "CHECKING RANGES in topoana ", nsb,nsw
     write(*,*) "  N,nsb,nsw ", N,nsb,nsw


    if (mini*nsb-nsw+1<1) then 
       write(*,*) "fixing bad mini ", mini*nsb-nsw+1, mini
       mini=mini+1
       if (mini*nsb-nsw+1<1) then 
          write(*,*) "irreparable mini ", mini*nsb-nsw+1, mini
          STOP
       endif
    endif
    if (maxi*nsb+nsw>N) then 
       write(*,*) "fixing bad maxi ", maxi*nsb+nsw , maxi
       maxi=maxi-1
       if (maxi*nsb+nsw>N) then 
          write(*,*) "irrep. maxi ", maxi*nsb+nsw , maxi
          STOP 
      endif
    endif
    if (minj*nsb-nsw+1<1) then 
       write(*,*) "fixing bad minj ", minj*nsb-nsw+1
       minj=minj+1
       if (minj*nsb-nsw+1<1) then 
          write(*,*) "irrep. minj ", minj*nsb-nsw+1
          STOP
       endif
    endif
    if (maxj*nsb+nsw>N) then 
       write(*,*) "fixing bad maxj ", maxj*nsb+nsw , maxi
       maxj=maxj-1
       if (maxj*nsb+nsw>N) then 
          write(*,*) "irrep. maxj ", maxj*nsb+nsw , maxi
          STOP
       endif
    endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do i=1,nsw+1 
     xrt(i) = 1.*(i-1)
  end do
  xmn=0.5*(xrt(1)+xrt(nsw+1))
  xvr=sum( (xrt-xmn)**2 )/nsw+1

  do j=minj,maxj
     suby = Y(j*nsb-nsw  :j*nsb+nsw )
     ys(j)= sum( suby )/size(suby,1)
     do i=mini,maxi
        subx = x(i*nsb-nsw  :i*nsb+nsw )
        xs(i)= sum( subx )/size(subx,1)
#ifdef DEBUGRIDGE
        if ( (ys(j)>=2600).and.(ys(j)<=3000).and. &
             (xs(i)>= 600).and.(xs(i)<=1300) ) then
#endif
        suba = AA( i*nsb-nsw  :i*nsb+nsw ,  j*nsb-nsw  :j*nsb+nsw )

#ifdef SPEEDUPLOOP
        if( (maxval(suba)-minval(suba))>1.0 ) then
#endif
        subarw = AARAW( i*nsb-nsw  :i*nsb+nsw ,  j*nsb-nsw  :j*nsb+nsw )

        basmn  =  sum( sum( suba(ns0:ns1-1,ns0:ns1-1) , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
        basvar =  sum( sum( (suba(ns0:ns1-1,ns0:ns1-1)-basmn)**2 , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
        bsvar(i,j,ip) = basvar

        vvaa(:)=0
        qual(:)=0.
        npkx(:)=0.
        do l=1,nang
                                           ! Rotate 2D topography by ANG
           ang = (L-1)*(180./nang)
           rt  = rotby3( suba, 2*nsw+1 , ang )
           rtrw  = rotby3( subarw, 2*nsw+1 , ang )  ! Raw topo rotation

                                           ! Take "Y" (and "X")-average of rotated topography.
                                           ! Yields topo profile in X ==> RTX
           rtx = sum( rt(ns0:ns1-1,ns0:ns1-1) , 2 ) /( ns1-ns0 ) ! Y-average 
           rty = sum( rt(ns0:ns1-1,ns0:ns1-1) , 1 ) /( ns1-ns0 ) ! X-average
           rtrwx = sum( rtrw(ns0:ns1-1,ns0:ns1-1) , 2 ) /( ns1-ns0 ) ! Y-average of Raw topo

#if 1 
                                          ! "Silhouettes" in x and y
           do m=ns0,ns1-1
              silux(m-ns0+1) = maxval(  rt( m ,ns0:ns1-1) )
              sillx(m-ns0+1) = minval(  rt( m ,ns0:ns1-1) )
              siluy(m-ns0+1) = maxval(  rt( ns0:ns1-1, m) )
              silly(m-ns0+1) = minval(  rt( ns0:ns1-1, m) )
           end do
#endif

                                           
                                           ! Mean elevation
           mnt = sum( rtx )/( ns1-ns0 ) 
           mn2 = sum( rty )/( ns1-ns0 ) 



                  ! Mean slope BETA (and intercept ALPH) of RTX
           beta(L) = sum( (rtx-mnt)*(xrt-xmn) )/(nsw*xvr)           
           alph    = mnt - beta(L)*xmn
               ! subtract linear slope
               do ii=1,nsw+1
                  rtx_dt(ii) = rtx(ii) - ( alph+beta(L)*xrt(ii) )
               end do


           rtxslp(1:nsw)=abs(rtx(2:nsw+1)-rtx(1:nsw) )                         
  

                 ! count actual peaks and valleys in RTX cross section
           pklc = 0.
           do i2=2,nsw
              if ( ( rtx(i2-1)<rtx(i2) ).and.( rtx(i2+1)<rtx(i2) ) ) pklc(i2)=1.
           end do
           npkx(L)=sum(pklc)

           pklc = 0.
           do i2=2,nsw
             if ( ( rtx(i2-1)>rtx(i2) ).and.( rtx(i2+1)>rtx(i2) ) ) pklc(i2)=1.
           end do
           nvlx(L)=sum(pklc)


                ! Accumulate rising and falling segments
           risex(L)=0.
           fallx(L)=0.
           do i2=2,nsw+1
              if ( rtx(i2-1)<=rtx(i2) ) risex(L) = risex(L) + ( rtx(i2) - rtx(i2-1) )
              if ( rtx(i2-1)>=rtx(i2) ) fallx(L) = fallx(L) + ( rtx(i2) - rtx(i2-1) )
           end do
           do i2=2,nsw+1
              dertx(i2-1) = rtx(i2) - rtx(i2-1)
              derty(i2-1) = rty(i2) - rty(i2-1)
           end do


                ! Record actual max and min elevations in RTX and Raw topo profile (RTRWX)
           pkht(L)=maxval(RTX)
           vldp(L)=minval(RTX)
           rwpk(L)=maxval(RTRWX)
           rwvl(L)=minval(RTRWX)






           var = sum( (rtx-mnt)**2 )/( ns1-ns0 ) 
           vvaa(L) = var
           dex(L)  = MAXVAL(RTX)-MINVAL(RTX)
           dex_dt(L)  = MAXVAL(RTX_dt)-MINVAL(RTX_dt)


           var2 = sum( (rty-mn2)**2 )/( ns1-ns0 ) 
           vva2(L) = var2

           Lhgts = ( ( rtx - minval(rtx) ) > 0.5*maxval( rtx - minval(rtx) ) )
           !hwdx(L) = 1.0*count( Lhgts )


              ! flats = where slope is less than 33% of max
           Lflats = ( abs(dertx) < 0.33*maxval( abs(dertx ) ) ) 
           hwdx(L) = 1.0*(nsw+1 - count( Lflats ) )

              ! sides = where sides are steeper than 20% of max
           Lsides = ( abs(derty) > 0.2*maxval( abs(dertx ) ) ) 
           lngth(L) = 1.0*(nsw+1 - count( Lsides ) )



#if 0
           if (dex(L) > 1.0) then 
              !hwdx(L) = sum( rtxslp**2  , 1) / (dex(L)**2)            
              hwdx(L) = ( dex(L)**2 ) / sum( rtxslp**2  , 1)            
           else
              hwdx(L) = -1.
           endif
#endif

           ipkh    = MAXLOC( RTX ) ! index of peak height in rotated topo avg cross-section
           xpkh(L) = XRT( ipkh(1) )-xmn
           ivld    = MINLOC( RTX ) ! index of valley depth in rotated topo avg cross-section
           xvld(L) = XRT( ivld(1) )-xmn


           ! if ipkh is at or close to the edge of the profile
           ! i.e., 1 or nsw, then is likely sloping terrain not
           ! a real peak-y feature. Flag it by zeroing out npkx. 
           ! Risks missing a true peak in the interior, but redundant
           ! subsamples should help. Bigger risk is getting fooled
           ! small (spurious) interior peak. 
           if( ( ipkh(1) <= max( nsw/2-nsb,1) ).or. &
               ( ipkh(1) >= min( nsw/2+nsb,nsw) ) ) then
               npkx(L)=0.0
           end if
           
 
           ipksv = ipkh

           ift0    = MINLOC( RTX( 1 : ipkh(1) ) )
           ift1    = MINLOC( RTX( ipkh(1) :  ) )+ipkh(1)-1

           dex0(L)  = MAXVAL(RTX)-MINVAL( RTX( 1 : ipkh(1) ) )
           dex1(L)  = MAXVAL(RTX)-MINVAL( RTX( ipkh(1) :  ) )

           xft0(L)  =  XRT( ift0(1) ) -XRT( ipkh(1) )
           xft1(L)  =  XRT( ift1(1) ) -XRT( ipkh(1) )

 !!        basmn  =  sum( sum( suba(ns0:ns1-1,ns0:ns1-1) , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
 !!        basvar =  sum( sum( (suba(ns0:ns1-1,ns0:ns1-1)-basmn)**2 , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
           rotmn  =  sum( sum( rt(ns0:ns1-1,ns0:ns1-1) , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
           rotvar =  sum( sum( (rt(ns0:ns1-1,ns0:ns1-1)-rotmn)**2 , 1 ), 1) /(( ns1-ns0 )*(ns1-ns0))
           if (rotvar>0.) qual(L) = var/rotvar
           if (rotvar>0.) vvaa(L) = vvaa(L)*(basvar/rotvar)
           if (rotvar>0.) vva2(L) = vva2(L)*(basvar/rotvar)


        !++++ Analysis using triangular profiles in X and Y
           do i2=1,ipkh(1)
              cusp(i2) = rtx(ipkh(1)) + (rtx(1)-rtx(ipkh(1)))*(i2-ipkh(1))/(1-ipkh(1)) 
           end do
           do i2=ipkh(1),nsw+1
              cusp(i2) = rtx(ipkh(1)) + (rtx(nsw)-rtx(ipkh(1)))*(i2-ipkh(1))/(nsw-ipkh(1)) 
           end do
 
           ipkh    = MAXLOC( RTY ) ! index of ridge "FACE"
           do i2=1,ipkh(1)
              face(i2) = rty(ipkh(1)) + (rty(1)-rty(ipkh(1)))*(i2-ipkh(1))/(1-ipkh(1)) 
           end do
           do i2=ipkh(1),nsw 
              face(i2) = rty(ipkh(1)) + (rty(nsw)-rty(ipkh(1)))*(i2-ipkh(1))/(nsw-ipkh(1)) 
           end do



           !if ( (ip == 2) .and. ( i >= 39 ).and.( i <= 39) .and. ( j>=67 ) .and. (j<=67) ) then
           !if ( (ip == 4) .and. ( i == 33 ).and.( j== 75 ) ) then
           if ( (ip == 4) .and. (( i == 129 )).and.(( j== 341 ).or.(j==340)) ) then
               write(32) i,j,ip,nsw
               write(32) ang, rt,rtx,subx,suby,rtrw,rtrwx,rty,dertx,derty,cusp,face,suba,subarw & 
                        ,silux,sillx,siluy,silly,rtx_dt
               write(32) vvaa(L),dex(L),risex(L),fallx(L),hwdx(L),vva2(L),qual(L),basvar,rotvar & 
                        ,npkx(L),nvlx(L),pkht(L),vldp(L),lngth(L),xpkh(L),xvld(L),xmn,ipksv,ivld
           endif 




        end do ! LOOP over angles



        iorn       = MAXLOC( vvaa )

        mxvrx(i,j,ip) = vvaa( iorn(1) ) !MAXVAL( vvaa )
        mxvry(i,j,ip) = vva2( iorn(1) ) !MAXVAL( vvaa )
        mxdis(i,j,ip) = dex( iorn(1) )
        mxds2(i,j,ip) = dex_dt( iorn(1) )
        hwdth(i,j,ip) = hwdx( iorn(1) )
        npks(i,j,ip)  = npkx( iorn(1) )
        nvls(i,j,ip)  = nvlx( iorn(1) )
        aniso(i,j,ip) = qual( iorn(1) )
        anglx(i,j,ip) = (iorn(1)-1)*(180./nang)
        mnslp(i,j,ip) = beta( iorn(1) )

        pkhts(i,j,ip) = pkht( iorn(1) )
        vldps(i,j,ip) = vldp( iorn(1) )
 
        rwpks(i,j,ip) = rwpk( iorn(1) )
        rwvls(i,j,ip) = rwvl( iorn(1) )

        ang00      =  (iorn(1)-1)*(180./nang)*(PI/180.)

        xspk(i,j,ip)  = xs(i)  + xpkh(iorn(1))*cos( ang00 )
        yspk(i,j,ip)  = ys(j)  - xpkh(iorn(1))*sin( ang00 )

        mxds0(i,j,ip) = dex0( iorn(1) )
        mxds1(i,j,ip) = dex1( iorn(1) )
        sft0(i,j,ip)  = xft0( iorn(1) )
        sft1(i,j,ip)  = xft1( iorn(1) )

        riseq(i,j,ip) = risex( iorn(1) )
        fallq(i,j,ip) = fallx( iorn(1) )
        clngth(i,j,ip)= lngth( iorn(1) )


#ifdef SPEEDUPLOOP
      END IF 
#endif
#ifdef DEBUGRIDGE
      END IF 
#endif
     end do
         write(*,*) " J", j," out of ", N/nsb
  end do

if ( .not.(KEEP_CUESTAS)) then
! Clean up "goofy" results for inclined planes
  where ( (nvls == 0) .and. (npks == 0) )
        mxdis = 0.
        aniso = 0.
        mxvrx = 0.
        mxvry = 0.
  endwhere
else
  write(*,*) " Keeping goofy 'Cuestas' "
endif



#if 0
  deallocate( suba )
  deallocate( rt )
  deallocate( subx )
  deallocate( suby )
  deallocate( rtx )
  deallocate( rty )
  deallocate( dertx )
  deallocate( derty )
#endif



  deallocate( suba)
  deallocate( subarw)
  deallocate( rt)
  deallocate( rtrw)
  deallocate( subx)
  deallocate( suby)
  deallocate( rtx)
  deallocate( rty)
  deallocate( cusp)
  deallocate( dertx)
  deallocate( derty)
  deallocate( rtrwx)
  deallocate( lhgts)
  deallocate( rtxslp)
  deallocate( pklc)
  deallocate( rtxmn)
  deallocate( xrt)
  deallocate( face)
  deallocate( lflats)
  deallocate( lsides)
  deallocate( silux)
  deallocate( siluy)
  deallocate( sillx)
  deallocate( silly)


write(*,*) "Got though ANISO ANA"

end subroutine ANISO_ANA

!====================================
   subroutine remapridge2target(area_target,target_center_lon,target_center_lat,  &
         weights_eul_index_all,weights_lgr_index_all,weights_all,ncube,jall,&
         nreconstruction,ntarget,nhalo,nsb,nsw,nsmcoarse,nsmfine,lzerovalley)
      use shr_kind_mod, only: r8 => shr_kind_r8
      use remap
      use reconstruct !, only : EquiangularAllAreas
      implicit none
      real(r8), intent(in) :: weights_all(jall,nreconstruction)
      integer , intent(in) :: weights_eul_index_all(jall,3),weights_lgr_index_all(jall)
      integer , intent(in) :: ncube,jall,nreconstruction,ntarget,nhalo,nsb,nsw,nsmcoarse,nsmfine
      real(r8), intent(in) :: area_target(ntarget),target_center_lon(ntarget),target_center_lat(ntarget)
      logical, intent(in)  :: lzerovalley
      real(r8):: f(ntarget)
  
      REAL  ,                                                          &
         DIMENSION(1-nhalo:ncube+nhalo )                          :: xv,yv,alph,beta
      
      integer :: alloc_error

      integer :: i,ix,iy,ip,ii,counti,norx,nory,i_last,isubr,iip,j
      real(r8):: wt
      real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: tmpx6
      real(KIND=dbl_kind), dimension(ncube*ncube*6) :: mxdisC , anglxC, anisoC, hwdthC
      real(KIND=dbl_kind), dimension(ncube*ncube*6) :: mxvrxC , mxvryC, bsvarC, clngtC
      real(KIND=dbl_kind), dimension(ncube*ncube*6) :: cwghtC , itrgtC, fallqC, riseqC, rwpksC
      real(KIND=dbl_kind), dimension(ncube*ncube)   :: dA

      CHARACTER(len=1024) :: ofile$
!----------------------------------------------------------------------------------------------------

    !!allocate ( dA(ncube,ncube),stat=alloc_error )
    CALL EquiangularAllAreas(ncube, dA)

    DO i=1-nhalo,ncube+nhalo
       xv(i)=1.*i
       yv(i)=1.*i
    END DO
 
    allocate (wghts_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for wghts_target'; stop; endif
    wghts_target = 0.
    allocate (mxdis_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxdis_target'; stop; endif
    mxdis_target = 0.
    allocate (anglx_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for anglx_target'; stop; endif
    anglx_target = 0.
    allocate (aniso_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for aniso_target'; stop; endif
    aniso_target = 0.
    allocate (anixy_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for anixy_target'; stop; endif
    anixy_target = 0.
    allocate (mxvrx_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxvrx_target'; stop; endif
    mxvrx_target = 0.
    allocate (mxvry_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxvry_target'; stop; endif
    mxvry_target = 0.
    allocate (bsvar_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for bsvar_target'; stop; endif
    bsvar_target = 0.
    allocate (hwdth_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for hwdth_target'; stop; endif
    hwdth_target = 0.
    allocate (clngt_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for clngt_target'; stop; endif
    clngt_target = 0.
    allocate (cwght_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for cwght_target'; stop; endif
    cwght_target = 0.
    allocate (count_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for count_target'; stop; endif
    count_target = 0.
    allocate (riseq_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for riseq_target'; stop; endif
    riseq_target = 0.
    allocate (fallq_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for fallq_target'; stop; endif
    fallq_target = 0.

     norx=size(mxdis,1)
     nory=size(mxdis,2)

     itrgtC = 0.



    if(lzerovalley)then
      write(*,*) " will ZERO out negative peaks and 'Cuestas' "
      do ip=1,6
      do j=1,nory
      do i=1,norx
         if( (pkhts(i,j,ip)<0.1*mxdis(i,j,ip)).or.(npks(i,j,ip)<1.0) ) then
           mxdis(i,j,ip)  = 0.
           riseq(i,j,ip)  = 0.
           fallq(i,j,ip)  = 0.
           clngth(i,j,ip) = 0.
         endif
      end do
      end do
      end do
    else
      write(*,*) " Leave negative peaks ALONE "
    endif

 



#if 0
     tmpx6 = mapridge2cube ( mxdis , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxdisC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( anglx , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     anglxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( aniso , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     anisoC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( mxvrx , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxvrxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( mxvry , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxvryC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( bsvar , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     bsvarC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( hwdth , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     hwdthC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     clngtC = 0. 
     cwghtC = 0. 
#else
     !where(pkhts < 0.)
     !   mxdis = 0.
     !endwhere

        write(*,*) " about to call paintridge2cube "
     tmpx6 = paintridge2cube ( mxdis ,  ncube,nhalo,nsb,nsw,lzerovalley )
     mxdisC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( anglx ,  ncube,nhalo,nsb,nsw,lzerovalley )
     anglxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( aniso ,  ncube,nhalo,nsb,nsw,lzerovalley )
     anisoC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( mxvrx ,  ncube,nhalo,nsb,nsw,lzerovalley )
     mxvrxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( mxvry ,  ncube,nhalo,nsb,nsw,lzerovalley )
     mxvryC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( bsvar ,  ncube,nhalo,nsb,nsw,lzerovalley )
     bsvarC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( hwdth ,  ncube,nhalo,nsb,nsw,lzerovalley )
     hwdthC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( clngth ,  ncube,nhalo,nsb,nsw,lzerovalley, crest_length=.true. )
     clngtC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( clngth ,  ncube,nhalo,nsb,nsw,lzerovalley, crest_weight=.true. )
     cwghtC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( fallq ,  ncube,nhalo,nsb,nsw,lzerovalley )
     fallqC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = paintridge2cube ( riseq ,  ncube,nhalo,nsb,nsw,lzerovalley )
     riseqC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
#endif

       !mxdis_target = remap_field(tmp,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       !weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)





    i_last = -9999
!++jtb 
! 
!      In the following loop "count" is the index of a piece of 
!      the "exchange grid" - created by cutting the cubed-sphere topo
!      and target grid into each other.  
    do counti=1,jall
     
      i    = weights_lgr_index_all(counti)
      !if( itracker_target(i) == 0 ) then ! If true this is first time target cell i has been visited
      !   itracker_target(i)=1
      !end if

      ix  = weights_eul_index_all(counti,1)
      iy  = weights_eul_index_all(counti,2)
      ip  = weights_eul_index_all(counti,3)
      !
      ! convert to 1D indexing of cubed-sphere
      !
      ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
      wt = weights_all(counti,1)

      iip=(iy-1)*ncube+ix

      itrgtC(ii) = i

      isubr = INT( anglxC(ii) * nsubr/180. ) + 1

      if ( (isubr >= 1).and.(isubr <= nsubr) ) then
      wghts_target( i , isubr ) = wghts_target( i , isubr ) + wt
      hwdth_target( i , isubr ) = hwdth_target( i , isubr ) + wt*hwdthC(ii) *cwghtC(ii)
      mxvrx_target( i , isubr ) = mxvrx_target( i , isubr ) + wt*mxvrxC(ii) *cwghtC(ii)
      mxvry_target( i , isubr ) = mxvry_target( i , isubr ) + wt*mxvryC(ii) *cwghtC(ii)
      bsvar_target( i , isubr ) = bsvar_target( i , isubr ) + wt*bsvarC(ii) *cwghtC(ii)
      mxdis_target( i , isubr ) = mxdis_target( i , isubr ) + wt*mxdisC(ii) *cwghtC(ii)
      fallq_target( i , isubr ) = fallq_target( i , isubr ) + wt*fallqC(ii) *cwghtC(ii)
      riseq_target( i , isubr ) = riseq_target( i , isubr ) + wt*riseqC(ii) *cwghtC(ii)
      aniso_target( i , isubr ) = aniso_target( i , isubr ) + wt*anisoC(ii) *cwghtC(ii)
      anglx_target( i , isubr ) = anglx_target( i , isubr ) + wt*anglxC(ii) *cwghtC(ii)
      !! clngt_target( i , isubr ) = clngt_target( i , isubr ) + wt*clngtC(ii) *cwghtC(ii)
      clngt_target( i , isubr ) = clngt_target( i , isubr ) + wt*cwghtC(ii)/dA(iip)
      cwght_target( i , isubr ) = cwght_target( i , isubr ) + wt*cwghtC(ii) 
      count_target( i , isubr ) = count_target( i , isubr ) + wt/dA(iip)
      endif

      i_last = i
    end do       

    ! change width (and length) to km
    hwdth_target = hwdth_target * grid_length_scale
    clngt_target = clngt_target * grid_length_scale

    ! Make fallq positive
    fallq_target = -1.*fallq_target

#if 0
    where( wghts_target > 1.e-15 )
        aniso_target = aniso_target / wghts_target
        anglx_target = anglx_target / wghts_target
        mxdis_target = mxdis_target / wghts_target
        hwdth_target = hwdth_target / wghts_target
        mxvrx_target = mxvrx_target / wghts_target
        mxvry_target = mxvry_target / wghts_target
        bsvar_target = bsvar_target / wghts_target
     elsewhere
        aniso_target = 0.
        anglx_target = -9000.
        mxdis_target = 0.
        hwdth_target = 0.
        mxvrx_target = 0.
        mxvry_target = 0.
        bsvar_target = 0.
     end where
#else
     where( cwght_target > 1.e-15 )
        !clngt_target = clngt_target / cwght_target
        mxdis_target = mxdis_target / cwght_target
        fallq_target = fallq_target / cwght_target
        riseq_target = riseq_target / cwght_target
        aniso_target = aniso_target / cwght_target
        anglx_target = anglx_target / cwght_target
        hwdth_target = hwdth_target / cwght_target
        mxvrx_target = mxvrx_target / cwght_target
        mxvry_target = mxvry_target / cwght_target
        bsvar_target = bsvar_target / cwght_target
     elsewhere      
        clngt_target = 0.
        mxdis_target = 0.
        riseq_target = 0.
        fallq_target = 0.
        aniso_target = 0.
        anglx_target = -9000.
        hwdth_target = 0.
        mxvrx_target = 0.
        mxvry_target = 0.
        bsvar_target = 0.
     end where
#endif


#if 1
    anixy_target = mxvrx_target /( mxvrx_target + mxvry_target + 0.0001 )
    !!anixy_target = 2.*( ( mxvrx_target /( mxvrx_target + mxvry_target + 0.0001) )-0.5 )
#endif
     

     call importancesort (ntarget)

     call latlonangles (target_center_lon,target_center_lat,ntarget)


       write( ofile$ , &
       "('./output/remap_nc',i0.4, '_Nsw',i0.3,'_Nrs',i0.3  &
       '_Co',i0.3,'_Fi',i0.3)" ) & 
        ncube, nsw, nsb, nsmcoarse, nsmfine
       ofile$= trim(ofile$)//'.dat'

       OPEN (unit = 911, file= trim(ofile$) ,form="UNFORMATTED" )

write(911) ncube
write(911) mxdisC
write(911) mxvrxC
write(911) mxvryC
write(911) anglxC
write(911) hwdthC
write(911) cwghtC
write(911) clngtC
write(911) itrgtC
write(911) fallqC
write(911) riseqC

close(911)


      write(*,*) " GOT OUT OF remapridge2target "

 
   end subroutine remapridge2target


!================================================================
  subroutine latlonangles (target_center_lon,target_center_lat,ntarget)
!-----------------------------------

      real(r8), intent(in) :: target_center_lon(ntarget),target_center_lat(ntarget)
      integer,  intent(in) :: ntarget

  real(r8) :: lat22,lon22,a22,b22,dx2,dy2,cosll
  real(r8) :: lat22s,lon22s,a22s,b22s

  integer  :: i,j,isubr,ipanel22
  integer  :: alloc_error


    allocate (ang22_target(ntarget,nsubr),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for ang22_target'; stop; endif
    ang22_target = 0.


! Calculate angles in lat-lon system

  do i=1,ntarget
     lon22 = target_center_lon(i) * PI/180.
     lat22 = target_center_lat(i) * PI/180.
     call CubedSphereABPFromRLL(lon22, lat22, a22, b22, ipanel22 , .true. )

     do isubr=1,nsubr
        if ( ANGLX_TARGET(i,isubr) > -9000. ) then
          a22s = a22 + 0.01*SIN( ANGLX_TARGET(i,isubr)*PI/180. )   
          b22s = b22 + 0.01*COS( ANGLX_TARGET(i,isubr)*PI/180. )
          call CubedSphereRLLFromABP(a22s, b22s , ipanel22, lon22s, lat22s )
          dx2 = COS( lat22 )*(lon22s-lon22 )
          dy2 = ( lat22s-lat22 )
          if ( dx2 < 0.0 ) dy2  = -1.*dy2  ! 
          COSLL  = dy2 /sqrt( dx2**2 + dy2**2 )
          ANG22_TARGET(i,isubr) = ACOS( COSLL )*180./PI
        else
          ANG22_TARGET(i,isubr) = -9000.
        end if   
     end do
  end do

 end subroutine latlonangles






!================================================================
  subroutine importancesort (ntarget)
!-----------------------------------

     integer,  intent(in) :: ntarget


  real(r8), allocatable, dimension(:) :: imprtnc,tmp2
  integer,  allocatable, dimension(:) :: insrt

  
  real(r8) :: tmp
  integer  :: i,ii,jj,kk,itmp
  integer  :: alloc_error

 
    allocate( insrt(nsubr) , stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for INSRT '; stop; endif
    allocate( imprtnc(nsubr) , stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for IMPRTNC '; stop; endif
    allocate( tmp2(nsubr) , stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for TMP2 '; stop; endif



  do i=1,ntarget

     do JJ = 1,nsubr
        insrt(JJ)   = JJ
        imprtnc(JJ) = mxvrx_target(i,JJ) * wghts_target(i,JJ) * aniso_target(i,JJ)
     end do

     do JJ = 1,nsubr-1
        do KK = JJ+1,nsubr
           if ( imprtnc(JJ) < imprtnc(KK) ) then
              tmp = imprtnc(KK)
              imprtnc(KK) = imprtnc(JJ)
              imprtnc(JJ) = tmp
              itmp = insrt(KK)
              insrt(KK)   = insrt(JJ)
              insrt(JJ)   = itmp
           end if
        end do
     end do


     do JJ = 1,nsubr
        tmp2(JJ) =  mxvrx_target(i,JJ)
     end do
     do JJ = 1,nsubr
        mxvrx_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  mxvry_target(i,JJ)
     end do
     do JJ = 1,nsubr
        mxvry_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  mxdis_target(i,JJ)
     end do
     do JJ = 1,nsubr
        mxdis_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  riseq_target(i,JJ)
     end do
     do JJ = 1,nsubr
        riseq_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  fallq_target(i,JJ)
     end do
     do JJ = 1,nsubr
        fallq_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  bsvar_target(i,JJ)
     end do
     do JJ = 1,nsubr
        bsvar_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  aniso_target(i,JJ)
     end do
     do JJ = 1,nsubr
        aniso_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  anixy_target(i,JJ)
     end do
     do JJ = 1,nsubr
        anixy_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  anglx_target(i,JJ)
     end do
     do JJ = 1,nsubr
        anglx_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  hwdth_target(i,JJ)
     end do
     do JJ = 1,nsubr
        hwdth_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  wghts_target(i,JJ)
     end do
     do JJ = 1,nsubr
        wghts_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  cwght_target(i,JJ)
     end do
     do JJ = 1,nsubr
        cwght_target(i,JJ) = tmp2( insrt(JJ) )
     end do

     do JJ = 1,nsubr
        tmp2(JJ) =  clngt_target(i,JJ)
     end do
     do JJ = 1,nsubr
        clngt_target(i,JJ) = tmp2( insrt(JJ) )
     end do

  end do


end subroutine importancesort





!======================================

function mapridge2cube ( a, norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb ) result( axc )
   
       integer, intent(in) :: norx,nory,ncube,nhalo,nsb
       real, intent(in), dimension(norx,nory,6) :: a
       real, intent(in), dimension(norx) :: xs
       real, intent(in), dimension(nory) :: ys
       real, intent(in), DIMENSION(1-nhalo:ncube+nhalo ) :: xv,yv

       real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: axc
       integer:: i,j,x0,x1,y0,y1,ip
!---------------------------------------------------

     do ip=1,6
       do j=1,nory
          y0 = int( ys(j) )+1 - NSB/2. 
          y1 = int( ys(j) )  +  NSB/2. 
          do i=1,norx
             x0 = int( xs(i) )+1 - NSB/2. 
             x1 = int( xs(i) )  +  NSB/2.
             if( (x0>-900).and.(x1>-900).and.(y0>-900).and.(y1>-900) ) then
                AXC( x0:x1, y0:y1, ip ) = A( i , j , ip )
             endif
          end do
       end do
     end do

end function mapridge2cube

!======================================

function paintridge2cube ( axr, ncube,nhalo,nsb,nsw, lzerovalley, crest_length, crest_weight ) result( axc )
   
       integer, intent(in) :: ncube,nhalo,nsb,nsw
       real, intent(in), dimension( size(xs) , size(ys) ,6) :: axr
       logical, intent(in) :: lzerovalley
       logical, optional, intent(in) :: crest_length
       logical, optional, intent(in) :: crest_weight

       real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: axc
       real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: qc
       real, dimension(-nsw:nsw,-nsw:nsw) :: suba,sub1
       real, dimension(-nsw:nsw,-nsw:nsw) :: subr,subq,subdis
       real, dimension(-nsw:nsw)          :: xq,yq
       real :: rotangl,dsq,ssq
       integer :: i,j,x0,x1,y0,y1,ip,ns0,ns1,ii,jj,norx,nory,nql,ncl,nhw
       logical :: lcrestln,lcrestwt
!---------------------------------------------------

    if(present(crest_length)) then
      lcrestln = crest_length
    else
       lcrestln = .false.
    endif
    if(present(crest_weight)) then
      lcrestwt = crest_weight
    else
       lcrestwt = .false.
    endif

    DO i=-nsw,nsw
       xq(i)=i
       yq(i)=i
    END DO
 

  ns0=nsw/2+1
  ns1=ns0+nsw

     norx=size(mxdis,1)
     nory=size(mxdis,2)



  if ( .not.(Lcrestln) ) then
   suba(:,:)=0.
   suba( -nsw/4-1:nsw/4+1 , -nsw/2-1:nsw/2+1 ) = 1.
   !else 
   ! suba( 0 , -nsw/2-1:nsw/2+1 ) = 1.
  endif  

  sub1(:,:)=0.
  nql = INT( nsw/1.0 )
    DO j=-nsw,nsw
    DO i=-nsw,nsw
       ssq= xq(i)**2+yq(j)**2
       if (SSQ < nql) sub1(i,j)=1.0
    END DO
    END DO
 
  axc = 0.
  qc  = 0.

#ifdef DEBUGRIDGE
     do ip=4,4
#else
     do ip=1,6
#endif
       do j=1,nory
         do i=1,norx
#ifdef DEBUGRIDGE
        if ( (ys(j)>=2600).and.(ys(j)<=3000).and. &
             (xs(i)>= 600).and.(xs(i)<=1300) ) then
#endif
        if(mxdis(i,j,ip)>=1.0) then

             if(Lcrestwt) then
               suba(:,:) = 0.
               ncl  = MIN( INT(clngth(i,j,ip)/2) , nsw/2 )
               suba( 0 , -ncl:ncl ) = 1.        
               rotangl = - anglx(i,j,ip) 
               subr = rotby3( suba , 2*nsw+1, rotangl )
               subdis = subr  ! * axr(i,j,ip)
             end if
             if(Lcrestln) then
               suba(:,:) = 0.
               ncl  = MIN( INT(clngth(i,j,ip)/2) , nsw/2 )
               suba( 0 , -ncl:ncl ) = 1.        
               rotangl = - anglx(i,j,ip) 
               subr = rotby3( suba , 2*nsw+1, rotangl )
               subdis = subr * axr(i,j,ip)
             end if
             if( (.not.(Lcrestln)).and.(.not.(Lcrestwt)) ) then
               suba(:,:) = 0.
               ncl  = MIN( INT(clngth(i,j,ip)/2) , nsw/2 )
               nhw  = MIN( INT(hwdth(i,j,ip)/2) , nsw/2 )
               !suba( -nhw:nhw , -ncl:ncl ) = 1.        
               suba( 0 , -ncl:ncl ) = 1.        
               rotangl = - anglx(i,j,ip) 
               subr = rotby3( suba , 2*nsw+1, rotangl )
               subdis = subr * axr(i,j,ip)
             endif



             dsq    = 1.0 - SQRT( (xs(i)-xspk(i,j,ip))**2 + (ys(j)-yspk(i,j,ip))**2 )/nsw
             subq   = sub1 * dsq

                do jj = -NSW/2,NSW/2
                do ii = -NSW/2,NSW/2
                    x0 = INT( xspk(i,j,ip) ) + 1
                    y0 = INT( yspk(i,j,ip) ) + 1
                    if ( (x0+ii>=1-nhalo).and.(x0+ii<=ncube+nhalo).AND.(Y0+ii>=1-nhalo).and.(Y0+ii<=ncube+nhalo) ) then
                       !if ( AXC( x0+ii, y0+jj, ip ) < subdis(ii,jj) )  AXC( x0+ii, y0+jj, ip ) = subdis(ii,jj)
                       if ( QC( x0+ii, y0+jj, ip ) <= subq(ii,jj) )  AXC( x0+ii, y0+jj, ip ) = subdis(ii,jj)
                       if ( QC( x0+ii, y0+jj, ip ) <= subq(ii,jj) )  QC( x0+ii, y0+jj, ip )  = subq(ii,jj)
                    endif
                end do
                end do
       end if
 
#ifdef DEBUGRIDGE
       end if
#endif
          end do
       end do
     end do

     write(*,*) " finished paintridge2cube "
     
  end function paintridge2cube


!======================================
subroutine paintridgeoncube ( ncube,nhalo,nsb,nsw , terr  )
   
       integer, intent(in) :: ncube,nhalo,nsb,nsw
       REAL (KIND=dbl_kind), &
            DIMENSION(ncube,ncube,6), INTENT(IN) :: terr

       real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: axc
       logical :: lzero=.false.
!---------------------------------------------------

     mxdsp = mxdis
     where( pkhts < 0)
       mxdsp=0.
     end where
     axc =  paintridge2cube ( mxdsp , ncube,nhalo,nsb,nsw,lzero )

     write(411) ncube
     write(411) axc(1:ncube,1:ncube,:)
     write(411) terr

  end subroutine paintridgeoncube

!====================================
   subroutine remapridge2tiles(area_target,target_center_lon,target_center_lat,  &
         weights_eul_index_all,weights_lgr_index_all,weights_all,ncube,jall,&
         nreconstruction,ntarget,nhalo,nsb)

      use shr_kind_mod, only: r8 => shr_kind_r8
      use remap
      implicit none
      real(r8), intent(in) :: weights_all(jall,nreconstruction)
      integer , intent(in) :: weights_eul_index_all(jall,3),weights_lgr_index_all(jall)
      integer , intent(in) :: ncube,jall,nreconstruction,ntarget,nhalo,nsb
      real(r8), intent(in) :: area_target(ntarget),target_center_lon(ntarget),target_center_lat(ntarget)
      real(r8):: f(ntarget)
  
      REAL  ,                                                          &
         DIMENSION(1-nhalo:ncube+nhalo )                          :: xv,yv,alph,beta
      
      integer :: alloc_error

      integer :: i,ix,iy,ip,ii,counti,norx,nory,i_last,isubr,itile,itili,current_uqrid
      real(r8):: wt
      real(KIND=dbl_kind), dimension(1-nhalo:ncube+nhalo,1-nhalo:ncube+nhalo ,6) :: tmpx6
      real(KIND=dbl_kind), dimension(ncube*ncube*6) :: mxdisC , anglxC, anisoC, hwdthC
      real(KIND=dbl_kind), dimension(ncube*ncube*6) :: mxvrxC , mxvryC, bsvarC, xspkC, yspkC
      integer, dimension(ncube*ncube*6)  :: uqridC 

!----------------------------------------------------------------------------------------------------

        write(*,*) " YOU ASKED FOR TILES !!!!!!!!!!!! "


    DO i=1-nhalo,ncube+nhalo
       xv(i)=1.*i
       yv(i)=1.*i
    END DO

    allocate (uqrid_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for uqrid_tiles'; stop; endif
    uqrid_tiles = -1

    allocate (numbr_tiles(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for numbr_tiles'; stop; endif
    numbr_tiles = 0
    allocate (error_tiles(ntarget),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for error_tiles'; stop; endif
    error_tiles = 0

    allocate (wghts_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for wghts_tiles'; stop; endif
    wghts_tiles = 0.
    allocate (mxdis_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxdis_tiles'; stop; endif
    mxdis_tiles = 0.
    allocate (anglx_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for anglx_tiles'; stop; endif
    anglx_tiles = 0.
    allocate (aniso_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for aniso_tiles'; stop; endif
    aniso_tiles = 0.
    allocate (mxvrx_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxvrx_tiles'; stop; endif
    mxvrx_tiles = 0.
    allocate (mxvry_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for mxvry_tiles'; stop; endif
    mxvry_tiles = 0.
    allocate (bsvar_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for bsvar_tiles'; stop; endif
    bsvar_tiles = 0.
    allocate (hwdth_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for hwdth_tiles'; stop; endif
    hwdth_tiles = 0.
    allocate (agnpk_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for agnpk_tiles'; stop; endif
    agnpk_tiles = 0.
    allocate (bgnpk_tiles(ntarget,ntiles),stat=alloc_error )
    if( alloc_error /= 0 ) then; print*,'Program could not allocate space for bgnpk_tiles'; stop; endif
    bgnpk_tiles = 0.

     norx=size(mxdis,1)
     nory=size(mxdis,2)

     tmpx6 = mapridge2cube ( mxdis , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxdisC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( anglx , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     anglxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( aniso , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     anisoC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( mxvrx , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxvrxC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( mxvry , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     mxvryC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( bsvar , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     bsvarC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( hwdth , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     hwdthC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( xspk , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     xspkC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )
     tmpx6 = mapridge2cube ( yspk , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     yspkC = reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) )


     tmpx6 = mapridge2cube ( 1.*uqrid , norx, nory,xs,ys,xv,yv,ncube,nhalo,nsb )
     uqridC = INT(reshape( tmpx6(1:ncube, 1:ncube, 1:6 ) , (/ncube*ncube*6/) ) )


       !mxdis_target = remap_field(tmp,area_target,weights_eul_index_all(1:jall,:),weights_lgr_index_all(1:jall),&
       !weights_all(1:jall,:),ncube,jall,nreconstruction,ntarget)



    do counti=1,jall
      i   = weights_lgr_index_all(counti)
      ix  = weights_eul_index_all(counti,1)
      iy  = weights_eul_index_all(counti,2)
      ip  = weights_eul_index_all(counti,3)
      !
      ! convert to 1D indexing of cubed-sphere
      !
      ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
      wt = weights_all(counti,1)

      current_uqrid = INT( uqridC( ii ) )
      ! find tile index for current ridge 
      itile=-1
      do itili = 1,ntiles
         if (uqrid_tiles(i,itili) == current_uqrid ) then 
            itile=itili
            numbr_tiles(i)=itili
         endif
      end do   
      ! if itile is still -1 then this is a new ridge for this target
      ! cell
      if (itile == -1 ) then 
      do itili = 1,ntiles
         if (uqrid_tiles(i,itili) == -1) then 
            itile=itili
            uqrid_tiles(i,itili) = current_uqrid
            numbr_tiles(i)=itili
            go to 611
         end if
      end do
      end if
 611  continue
      ! if itile is still -1 then we have more ridges 
      ! than allocated tiles - punt, stuff into last tile
      ! set error=1
      if (itile == -1 ) then
         itile=ntiles 
         error_tiles(i)=1
         numbr_tiles(i)=numbr_tiles(i)+1
      endif
      
      wghts_tiles( i , itile ) = wghts_tiles( i , itile ) + wt
      agnpk_tiles( i , itile ) = xspkC(ii)   
      bgnpk_tiles( i , itile ) = yspkC(ii)   
      !agnxg_tiles( i , itile ) = agnxg_tiles( i , itile ) + wt*agnom(ix)
      !bgnxg_tiles( i , itile ) = bgnxg_tiles( i , itile ) + wt*bgnom(iy)
      bsvar_tiles( i , itile ) = bsvarC(ii)
      mxdis_tiles( i , itile ) = mxdisC(ii)
      mxvrx_tiles( i , itile ) = mxvrxC(ii)
      mxvry_tiles( i , itile ) = mxvryC(ii)
      aniso_tiles( i , itile ) = anisoC(ii)
      hwdth_tiles( i , itile ) = hwdthC(ii)
      anglx_tiles( i , itile ) = anglxC(ii)
      !rwpks_tiles( i , itile ) = rwpks(ii)
      !rwvls_tiles( i , itile ) = rwvls(ii)
      !npks_tiles( i , itile )  = npks(ii)
      !nvls_tiles( i , itile )  = nvls(ii)

    enddo


       write(*,*) " Max Tiles Needed ",maxval( numbr_tiles ) 

    where( wghts_tiles > 0.)
     !agnxg_tiles = agnxg_tiles/wghts_tiles
     !bgnxg_tiles = bgnxg_tiles/wghts_tiles
    elsewhere
                                 !where( wghts_tiles <= 0.)
     agnpk_tiles = -9999._r8
     bgnpk_tiles = -9999._r8
     !agnxg_tiles = -9999._r8
     !bgnxg_tiles = -9999._r8
    end where

    !OPEN (unit = 611, file= trim(tfile$) ,form="UNFORMATTED" )
    OPEN (unit = 611, file= 'RidgeTile.dat' ,form="UNFORMATTED" )
    write( 611 ) ntarget,ntiles
    write( 611 ) error_tiles
    write( 611 ) uqrid_tiles
    write( 611 ) wghts_tiles
    write( 611 ) mxdis_tiles
    write( 611 ) anglx_tiles
    write( 611 ) agnpk_tiles
    write( 611 ) bgnpk_tiles
    !write( 611 ) agnxg_tiles
    !write( 611 ) bgnxg_tiles
    write( 611 ) mxvrx_tiles
    write( 611 ) mxvry_tiles
    write( 611 ) aniso_tiles
    write( 611 ) hwdth_tiles



      write(*,*) " GOT OUT OF remapridge2tiles OK "

 
   end subroutine remapridge2tiles

!==================================================================


end module ridge_ana
