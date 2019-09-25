C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/rtrdis.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 5.4 $
C     created:   $Date: 2010/07/07 21:10:53 $

      SUBROUTINE RTRDIS(ncol, nlay, ngpt,  optProp, 
     &   t_lev, sfc_t, sfc_emis, planckFrac, planckLev,
     &   TOTUFLUX, TOTDFLUX, top_at_1)
      use mo_rte_kind,           only: wp
      use mo_optical_props,       only: ty_optical_props_arry, 
     &    ty_optical_props_1scl, ty_optical_props_2str
C       use DISORT_LW_mod,          only: DISORT_LW  
      IMPLICIT NONE
      integer, PARAMETER :: MXLAY=603
      integer, PARAMETER :: MG = 16
      integer, PARAMETER :: NBANDS = 16
      integer, PARAMETER :: MXANG = 4
      integer, PARAMETER :: MCMU = 32, MUMU = 32, MPHI = 3
      integer, PARAMETER :: MXSTR = 16

      integer, intent(in) :: ncol, nlay, ngpt
      class(ty_optical_props_arry),    intent(in) :: optProp
      real(wp), dimension(ncol,nlay+1),intent(in) :: t_lev
      real(wp), dimension(ncol),       intent(in) :: sfc_t  ! block_size, nblocks (emissivity is spectrally constant)
      real(wp), dimension(NBANDS,ncol),intent(in) :: sfc_emis  ! block_size, nblocks (emissivity is spectrally constant)
      real(wp), dimension(ncol,nlay,ngpt),intent(in) :: planckFrac  ! block_size, nblocks (emissivity is spectrally constant)
      real(wp), dimension(ncol,nlay+1,ngpt),intent(in) :: planckLev

      REAL(wp), DIMENSION(0:nlay,ncol),intent(inout)   :: TOTDFLUX 
      REAL(wp), DIMENSION(0:nlay,ncol),intent(inout)   :: TOTUFLUX 
      LOGICAL                                          :: top_at_1
C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and all Planck function
C     information.  First-order "numerical" quadrature is used for the 
C     angle integration, i.e. only one exponential is computed per layer
C     per g-value per band.

      REAL, DIMENSION(MXLAY,NBANDS) :: TAUCLOUD 
      REAL, DIMENSION(MXLAY,NBANDS) :: SSACLOUD 
      REAL, DIMENSION(0:16,MXLAY,NBANDS) :: XMOM 

      REAL, DIMENSION(0:MXLAY)      :: FNET 
      REAL, DIMENSION(0:MXLAY)      :: HTR 

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
C     (ORIGINAL RRTM VALUE) DATA HEATFAC /8.4391/

C     Calculated value:
C     (grav) x (#sec/day) / (specific heat of dry air at const. p x 1.e2)
C     Here, cpdair is in units of J g-1 K-1; grav is cm sec-2; 
C     and a constant (1.e2) converts mb to Pa when heatfac 
C     is multiplied by W m-2 mb-1. 
      REAL, PARAMETER :: PLANCK = 6.62606876E-27 
      REAL, PARAMETER :: BOLTZ  = 1.3806503E-16 
      REAL, PARAMETER :: CLIGHT = 2.99792458E+10 
      REAL, PARAMETER :: AVOGAD = 6.02214199E+23 
      REAL, PARAMETER :: ALOSMT = 2.6867775E+19 
      REAL, PARAMETER :: GASCON = 8.314472E+07
      REAL, PARAMETER :: RADCN1 = 1.191042722E-12
      REAL, PARAMETER :: RADCN2 = 1.4387752    
      REAL, PARAMETER :: GRAV   = 9.80665E+02
      REAL, PARAMETER :: CPDAIR = 1.00464
      REAL, PARAMETER :: AIRMWT = 28.964  
      REAL, PARAMETER :: SECDY  = 8.64E+04

      REAL, PARAMETER :: PI = ACOS(-1.0)
      REAL, PARAMETER :: HEATFAC = 1.0E-7*(GRAV * SECDY)/(CPDAIR)
      REAL, PARAMETER :: ONEMINUS = 1. - 1.E-6
      REAL, PARAMETER :: FLUXFAC = PI * 2.E4  
      INTEGER         :: IBAND
      INTEGER         :: ICLD
      INTEGER         :: ISTART
      INTEGER         :: IEND
      INTEGER         :: IG
      INTEGER         :: IQ
      INTEGER         :: IREFLECT
      INTEGER         :: LAY, LEV, K
      INTEGER         :: NLAYERS

      REAL, DIMENSION(0:MXLAY)      :: TZ, PZ 
      REAL, DIMENSION(NBANDS)       :: SEMISS 
      REAL, DIMENSION(NBANDS)       :: WAVENUM1, WAVENUM2, DELWAVE
      REAL, DIMENSION(MXLAY,MG)     :: FRACS, TAUG
      REAL, DIMENSION(NBANDS)       :: NG, NSPA, NSPB
      REAL, DIMENSION(MXLAY)        :: PLANKLAY
      REAL                          :: PLANKBND

      CHARACTER(len=18)  ::   HNAMRDS,HVRRDS
      CHARACTER(len=127) ::   HEADER

                                       
      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, 
     &          NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0 
      LOGICAL   PRNT( 7 )
      REAL      ALBMED( MUMU ), DFDT( MXLAY ), TAUREV(MXLAY),
     &          FLUP( MXLAY ), HL( 0:MCMU ), PHI( MPHI ),
     &          PMOM( 0:MCMU, MXLAY ), RFLDIR( MXLAY ),
     &          RFLDN( MXLAY ), SSALB( MXLAY ), TEMPER( 0:MXLAY ),
     &          TRNMED( MUMU ), U0U( MUMU, MXLAY ), UAVG( MXLAY ),
     &          UMU( MUMU ), UTAU( MXLAY ),
     &          UU( MUMU, MXLAY, MPHI ),fldir(mxlay),fldn(mxlay), 
     &          TZREV(0:MXLAY),FRACSREV(MXLAY)

      REAL ::  PLNKLEV(MXLAY) 
      REAL ::  PHASERAY(0:MXSTR)
      REAL ::  gasm
      REAL ::  TBOUND
      REAL ::  WAVENUMHI
      REAL ::  WAVENUMLO
      INTEGER :: jj, ICOL
      character(len=100) :: fname
      integer, parameter :: fid = 11000


      DATA PRNT /.false.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,
     &     .FALSE.,.FALSE./

c     Ensure all cloud properties are equal to 0.0
      IF (ICLD .EQ. 0) THEN
         TAUCLOUD(:,:) = 0.0
         SSACLOUD(:,:) = 0.0
         XMOM(:,:,:) = 0.0
      ENDIF
         
      HVRRDS = '$Revision: 5.4 $'
      HEADER = ''
C       convertion to local variables
      NSTR = 4*2
      print *,'NSTR:', NSTR
      NLAYERS=nlay
! TOP layer
      TEMIS = 0.0
! bottom layer

      USRTAU = .FALSE.
      USRANG = .FALSE.
      NPHI = 0
      IBCND = 0
      PHI0 = 0.
      UMU0 = 0.0
      FBEAM = 0.0
      FISOT = 0.0

C     surface reflectivity
      IREFLECT = 0

      IF (IREFLECT .EQ. 0) THEN
         LAMBER = .TRUE.
      ELSE
         LAMBER = .FALSE.
      ENDIF

      PLANK = .TRUE.
      ONLYFL = .TRUE.
      ACCUR = 0.0001
      MAXCLY = MXLAY
      MAXULV = MXLAY
      MAXUMU = MUMU
      MAXCMU = MCMU
      MAXPHI = MPHI
      TOTDFLUX=0.
      TOTUFLUX=0.
      SSALB(1:NLAYERS) = 0.
      PMOM(0:NSTR,1:NLAYERS) = 0.

C *** Loop over frequency bands.
      ISTART=1
      IEND = 16
      DO ICOL=1,2!ncol
      BTEMP  = sfc_t(ICOL)
      TTEMP  = t_lev(ICOL, nlay+1)
      if (top_at_1)  then
        TZREV(0:NLAYERS) = t_lev(ICOL, 1:NLAYERS+1)
      else
        DO LAY = 0, NLAYERS
           TZREV(NLAYERS-LAY) = t_lev(ICOL, LAY+1)
        enddo
      endif
      DO 6000 IBAND = ISTART, IEND
        

c  set albedo for this band
         ALBEDO = 1. - sfc_emis(IBAND,ICOL)
         WAVENUMLO = optProp%band_lims_wvn(1,IBAND)
         WAVENUMHI = optProp%band_lims_wvn(2,IBAND)

C          if (iband .eq. 16 .and. istart .ne. 16) 
C      &        wavenumhi = 5000.
C ***    Loop over g-channels.
         DO IG= optProp%band2gpt(1,IBAND), optProp%band2gpt(2,IBAND)
C          print *, 'BAND', IBAND, 'IG', IG , wavenumlo,wavenumhi

C ***    Downward radiative transfer.

        !  layers
          select type (optProp)
            type is (ty_optical_props_1scl) ! two-stream calculation
        !  layers
              if (top_at_1)  then
                TAUREV(1:nlay) = optProp%tau(ICOL,1:nlay,IG)
                SSALB (1:nlay) = 0.
              else  
                TAUREV(nlay:1:-1) = optProp%tau(ICOL,1:nlay,IG)
                SSALB(nlay:1:-1)  = 0.
              endif
              PMOM=0.
              PMOM(0,:)=1.
            type is (ty_optical_props_2str) ! two-stream calculation

              if (top_at_1)  then
                TAUREV(1:nlay) = optProp%tau(ICOL,1:nlay,IG)
                SSALB (1:nlay) = optProp%ssa(ICOL,1:nlay,IG)
                do jj = 1, nlay
                  gasm = optProp%g(ICOL,jj,IG)
                  CALL  GETMOMloc( 3, gasm, NSTR, PMOM(:,jj) )
                enddo  
              else  
                TAUREV(nlay:1:-1) = optProp%tau(ICOL,1:nlay,IG)
                SSALB(nlay:1:-1)  = optProp%ssa(ICOL,1:nlay,IG)
                do jj = 1, nlay
                  gasm = optProp%g(ICOL,jj,IG)
                  CALL  GETMOMloc( 3, gasm, NSTR, PMOM(:,nlay-jj+1) )
                enddo  
              endif

            class default
                print *,'TYPE IS NOT SUPPORTED'
                call exit(200)
          end select

          NTAU = nlay+1
          if (top_at_1)  then
            FRACSREV(1:nlay) = planckFrac(ICOL, 1:nlay, IG)
            PLNKLEV (1:nlay) = planckLev(ICOL, 1:NTAU, IG)
          else  
            FRACSREV(nlay:1:-1) = planckFrac(ICOL, 1:nlay, IG)
            PLNKLEV (NTAU:1:-1) = planckLev(ICOL, 1:NTAU, IG)
          endif

         if (IG==1 .and. ICOL==1) THEN
          select type (optProp)
            type is (ty_optical_props_1scl) ! two-stream calculation
            print *,'clear'
          type is (ty_optical_props_2str) ! two-stream calculation
              print *,'2str'
          end select

            do lay=1, NLAYERS
                print *, lay, TAUREV(LAY), SSALB(lay), PMOM(1,LAY)
            enddo
         endif

C       do jj=1, NLAY
C           print *, FRACSREV(jj), TAUREV(jj), SSALB(jj), PMOM(0:2, jj)
C       enddo

         CALL DISORT_LW( NLAYERS, FRACSREV, PLNKLEV,
     &        TAUREV, SSALB, NSTR, PMOM, 
     &        TZREV, wavenumlo,wavenumhi,
     &        USRTAU, NTAU, UTAU, NSTR, 
     &        USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, 
     &        UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP,
     &        TEMIS, PLANK, ONLYFL, ACCUR, PRNT, 
     &        HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXCMU,
     &        RFLDIR, RFLDN, FLUP, 
     &        DFDT, UAVG, UU, ALBMED, TRNMED)

         DO LEV = NLAYERS, 0, -1
          TOTUFLUX(LEV,ICOL) = TOTUFLUX(LEV,ICOL) + FLUP(NLAYERS-LEV+1)
          TOTDFLUX(LEV,ICOL) = TOTDFLUX(LEV,ICOL)+RFLDIR(NLAYERS-LEV+1)
     &             + RFLDN(NLAYERS-LEV+1) 
         ENDDO

         IF (RFLDIR(1) .GT. 1.e-5 .OR. RFLDN(1) .GT. 1.e-5) 
     &        WRITE(*,9000) IBAND, IG, RFLDN(1)

      ENDDO
            
 6000 CONTINUE
      ENDDO
C       FNET(0:NLAYERS) = TOTUFLUX(0:NLAYERS) - TOTDFLUX(0:NLAYERS)
C C       WRITE(*,'(//,A,/,A,/,A)') 
C C      &  '                  <-------------- FLUXES -------------->', 
C C      &  '    Optical       Downward       Downward         Upward'// 
C C      &  '    d(Net Flux)', 
C C      &  '      Depth         Direct        Diffuse        Diffuse'// 
C C      &  '    / d(Op Dep)'
C C        DO  JJ = 0, NLAYERS

C C        WRITE(*,'(0P,26X,1P,4E15.4)') 
C C      &          TOTDFLUX(JJ), TOTUFLUX(JJ), FNET(JJ)

C C        ENDDO

C       HTR(NLAYERS) = 0.
C       DO LEV = NLAYERS-1, 0, -1
C          HTR(LEV) = HEATFAC * (FNET(LEV) -FNET(LEV+1)) /
C      &        (PZ(LEV) - PZ(LEV+1))
C       ENDDO

 9000 FORMAT('DOWNWARD FLUX AT TOA GTR THAN 0. IN BAND ',i2,
     & 'AT IG =',i2,'. POSSIBLE',/,
     &'INSTABILITY IN DISORT, TRY INCREASING NUMBER OF STREAMS.',
     &     e15.7)

      RETURN
      END   

      SUBROUTINE  GETMOMloc( IPHAS, GG, NMOM, PMOM )

!        Calculate phase function Legendre expansion coefficients
!        in various special cases


!       INPUT: IPHAS   Phase function options
!                      1 : Isotropic
!                      2 : Rayleigh
!                      3 : Henyey-Greenstein with asymmetry factor GG
!                      4 : Haze L as specified by Garcia/Siewert
!                      5 : Cloud C.1 as specified by Garcia/Siewert

!              GG      Asymmetry factor for Henyey-Greenstein case

!              NMOM    Index of highest Legendre coefficient needed
!                        ( = number of streams 'NSTR'  chosen
!                         for the discrete ordinate method)

!      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
!                         (be sure to dimension '0:maxval' in calling
!                          program)

!      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results
!                     in Radiative Transfer, Transp. Theory and Stat.
!                     Physics 14, 437-484, Tables 10 And 17
! ------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   IPHAS, NMOM
      REAL      GG
!     ..
!     .. Array Arguments ..

      REAL      PMOM( 0:NMOM )
!     ..
!     .. Local Scalars ..

      INTEGER   K
!     ..
!     .. Local Arrays ..

      REAL      CLDMOM( 299 ), HAZELM( 82 )
!     ..
!     .. External Subroutines ..

      EXTERNAL  ERRMSG
!     ..
!     .. INTRINSIC  Functions ..

      INTRINSIC  MIN
!     ..

   
      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
         PMOM(K) = 0.0
   10 CONTINUE


      IF ( IPHAS.EQ.2 )  THEN
!                                       ** Rayleigh phase function
         PMOM(2) = 0.1

      ELSE IF ( IPHAS.EQ.3 ) THEN
!                                       ** Henyey-Greenstein phase fcn
         DO  20  K = 1, NMOM
            PMOM(K) = GG**K
   20    CONTINUE

      END IF

      END SUBROUTINE  GETMOMloc
