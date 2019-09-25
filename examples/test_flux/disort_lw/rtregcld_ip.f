C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/rtregcld.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 3.4 $
C     created:   $Date: 2010/07/07 21:10:53 $
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER). |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------

      SUBROUTINE RTREGCLDIP(ncol, nlay, ngpt,  NUMANG, optProp, 
     &   sfc_emis, planckFrac, planckLay, planckLev,
     &   planckSfc, 
     &   lay_source, lev_source_inc, lev_source_dec,
     &   TOTUFLUXout, TOTDFLUXout, top_at_1)
      use mo_rte_kind,           only: wp
      use mo_optical_props,       only: ty_optical_props_arry, 
     &    ty_optical_props_1scl, ty_optical_props_2str
      use mo_optical_props_add,   only: ty_optical_props_tip

C       use DISORT_LW_mod,          only: DISORT_LW  
      IMPLICIT NONE
      real,    PARAMETER :: TBLINT = 10000.0
      integer, PARAMETER :: NTBL = 10000
      integer, PARAMETER :: MXLAY=603
      integer, PARAMETER :: MG = 16
      integer, PARAMETER :: NBANDS = 16
      integer, PARAMETER :: MXANG = 4
      integer, PARAMETER :: MCMU = 32, MUMU = 32, MPHI = 3
      integer, PARAMETER :: MXSTR = 16

      integer, intent(in) :: ncol, nlay, ngpt
      integer, intent(in) :: NUMANG

      class(ty_optical_props_arry),    intent(in) :: optProp
      real(wp), dimension(NBANDS,ncol),intent(in) :: sfc_emis  ! block_size, nblocks (emissivity is spectrally constant)
      real(wp), dimension(ncol,nlay  ,ngpt),intent(in) :: planckFrac  ! block_size, nblocks (emissivity is spectrally constant)
      real(wp), dimension(ncol,nlay+1,ngpt),intent(in) :: planckLev
      real(wp), dimension(ncol,nlay,  ngpt),intent(in) :: planckLay
      real(wp), dimension(ncol,       ngpt),intent(in) :: planckSfc

      real(wp), dimension(ncol,nlay,  ngpt),intent(in) :: lev_source_inc
      real(wp), dimension(ncol,nlay,  ngpt),intent(in) :: lev_source_dec
      real(wp), dimension(ncol,nlay,  ngpt),intent(in) :: lay_source

      REAL(wp), DIMENSION(0:nlay,ncol),intent(inout)   :: TOTDFLUXout 
      REAL(wp), DIMENSION(0:nlay,ncol),intent(inout)   :: TOTUFLUXout    
      LOGICAL                                          :: top_at_1

C *** This program calculates the upward fluxes, downward fluxes, and
C     heating rates for an arbitrary cloudy atmosphere.  The input
C     to this program is the atmospheric profile, including cloud
C     properties,  and all needed Planck function information.  First  
C     order standard Gaussian quadrature is used for the angle 
C     integration.

C     Clouds are treated with random overlap scheme.


C       IMPLICIT DOUBLE PRECISION (V)    

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
      REAL, PARAMETER :: FLUXFAC = PI * 2.

      REAL, DIMENSION(0:MXLAY)      :: TZ, PZ 
      REAL, DIMENSION(NBANDS)       :: SEMISS 
      REAL, DIMENSION(NBANDS)       :: WAVENUM1, WAVENUM2
      REAL   :: DELWAVE
      REAL, DIMENSION(MXLAY)        :: TAUG
      REAL, DIMENSION(MXLAY)        :: FRACS
      REAL, DIMENSION(NBANDS)       :: NG, NSPA, NSPB
      REAL                          :: PLNKBND
      REAL, DIMENSION(MXLAY)        :: PLNKLAY
      REAL, DIMENSION(0:MXLAY)      :: PLNKLEV
      REAL ::  WAVENUMHI
      REAL ::  WAVENUMLO


      REAL, DIMENSION(MXLAY)        :: tauLay
      REAL, DIMENSION(MXLAY)        :: tranLay
      REAL, DIMENSION(MXLAY)        :: source_dn, source_up
      REAL, DIMENSION(MXLAY)        :: lev_dn, lev_up

      INTEGER         :: IBAND
      INTEGER         :: ICLD, ICOL, ITR
      INTEGER         :: ISTART
      INTEGER         :: IEND
      INTEGER         :: IG
      INTEGER         :: IQ
      INTEGER         :: IREFLECT=0
      INTEGER         :: LAY, LEV, K
      INTEGER         :: NLAYERS
      INTEGER         :: NCBANDS=-99

      INTEGER         :: iScheme = 4


      REAL, DIMENSION(MXLAY,NBANDS) :: TAUCLOUD 
      REAL, DIMENSION(MXLAY) :: SSACLOUD 
      REAL, DIMENSION(0:16,MXLAY,NBANDS) :: XMOM 

      REAL, DIMENSION(0:MXLAY)      :: FNET 
      REAL, DIMENSION(0:MXLAY)      :: HTR 
      REAL, DIMENSION(MXLAY)        :: CLDFRAC 

C       COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
C       COMMON /CONTROL/  NUMANGS, ISCAT, NSTR, 
C      &                  IOUT, ISTART, IEND, ICLD


C       COMMON /PLANKG/    FRACS(MXLAY,MG)                                       
C       COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      REAL :: TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY)

      REAL :: BPADE, TAUTBL(0:NTBL),TRANS(0:NTBL),TF(0:NTBL)

      REAL:: ATRANS(MXLAY,MXANG),BBUGAS(MXLAY,MXANG)
      REAL:: BBDGAS(MXLAY,MXANG)
      Real:: ATOT(MXLAY,MXANG),ODCLD(MXLAY,NBANDS,MXANG)
      Real:: UFLUX(0:MXLAY),DFLUX(0:MXLAY),BBUTOT(MXLAY,MXANG)
      Real:: DRAD(0:MXLAY-1,MXANG),URAD(0:MXLAY,MXANG)
      Real:: dradg(0:mxlay,mxang),uradg(0:mxlay,mxang) !tang
      Real:: SECANG(MXANG),ANGWEIGH(MXANG),RAD(MXANG)
      Real:: SECREG(MXANG,MXANG),WTREG(MXANG,MXANG)
      Real:: EFCLFRAC(MXLAY,NBANDS,MXANG)
      Real:: ABSCLD(MXLAY,NBANDS,MXANG)

      real :: aaa, bbb, ccc, BBD, BBDTOT, BLAY, dcoefff, ttt
      real :: DPLANKDN, DPLANKUP, GASSRC
      real :: TAUSFAC, TBLIND, TFACGAS, TFN, TRANSC, TRANSCLD

      integer :: IANG, IB, ITGAS, ITTOT, jj

C Dimensions for cloud 
      real :: ODEPTH, ODEPTH_REC, ODTOT, ODTOT_REC, PADE, PLFRAC, RAD0
      real :: RADLD, RADLD1, RADLU, RADSUM, REFLECT, TFACTOT

      real :: xx, yy, xx1, xx2, xx3


C *** When standard first-order Gaussian quadrature is chosen as
C     the method to approximate the integral over angles that yields
C     flux from radiances, then SECREG(I,J) is the secant of the Ith  
C     (out of a total of J angles) and WTREG(I,J) is the corresponding
C     weight.
      DATA SECREG(1,1) / 1.5/
      DATA SECREG(2,2) / 2.81649655/, SECREG(1,2) / 1.18350343/
      DATA SECREG(3,3) / 4.70941630/, SECREG(2,3) / 1.69338507/
      DATA SECREG(1,3) / 1.09719858/
      DATA SECREG(4,4) / 7.15513024/, SECREG(3,4) / 2.40148179/
      DATA SECREG(2,4) / 1.38282560/, SECREG(1,4) / 1.06056257/
      DATA WTREG(1,1) / 0.50000000/
      DATA WTREG(2,2) /0.1819586183/, WTREG(1,2) /0.3180413817/
      DATA WTREG(3,3) /0.0698269799/, WTREG(2,3) /0.2292411064/
      DATA WTREG(1,3) /0.2009319137/
      DATA WTREG(4,4) /0.0311809710/, WTREG(3,4) /0.1298475476/
      DATA WTREG(2,4) /0.2034645680/, WTREG(1,4) /0.1355069134/
      real, parameter :: REC_6 = 0.166667

      RADSUM = 0.
      CLDFRAC = 1.0
!/tang
      if (iScheme.eq.0.or.iScheme.eq.1.or.iScheme.eq.3) then
         dcoefff = 0.0
      elseif (iScheme.eq.2) then
         dcoefff = 0.3
      elseif (iScheme.eq.4) then
         dcoefff = 0.4
      endif
       dcoefff = 0.5
!\tang
C *** Load angle data in arrays depending on angular quadrature scheme.
      DO 100 IANG = 1, NUMANG
         SECANG(IANG)   = SECREG(IANG,NUMANG)
         ANGWEIGH(IANG) = WTREG(IANG,NUMANG)
 100  CONTINUE

C  Compute lookup tables for transmittance, tau transition function,
C  and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
C  computed as a function of the tau transition function, transmittance 
C  is calculated as a function of tau, and the tau transition function 
C  is calculated using the linear in tau formulation at values of tau 
C  above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
C  are computed at intervals of 0.001.  The inverse of the constant used
C  in the Pade approximation to the tau transition function is set to b.
C  These values are not necessary when using DISORT as the RT solver.

            TAUTBL(0) = 0.0
            TAUTBL(NTBL) = 1.E10
            TRANS(0) = 1.0
            TRANS(NTBL) = 0.0
            TF(0) = 0.0
            TF(NTBL) = 1.0
            PADE  = 0.278
            BPADE = 1.0/PADE
            DO 500 ITR = 1,NTBL-1
               TFN = ITR/FLOAT(NTBL)
               TAUTBL(ITR) = BPADE*TFN/(1.-TFN)
               TRANS(ITR) = EXP(-TAUTBL(ITR))
               IF (TAUTBL(ITR) .LT. 0.06) THEN
                  TF(ITR) = TAUTBL(ITR)/6.
               ELSE
                  TF(ITR) = 1.-
     &                 2.*((1./TAUTBL(ITR))-
     &                 (TRANS(ITR)/(1.-TRANS(ITR))))
               ENDIF
 500        CONTINUE


C *** Loop over frequency bands.
      ISTART=1
      IEND = 16
      NLAYERS = nlay

      DO ICOL=1, ncol
        TOTUFLUX(0:NLAYERS) = 0.0
        TOTDFLUX(0:NLAYERS) = 0.0
        URAD(0:NLAYERS,1:NUMANG) = 0.
        DRAD(0:NLAYERS,1:NUMANG) = 0.

C *** Loop over frequency bands.
          RADSUM=0.

      DO 6000 IBAND = ISTART, IEND
          IB=IBAND
         WAVENUMLO = optProp%band_lims_wvn(1,IBAND)
         WAVENUMHI = optProp%band_lims_wvn(2,IBAND)
         DELWAVE = abs(1e4/WAVENUMHI-1e4/WAVENUMLO)

         
C ***    Loop over g-channels.
         DO IG= optProp%band2gpt(1,IBAND), optProp%band2gpt(2,IBAND)
         PLNKBND = planckSfc(ICOL, IG)
         if (top_at_1) then
           PLNKLAY(1:nlay)=planckLay(ICOL,  nlay:1-1, IG)
           PLNKLEV(0:nlay)=planckLev(ICOL,  nlay+1:1-1, IG)
           FRACS(1:nlay) = planckFrac(ICOL, nlay:1-1, IG)
         else
           PLNKLAY(1:nlay)=planckLay(ICOL, 1:nlay, IG)
           PLNKLEV(0:nlay)=planckLev(ICOL, 1:nlay+1, IG)
           FRACS(1:nlay) = planckFrac(ICOL, 1:nlay, IG)
         endif  
         DO LAY = 1, NLAYERS

              if (top_at_1)  then
                jj=NLAYERS-LAY+1
              else
                jj=LAY
              endif
              ODEPTH=optProp%tau(ICOL,jj,IG)

              TAUG(LAY)=ODEPTH
                !  layers
              select type (optProp)
                type is (ty_optical_props_1scl) ! two-stream calculation
                    SSACLOUD(LAY)=0.
                    ODEPTH=0.

                type is (ty_optical_props_2str) ! two-stream calculation
                    SSACLOUD(LAY)=optProp%ssa(ICOL,jj,IG)
                    ODEPTH=ODEPTH*optProp%ssa(ICOL,jj,IG)
        
                type is (ty_optical_props_tip)
                    SSACLOUD(LAY)=optProp%ssa(ICOL,jj,IG)
                    ODEPTH=ODEPTH*optProp%ssa(ICOL,jj,IG)
              end select
              IF (ODEPTH < 1.E-6) SSACLOUD(LAY)=0.
          ENDDO


         
C ***    Loop over each angle for which the radiance is to be computed.
         DO 3000 IANG = 1, NUMANG
C ***       Radiative transfer starts here.
            RADLD = 0.
            dradg(NLAYERS,iang) = radld

C ***       Downward radiative transfer.  
            DO 2500 LEV = NLAYERS, 1, -1
               BLAY     = PLNKLAY(LEV)
               ODEPTH = SECANG(IANG) * TAUG(LEV)
               IF (ODEPTH .LT. 0.0) ODEPTH = 0.0
               call getTFandT(ODEPTH, ATRANS(LEV,IANG),TAUSFAC,TRANS,TF)
               PLFRAC = FRACS(LEV)*ATRANS(LEV,IANG)

               BBD = PLFRAC*(BLAY+(PLNKLEV(LEV-1) - BLAY)*TAUSFAC)

               BBDGAS(LEV,IANG) = BBD

               BBUGAS(LEV,IANG) = PLFRAC*
     &                    (BLAY+(PLNKLEV(LEV)   - BLAY)*TAUSFAC)

               RADLD = RADLD*(1.-ATRANS(LEV,IANG)) + BBD
               dradg(lev-1,iang) = radld !tang
 2500       CONTINUE
            RAD(IANG) = RADLD
            RADSUM = RADSUM + ANGWEIGH(IANG) * RADLD
 3000    CONTINUE 

         RAD0 = FRACS(1) * PLNKBND*sfc_emis(IBAND, ICOL)
         REFLECT = 1. - sfc_emis(IBAND, ICOL)
         DO 4000 IANG = 1, NUMANG
C           Add in reflection of surface downward radiance.
            IF (IREFLECT .EQ. 1) THEN
C              Specular reflection.
               RADLU = RAD0 + REFLECT * RAD(IANG)
            ELSE
C              Lambertian reflection.
               RADLU = RAD0 + 2. * REFLECT * RADSUM
            ENDIF
            uradg(0,iang) = radlu !tang
            URAD(0,IANG) = URAD(0,IANG) + RADLU

C ***       Upward radiative transfer.
            DO 2600 LEV = 1, NLAYERS
               RADLU = RADLU*(1.-ATRANS(LEV,IANG)) + BBUGAS(LEV,IANG)
               radld1 = ssacloud(lev)
               if (radld1 > epsilon(1.)) then

                 radld1 = dcoefff*radld1*(
     &             dradg(lev,iang)*(1.-(1.-ATRANS(LEV,IANG))**2)
     &            -BBDGAS(LEV,IANG)*(1.-ATRANS(LEV,IANG))
     &            -BBUGAS(LEV,IANG))
              endif
              uradg(lev,iang) = radlu !tang
              RADLU = RADLU + radld1
              URAD(LEV,IANG) = URAD(LEV,IANG) + RADLU  
 2600       CONTINUE
 4000    CONTINUE 
!/tang
         DO 4001 IANG = 1, NUMANG
            radld = 0.
            DO 2601 LEV = NLAYERS,1,-1
               RADLD = RADLD*(1.-ATRANS(LEV,IANG)) +  BBDGAS(LEV,IANG)
               radld1 = ssacloud(lev)
               if (radld1 > epsilon(1.)) then
                  radld1 = dcoefff*radld1*(
     &              uradg(lev-1,iang)*(1.-(1.-ATRANS(LEV,IANG))**2)
     &              -BBUGAS(LEV,IANG)*(1.-ATRANS(LEV,IANG))
     &           -BBDGAS(LEV,IANG))
              endif    
              RADLD = RADLD + radld1
              DRAD(LEV-1,IANG) = DRAD(LEV-1,IANG) + RADLD
 2601       CONTINUE
 4001    CONTINUE
!\tang
         RADSUM = 0.
         ENDDO

C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = SUM(URAD(LEV,1:NUMANG)*ANGWEIGH(1:NUMANG))
            DFLUX(LEV) = SUM(DRAD(LEV,1:NUMANG)*ANGWEIGH(1:NUMANG))
            URAD(LEV,1:NUMANG) = 0.
            DRAD(LEV,1:NUMANG) = 0.
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + UFLUX(LEV) *FLUXFAC
            TOTDFLUX(LEV) = TOTDFLUX(LEV) + DFLUX(LEV) *FLUXFAC
 5000    CONTINUE
 6000 CONTINUE 
      TOTUFLUXout(0:NLAYERS, ICOL) = TOTUFLUX(0:NLAYERS)
      TOTDFLUXout(0:NLAYERS, ICOL) = TOTDFLUX(0:NLAYERS)
      ENDDO   
C       DO LEV = 0, NLAYERS
C          TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
C          TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
C          FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
C       ENDDO   

C C *** Calculate Heating Rates.
C       HTR(NLAYERS) = 0.0
C       DO LEV  = NLAYERS-1, 0, -1
C          LAY = LEV + 1
C          HTR(LEV)=HEATFAC*(FNET(LEV)-FNET(LAY))/(PZ(LEV)-PZ(LAY)) 
C       ENDDO   

      RETURN
      END   

  