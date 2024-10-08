MODULE trcopt_blingv0

   !!======================================================================
   !!                     ***  MODULE trcopt_bling  ***
   !! TOP :  BLINGv0 Compute the light availability in the water column
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   !! Available light calculation in BLINGv0
   !!-----------------------------------------------------------------------
   !! There are multiple types of light.
   !!   IRR_INST is the instantaneous irradiance field.
   !!   IRR_MIX  is the same, but with the irr_inst averaged throughout the  
   !! mixed layer (turbocline*) to account for mixing directly below the boundary 
   !! layer. This quantity is intended to represent the light to which phytoplankton 
   !! subject to turbulent transport in the mixed-layer would be exposed.
   !!   IRR_MEM  is a temporally smoothed field carried between timesteps, to 
   !! represent photoadaptation.
   !!-----------------------------------------------------------------------
   !! *The turbocline depth is the depth at which the
   !! vertical eddy diffusivity coefficient (resulting from the vertical physics
   !! alone, not the isopycnal part, see trazdf.F) fall below a given value
   !! defined locally (avt_c here taken equal to 5 cm/s2). Check subroutine
   !! OPA_SRC/ZDF/zdfmxl.F90 for further explanation.
   !!-----------------------------------------------------------------------

   USE oce_trc
   USE trc
!   USE prtctl_trc      ! Print control for debbuging
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_opt_bling_rb  ! called by trc_sms_bling
   PUBLIC trc_opt_bling_rgb ! called by trc_sms_bling
   PUBLIC trc_opt_bling_init

   INTEGER  :: nksrp
   REAL(wp) :: parlux = 0.43_wp / 3._wp
   REAL(wp), DIMENSION(3,61), PUBLIC ::   xkrgb  

#include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_opt_bling_rb (kt, Kmm, irr_inst, irr_mix)

      INTEGER , INTENT(in) :: kt
      INTEGER , INTENT(in) ::   Kmm   ! time level index
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: irr_mix, irr_inst

      INTEGER  :: ji, jj, jk, kblt          
      REAL(wp) :: zkr, zparr, zcoef, parfc, sumz_irrad_ML, sumz_hblt
      
!ALLOCATE ( irr_inst(jpi,jpj,jpk) )
!      ALLOCATE ( irr_mix(jpi,jpj,jpk) )

      IF( ln_timing == 1 )  CALL timing_start('trc_opt_bling_rb')

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_opt_bling : BLING optic-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~ '
      ENDIF

      ! light attenuation by water and phytoplankton
      zcoef=12*c2n/rcchl/rpig
      parfc=0.43d0/2.d0

      DO ji=1, jpi
         DO jj=1, jpj
               !One-wavelength biopt
               !!! AGT fix here!!!
               zkr   = xkr0 + xkrp * chl_bling(ji,jj,1)
               zparr = qsr(ji,jj)*0.43d0 
               zparr = zparr * EXP(-zkr*0.5d0*e3t(ji,jj,1,Kmm))
               irr_inst(ji,jj,1)=zparr
               irr_mix (ji,jj,1)=zparr               
               kblt=1
               !Irradiance at the surface
               sumz_irrad_ML=irr_mix(ji,jj,1)*e3t(ji,jj,1,Kmm)
               sumz_hblt    =e3t(ji,jj,1,Kmm)
               
               DO jk=2, jpk
               
                  zkr   = xkr0 + xkrp * chl_bling(ji,jj,jk-1)
                  zparr = zparr * EXP(-zkr*0.5d0*e3t(ji,jj,jk-1,Kmm))
                  zkr   = xkr0 + xkrp * chl_bling(ji,jj,jk)
                  zparr = zparr * EXP(-zkr*0.5d0*e3t(ji,jj,jk,Kmm))
                  irr_inst(ji,jj,jk)=zparr
                  irr_mix (ji,jj,jk)=irr_inst(ji,jj,jk)
               
                  ! irradiance sum within the MLD
                  IF (sumz_hblt < hmld(ji,jj)) THEN
                     kblt=kblt+1
                     sumz_irrad_ML=sumz_irrad_ML+irr_mix(ji,jj,jk)*e3t(ji,jj,jk,Kmm)
                     sumz_hblt    =sumz_hblt+e3t(ji,jj,jk,Kmm)
                  ENDIF
               
               END DO
               
               ! irradiance mean average within the MLD
               irr_mix(ji,jj,1:kblt)=sumz_irrad_ML / MAX(1.0e-6,sumz_hblt)
         END DO
      END DO

      ! Initialize memory irradiance
      IF (kt==nittrc000) irr_mem(:,:,:) = irr_inst(:,:,:)
      !WRITE(numout,*) '   ==>>>   AGT: irradiance memory initialised ' 
      ! Phytoplankton photoadaptation timescale
      ! Forward time-stepping for memory irradiance
      DO jk=1, jpk
         DO jj=1, jpj
            DO ji=1, jpi
               irr_mem(ji,jj,jk)=irr_mem(ji,jj,jk)                                                    &
                                 !+ (irr_mix (ji,jj,jk)-irr_mem(ji,jj,jk))*MIN(1.d0,gam_irr_mem*rfact)*tmask(ji,jj,jk)
                                 + (irr_inst (ji,jj,jk)-irr_mem(ji,jj,jk))*MIN(1.d0,gam_irr_mem*rfact)*tmask(ji,jj,jk)
               irr_mem(ji,jj,jk) = MAX(irr_mem(ji,jj,jk),0.1d0)                  
            END DO
         END DO
      END DO

!      DEALLOCATE ( irr_inst )
!      DEALLOCATE ( irr_mix )

      IF( ln_timing == 1 )  CALL timing_stop ('trc_opt_bling_rb')

   END SUBROUTINE trc_opt_bling_rb

   SUBROUTINE trc_opt_bling_rgb (kt, Kmm, irr_inst, irr_mix)

      INTEGER , INTENT(in) :: kt
      INTEGER , INTENT(in) ::   Kmm   ! time level index
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: irr_mix, irr_inst

      INTEGER  :: ji, jj, jk, kblt, irgb
      REAL(wp) :: zchl, zc1, zc2, zc3, zc1km1, zc2km1, zc3km1
      REAL(wp) :: sumz_irrad_ML, sumz_hblt
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zekg, zekr, zekb

!      ALLOCATE ( irr_inst(jpi,jpj,jpk) )
!      ALLOCATE ( irr_mix(jpi,jpj,jpk) )
      ALLOCATE ( zekg(jpi,jpj,jpk) )
      ALLOCATE ( zekr(jpi,jpj,jpk) )
      ALLOCATE ( zekb(jpi,jpj,jpk) )

      IF( ln_timing == 1 )  CALL timing_start('trc_opt_bling_rgb')

      IF( kt == nittrc000 ) THEN
         irr_mix (:,:,:)=0._wp
         irr_inst(:,:,:)=0._wp
         irr_mem (:,:,:)=0._wp

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_opt_bling : BLING optic-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~ '
      ENDIF

      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      DO jk = 1, jpkm1                         !  --------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
               zchl = ( chl_bling(ji,jj,jk) + rtrn )
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )

               !                                                         
               zekb(ji,jj,jk) = xkrgb(1,irgb) * e3t(ji,jj,jk,Kmm)
               zekg(ji,jj,jk) = xkrgb(2,irgb) * e3t(ji,jj,jk,Kmm)
               zekr(ji,jj,jk) = xkrgb(3,irgb) * e3t(ji,jj,jk,Kmm)
            END DO
         END DO
      END DO

      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            !Irradiance at the surface 
            zc1 = parlux * qsr(ji,jj) * EXP( -0.5 * zekb(ji,jj,1) )
            zc2 = parlux * qsr(ji,jj) * EXP( -0.5 * zekg(ji,jj,1) )
            zc3 = parlux * qsr(ji,jj) * EXP( -0.5 * zekr(ji,jj,1) )

            irr_inst (ji,jj,1) = zc1 + zc2 + zc3 
            irr_mix  (ji,jj,1) = zc1 + zc2 + zc3 

            kblt=1
            sumz_irrad_ML=irr_mix(ji,jj,1)*e3t(ji,jj,1,Kmm)  
            sumz_hblt    =e3t(ji,jj,1,Kmm)

            DO jk = 2, nksrp   !<jk

               zc1km1=zc1
               zc2km1=zc2
               zc3km1=zc3

               zc1 = zc1km1 * EXP( -0.5 * ( zekb(ji,jj,jk-1) + zekb(ji,jj,jk) ) )
               zc2 = zc2km1 * EXP( -0.5 * ( zekg(ji,jj,jk-1) + zekg(ji,jj,jk) ) )
               zc3 = zc3km1 * EXP( -0.5 * ( zekr(ji,jj,jk-1) + zekr(ji,jj,jk) ) )

               irr_inst(ji,jj,jk) = zc1 + zc2 + zc3
               irr_mix (ji,jj,jk) = irr_inst(ji,jj,jk)

               ! irradiance sum within the MLD
               IF (sumz_hblt < hmld(ji,jj)) THEN
                  kblt=kblt+1
                  sumz_irrad_ML=sumz_irrad_ML+irr_mix(ji,jj,jk)*e3t(ji,jj,jk,Kmm)
                  sumz_hblt    =sumz_hblt+e3t(ji,jj,jk,Kmm)
               ENDIF
            ENDDO !jk>

            ! irradiance mean average within the MLD 
            irr_mix(ji,jj,1:kblt)=sumz_irrad_ML / MAX(1.0e-6,sumz_hblt)

         END DO
      END DO
      
      ! Initialize memory irradiance
      IF (kt==nittrc000) irr_mem(:,:,:)=irr_mix(:,:,:)

      ! Phytoplankton photoadaptation timescale
      ! Forward time-stepping for memory irradiance
      DO jk=1, jpk
         DO jj=1, jpj
            DO ji=1, jpi
               irr_mem(ji,jj,jk)=irr_mem(ji,jj,jk)                                                    &
                                 + (irr_mix (ji,jj,jk)-irr_mem(ji,jj,jk))*MIN(1.d0,gam_irr_mem*rfact)*tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      !DEALLOCATE ( irr_inst )
      !DEALLOCATE ( irr_mix )
      DEALLOCATE ( zekg )
      DEALLOCATE ( zekb )
      DEALLOCATE ( zekr )

      IF( ln_timing == 1 )  CALL timing_stop ('trc_opt_bling_rgb')

   END SUBROUTINE trc_opt_bling_rgb

   SUBROUTINE trc_opt_bling_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_opt_bling_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!----------------------------------------------------------------------
      !
      IF( ln_timing == 1 )  CALL timing_start('trc_opt_bling_init')
      !
      CALL trc_oce_rgb( xkrgb )                  ! tabulated attenuation coefficients
      nksrp = trc_oce_ext_lev( r_si2, 0.33e2 )   ! max level of light extinction (Blue Chl=0.01)
      !
      IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksrp, ' ref depth = ', gdepw_1d(nksrp+1), ' m'

      ! 
      IF( ln_timing == 1 )  CALL timing_stop('trc_opt_bling_init')
      !
   END SUBROUTINE trc_opt_bling_init

END MODULE trcopt_blingv0
