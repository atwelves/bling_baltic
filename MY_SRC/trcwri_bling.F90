MODULE trcwri_bling

   USE oce_trc
   USE trc
   USE iom
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_bling

CONTAINS

   SUBROUTINE trc_wri_bling( Kmm )
      INTEGER, INTENT(in)     :: Kmm  ! time level indices
      INTEGER           :: jn
      CHARACTER(len=20) :: cltra

      DO jn=jp_blg0, jp_blg1
         cltra = TRIM( ctrcnm(jn) )  
         CALL iom_put( cltra, tr(:,:,:,jn,Kmm))
      END DO

      ! diagnostic tracers output
      CALL iom_put( "CHL_bling" , chl_bling(:,:,:) * tmask(:,:,:) )
      CALL iom_put( "BIOMASS_P" , (biomass_p(:,:,:)+biomass_p_diaz(:,:,:)) * tmask(:,:,:) )
      !CALL iom_put( "BIOMASS_P_DIAZ" , biomass_p_diaz(:,:,:) * tmask(:,:,:) )
      CALL iom_put( "IRR_MEM"   , irr_mem(:,:,:) * tmask(:,:,:) )
      !CALL iom_put( "co3_ion"   ,   co3_ion(:,:,:) * tmask(:,:,:) )
      !CALL iom_put( "htotal"    ,    htotal(:,:,:) * tmask(:,:,:) )

   END SUBROUTINE trc_wri_bling

END MODULE trcwri_bling
