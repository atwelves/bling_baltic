MODULE trcsms_blingv0

   !!======================================================================
   !!                     ***  MODULE trcsms_bling  ***
   !! TOP :  BLINGv0 model main routine. Computes the sources and sinks
   !!======================================================================
   !! History :  MClaret@McGill@04/2016. Carbon cycle included
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc_oce
   USE iom

   ! BLINGv0 specific modules
   USE trcopt_blingv0
   USE trcext_blingv0
   USE vars_bling
   USE FMS_ocmip2_co2calc_mod

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_bling       ! called by trcsms.F90 module
   PUBLIC   trc_sms_init_bling

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_biomass_init ! structure of input fields (file informations, fields read)

#include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_sms_bling( kt, Kbb, Kmm, Krhs )
      !!-------------------------------------------------------------------------------
      !!                     ***  trc_sms_bling  ***
      !!
      !! ** History : default leap-frog scheme(MY_TRC) changed to FORWARD scheme
      !!-------------------------------------------------------------------------------
      !! The prefixes "f" refers to a "field" and "j" to a volumetric rate
      !! j is followed by the currency unit (p for phosphate, ca for calcium carbonate)
      !!------------------------------------------------------------------------------- 
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level index

      INTEGER  :: ji, jj, jk, jn

      REAL(wp) :: fpo4, fdop, ffed, foxy
      REAL(wp) :: ztra
      ! Irradiance k
      REAL(wp) :: po4_up, thetamax_fe, alpha_chl 
      ! Production
      REAL(wp) :: pc_tot, biomass_p_ts, theta, chl_dia, mulamb0expkT
      ! Phosphorous
      REAL(wp) :: frac_pop, jp_dop, fe2p_up, fpopkm1
      REAL(wp) :: zzz, wsink, oxy_up

      !!! --- AGT: Add local variables for nitrogen cycle --- !!!
      !IF ( ln_nitro ) THEN
      REAL(wp) :: fno3, fdon
      REAL(wp) :: no3_up
      REAL(wp) :: pc_tot_diaz, biomass_p_ts_diaz, mulamb0expkT_diaz
      REAL(wp) :: jn_don, fponkm1
      !ENDIF
      !!! ------!!!

      ! Iron
      REAL(wp) :: jfe_pofe, fpofekm1
      REAL(wp) :: dum5, dum2, dum3

!DO_CARBON
      ! Carbonate system
      REAL(wp) :: s_over_p, co3_solubility, fcaco3km1
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jca_uptake, jca_reminp, fcaco3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jdic, jalk
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zremin_caco3, frac_lg
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: pco2_surf
!DO_CARBON

      ! global sum
      REAL(wp) ::  sumrecycle, sumremin, sumuptake, sumuptake2, &
                   sumremin2, sumpo4, sumdop, sumdic

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: rho
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: expkT, irrk, pc_m, mu
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: irr_inst, irr_mix
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: def_fe, feprime, kfe_eq_lig, fpofe
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jp_pop, fpop, zremin
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jp_uptake, jp_remin, jp_recycle
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jfe_uptake, jfe_remin, jfe_recycle
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jfe_ads_inorg, jfe_ads_org
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jpo4, jdop, jfed, joxy
      !!! --- AGT: Declare arrays for nitrogen cycle --- !!!
      !IF( ln_nitro ) THEN
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jno3, jdon
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jn_pon, fpon
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: jn_uptake, jn_fix, jn_remin, jn_recycle
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: pc_m_diaz, mu_diaz
      !ENDIF
      !!! ------ !!!
      REAL(wp), ALLOCATABLE, DIMENSION(:)     :: dum4
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: xnegtr

      !!----------------------------------------------------------------------
      !
      IF( ln_timing == 1 )  CALL timing_start('trc_sms_bling')
      !
      IF (kt==nittrc000) THEN ! Aug 28, 2015, xhu
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) ' trc_sms_bling:  BLINGv0 model'
          IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      ENDIF

      ! allocate matrix variables
      ALLOCATE ( rho(jpi,jpj,jpk) )
      ALLOCATE ( expkT(jpi,jpj,jpk) )
      ALLOCATE ( irr_inst(jpi,jpj,jpk) )
      ALLOCATE ( irr_mix(jpi,jpj,jpk) )
      ALLOCATE ( irrk(jpi,jpj,jpk) )
      ALLOCATE ( pc_m(jpi,jpj,jpk) )
      ALLOCATE ( mu(jpi,jpj,jpk) )
      ALLOCATE ( def_fe(jpi,jpj,jpk) )
      ALLOCATE ( feprime(jpi,jpj,jpk) )
      ALLOCATE ( kfe_eq_lig(jpi,jpj,jpk) )
      ALLOCATE ( fpofe(jpi,jpj,jpk) )
      ALLOCATE ( jp_pop(jpi,jpj,jpk) )
      ALLOCATE ( fpop(jpi,jpj,jpk) )
      ALLOCATE ( zremin(jpi,jpj,jpk) )
      ALLOCATE ( jp_uptake(jpi,jpj,jpk) )
      ALLOCATE ( jp_remin(jpi,jpj,jpk) )
      ALLOCATE ( jp_recycle(jpi,jpj,jpk) )
      ALLOCATE ( jfe_uptake(jpi,jpj,jpk) )
      ALLOCATE ( jfe_remin(jpi,jpj,jpk) )
      ALLOCATE ( jfe_recycle(jpi,jpj,jpk) )
      ALLOCATE ( jfe_ads_inorg(jpi,jpj,jpk) )
      ALLOCATE ( jfe_ads_org(jpi,jpj,jpk) )
      ALLOCATE ( jpo4(jpi,jpj,jpk) )
      ALLOCATE ( jdop(jpi,jpj,jpk) )
      ALLOCATE ( jfed(jpi,jpj,jpk) )
      ALLOCATE ( joxy(jpi,jpj,jpk) )
      !!! --- AGT: Allocate memory for nitrogen cycle --- !!! 
 !     IF ( ln_nitro ) THEN
              ALLOCATE ( jno3(jpi,jpj,jpk) )
              ALLOCATE ( jdon(jpi,jpj,jpk) )
              ALLOCATE ( jn_pon(jpi,jpj,jpk) )
              ALLOCATE ( fpon(jpi,jpj,jpk) )
              ALLOCATE ( jn_uptake(jpi,jpj,jpk) )
              ALLOCATE ( jn_fix(jpi,jpj,jpk) )
              ALLOCATE ( jn_remin(jpi,jpj,jpk) )
              ALLOCATE ( jn_recycle(jpi,jpj,jpk) )
              ALLOCATE ( pc_m_diaz(jpi,jpj,jpk) )
              ALLOCATE ( mu_diaz(jpi,jpj,jpk) )
      !!! ------ !!!
  !    ENDIF
      !DO_CARBON
      ALLOCATE ( jca_uptake(jpi,jpj,jpk) )
      ALLOCATE ( jca_reminp(jpi,jpj,jpk) )
      ALLOCATE ( fcaco3(jpi,jpj,jpk) )
      ALLOCATE ( jdic(jpi,jpj,jpk) )
      ALLOCATE ( jalk(jpi,jpj,jpk) )
      ALLOCATE ( zremin_caco3(jpi,jpj,jpk) )
      ALLOCATE ( frac_lg(jpi,jpj,jpk) )
      ALLOCATE ( pco2_surf(jpi,jpj) )
      ALLOCATE ( dum4(jpk) )
      ALLOCATE ( xnegtr(jpi,jpj,jpk) )

      IF ( ln_bling_mass) THEN      !   Write values for phosphate budget
        CALL trc_sms_bling_mass_conserv (kt, Kbb, Kmm, Krhs)
      ENDIF
  !    WRITE(numout,*) '   ==>>>   AGT: trc_sms_bling '
      rho(:,:,:)=rhd(:,:,:)*1035.e0+1035.e0

      jk=1
      CALL FMS_ocmip2_co2calc( 1,jpi,1,jpj,1,jpi,1,jpj,tmask(:,:,jk)-tmask(:,:,jk)+1.e0                &
                              ,ts(:,:,jk,jp_tem,Kmm),ts(:,:,jk,jp_sal,Kmm)                                   &
                              ,tr(:,:,jk,jpDIC_bling,Kmm)/rho(:,:,jk),tr(:,:,jk,jpPO4_bling,Kmm)/rho(:,:,jk) &
                              ,tr(:,:,jk,jpPO4_bling,Kmm)/rho(:,:,jk),tr(:,:,jk,jpalk_bling,Kmm)/rho(:,:,jk) &
                              ,htotal(:,:,jk)*htotal_scale_lo                                          &
                              ,htotal(:,:,jk)*htotal_scale_hi                                          &
                              ,htotal(:,:,jk)                                                          & 
                              ,co3_ion=co3_ion(:,:,jk)                                                 &
                              ,co2star=co2_csurf,alpha=co2_alpha                                       &
                              ,pCO2surf=pco2_surf                                                       )

      ! convert from mol/kg to mol/m3
      co2_alpha(:,:)  =co2_alpha(:,:)  *rho(:,:,1)
      co2_csurf(:,:)  =co2_csurf(:,:)  *rho(:,:,1)

      DO jk=2, jpk

          CALL FMS_ocmip2_co2calc( 1,jpi,1,jpj,1,jpi,1,jpj,tmask(:,:,jk)-tmask(:,:,jk)+1.e0                &
                                  ,ts(:,:,jk,jp_tem,Kmm),ts(:,:,jk,jp_sal,Kmm)                                   &
                                  ,tr(:,:,jk,jpDIC_bling,Kmm)/rho(:,:,jk),tr(:,:,jk,jpPO4_bling,Kmm)/rho(:,:,jk) &
                                  ,tr(:,:,jk,jpPO4_bling,Kmm)/rho(:,:,jk),tr(:,:,jk,jpalk_bling,Kmm)/rho(:,:,jk) &
                                  ,htotal(:,:,jk)*htotal_scale_lo                              &
                                  ,htotal(:,:,jk)*htotal_scale_hi                              &
                                  ,htotal(:,:,jk),co3_ion=co3_ion(:,:,jk)                       )

          ! MC note 04/02/16: I do not use the mask. NEMO advects fields w/o mask,
          ! and I'm not sure how masking co3_ion and htotal could affect DIC/ALK fields.

      ENDDO

      ! NOTE: htotal and co3_ion are left to mol/kg deliberately. 
      ! See zremin_caco3 computation for comment on co3_ion units
      !----------------------------------------------------------------------
      ! BLING model

      CALL trc_opt_bling_rb  (kt, Kmm, irr_inst, irr_mix)  ! optical model (red and blue wavelengths)
      !CALL trc_opt_bling_rgb (kt, irr_inst, irr_mix) ! optical model (red, blue, and green wavelengths)

      DO ji=1, jpi
          DO jj=1, jpj

              dum4(:) = 0.d0

              DO jk=1, jpk

               ! ----------------------------------------------------------
               ! negative trophic variables DO not contribute to the fluxes
               ! ----------------------------------------------------------

               ! [mol/m3]
               fpo4 = MAX( 0.e0, tr(ji,jj,jk,jpPO4_bling,Kmm) )
               fdop = MAX( 0.e0, tr(ji,jj,jk,jpDOP_bling,Kmm) )
               ffed = MAX( 0.e0, tr(ji,jj,jk,jpFed_bling,Kmm) )
               foxy = tr(ji,jj,jk,jpOxy_bling,Kmm)
               !!! --- AGT: Add fluxes for nitrogen cycle --- !!!
              ! IF ( ln_nitro ) THEN
                       fno3 = MAX( 0.e0, tr(ji,jj,jk,jpNO3_bling,Kmm) )
                       fdon = MAX( 0.e0, tr(ji,jj,jk,jpDON_bling,Kmm) )
               !!! ------ !!!
              ! ENDIF

               ! ----------------------------------------------------------
               ! TEMPERATURE DEPENDENCE
               ! NB: The temperature effect of Eppley (1972) is used instead 
               !     of that in Geider et al (1997) for both simplicity and 
               !     to incorporate combined effects on uptake, incorporation
               !     into organic matter and photorespiration.  Values of PCmax
               !     are normalized to 0C rather than 20C in Geider et al.(1997)
               ! ----------------------------------------------------------

               ! [no units]
               expkT(ji,jj,jk)=EXP(kappa_eppley*ts(ji,jj,jk,jp_tem,Kmm))

               ! ----------------------------------------------------------
               ! Phytoplankton are assumed to grow according to the general properties 
               ! described in Geider (1997). This formulation gives a biomass-specific 
               ! growthrate as a function of light, nutrient limitation, and 
               ! temperature. We modify this relationship slightly here, as described 
               ! below, and also use the assumption of steady state growth vs. loss to 
               ! derive a simple relationship between growth rate, biomass and uptake.
               ! ----------------------------------------------------------
               ! First, we calculate the limitation terms for PO4 and Fe, and the 
               ! Fe-limited Chl:C maximum.
               ! The light-saturated maximal photosynthesis rate term (pc_m) is simply 
               ! the product of a prescribed maximal photosynthesis rate (pc_0), the 
               ! Eppley temperature dependence, and a Liebig limitation (the minimum
               ! of Michaelis-Menton PO4-limitation, or iron-limitation).
               ! The iron limitation term has a lower limit of def_fe_min 
               ! and is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
               ! so that it approaches 1 as fed approaches infinity. Thus, 
               ! it's of comparable magnitude to the PO4 limitation term.
               !
               ! Fe limitation acts in two additional mechanisms:
               ! 1. By reducing the maximum achievable Chl:C ratio 
               ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
               ! approach a prescribed minimum Chl:C (thetamin) under extreme
               ! Fe-limitation.
               ! 2. By reducing alpha (the initial slope of the P-I curve) under Fe-
               ! limitation.
               ! ----------------------------------------------------------

               ! Iron uptake [no units]
               fe2p_up = fe2p_max * ffed / (kfe + ffed)  ![mol Fe/mol P]
               def_fe(ji,jj,jk)  = def_fe_min + (1.d0-def_fe_min) &
                                  *fe2p_up/(kfe2p_up+fe2p_up)*(kfe2p_up+fe2p_max)/fe2p_max 

               ! Phosphate uptake [no units]
               po4_up = fpo4 /( kpo4 + fpo4 )

               !!! --- AGT: Add nitrogen limitation --- !!!
              ! IF ( ln_nitro ) THEN
                       no3_up = fno3 /( kno3 + fno3 )
                       ! Maximum production (units of pc_0 [s-1])
                       pc_m(ji,jj,jk) = pc_0 * expkT(ji,jj,jk) * MIN(po4_up,no3_up,def_fe(ji,jj,jk))
                       ! add diazotrophs:
                       pc_m_diaz(ji,jj,jk) = pc_0 * expkT(ji,jj,jk) * MIN(po4_up,def_fe(ji,jj,jk))
              ! ELSE 
               !!! ------ !!!
               ! Maximum production (units of pc_0 [s-1])
              ! pc_m(ji,jj,jk) = pc_0 * expkT(ji,jj,jk) * MIN(po4_up,def_fe(ji,jj,jk)) 
              ! ENDIF

               ! Iron limitation on photosyntesis machinery
               thetamax_fe=thetamax_lo + (thetamax_hi - thetamax_lo)*def_fe(ji,jj,jk)
               alpha_chl  =alpha_min   + (alpha_max   - alpha_min  )*def_fe(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
               ! calculated. This requires a prescribed quantum yield, alpha.
               ! The iron deficiency term is included here as a multiplier of the 
               ! thetamax_fe to represent the importance of Fe in forming chlorophyll
               ! accessory antennae, which do not affect the Chl:C but still affect the
               ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
               ! 2004).
               !-----------------------------------------------------------------------
               !WRITE(numout,*) '   ==>>>   AGT: iron limitation complete '
               ! I_k [W/m2]
               irrk(ji,jj,jk) = pc_m(ji,jj,jk) / (epsln + alpha_chl*thetamax_fe) + irr_mem(ji,jj,jk)*0.5d0
               !WRITE(numout,*) '   ==>>>   AGT: light saturation calculated '
               !-----------------------------------------------------------------------
               ! Now we can calculate the carbon-specific photosynthesis rate, pc_tot.
               !-----------------------------------------------------------------------

               ![s-1]
               !pc_tot = pc_m(ji,jj,jk)*(1.d0-EXP(-irr_mix(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))
               pc_tot = pc_m(ji,jj,jk)*(1.d0-EXP(-irr_inst(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))
               !!! --- AGT --- !!!
              ! IF ( ln_nitro ) THEN
               pc_tot_diaz = pc_m_diaz(ji,jj,jk)*(1.d0-EXP(-irr_inst(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))
              ! ENDIF
               !!! ------ !!!
               !-----------------------------------------------------------------------
               ! Next, we account for the maintenance effort that phytoplankton must 
               ! exert in order to combat decay. This is prescribed as a fraction of the
               ! light-saturated photosynthesis rate, resp_frac. The result of this is 
               ! to set a level of energy availability below which net growth (and 
               ! therefore nutrient uptake) is zero, given by resp_frac * pc_m.
               !-----------------------------------------------------------------------

               ! Net total production [s-1]
               mu(ji,jj,jk) = MAX (0.d0,pc_tot-resp_frac*pc_m(ji,jj,jk))
               !!! --- AGT --- !!!
              ! IF ( ln_nitro ) THEN
               mu_diaz(ji,jj,jk) = MAX (0.d0,pc_tot_diaz-resp_frac*pc_m_diaz(ji,jj,jk))
              ! ENDIF
               !!! ------ !!!
               !-----------------------------------------------------------------------
               ! We now must convert this net carbon-specific growth rate to nutrient 
               ! uptake rates, the quantities we are interested in. Since we have no 
               ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
               ! to calculate an implicit biomass from the uptake rate through the  
               ! application of a simple idealized grazing law. This has the effect of 
               ! reducing uptake in low growth-rate regimes and increasing uptake in 
               ! high growth-rate regimes - essentially a non-linear amplification of 
               ! the growth rate variability. The result is:
               !-----------------------------------------------------------------------

               ! Biomass (units of pstar [mol P/m3])
               mulamb0expkT = mu(ji,jj,jk)/(lambda0*expkT(ji,jj,jk))  ![no units]
               biomass_p_ts = p_star*mulamb0expkT*(1.d0+(mulamb0expkT)**2)
               
               ! correct units!
               IF (kt==nittrc000) biomass_p(ji,jj,jk)=1.d3 * biomass_p(ji,jj,jk)

               biomass_p(ji,jj,jk) =   biomass_p(ji,jj,jk) &
                                    + (biomass_p_ts-biomass_p(ji,jj,jk))*MIN(1.d0,gam_biomass*rfact)!*tmask(ji,jj,jk)

               !!! --- AGT: this needs some more thought... --- !!!
              ! IF ( ln_nitro ) THEN
                       mulamb0expkT_diaz = mu_diaz(ji,jj,jk)/(lambda0*expkT(ji,jj,jk))  ![no units]
                       biomass_p_ts_diaz = p_star*mulamb0expkT_diaz*(1.d0+(mulamb0expkT_diaz)**2)
                       IF (kt==nittrc000) biomass_p_diaz(ji,jj,jk)=epsln

                       biomass_p_diaz(ji,jj,jk) =   biomass_p_diaz(ji,jj,jk) &
                                    + (biomass_p_ts_diaz-biomass_p_diaz(ji,jj,jk))*MIN(1.d0,gam_biomass*rfact)!*tmask(ji,jj,jk)
              ! ENDIF
               !!! ------ !!!

               !if (ji==80 .and. jj==60 .and. jk==1) write(*,'(I3,5(1X,E11.4))') kt, &
               !pc_0, expkT(ji,jj,jk), po4_up, def_fe(ji,jj,jk), pc_tot,mu(ji,jj,jk),biomass_p(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! We can now use the diagnostic biomass to calculate the chlorophyll
               ! concentration:
               !-----------------------------------------------------------------------

               ! Chl:C ration [g chl/g C]
               theta   = thetamax_fe / (1.d0 + (thetamax_fe*alpha_chl*irr_mem(ji,jj,jk))/(2.d0*pc_m(ji,jj,jk)+epsln) )

               ! Chl biomass [mg chl/m3]
               !!! --- AGT: include diazotrophs in chlorophyll calculation --- !!!
              ! IF ( ln_nitro ) THEN
                       chl_dia = (biomass_p(ji,jj,jk)+biomass_p_diaz(ji,jj,jk)) * c2p * 12.011e+3 * theta
              ! ELSE
               !!! ------ !!!
              ! chl_dia = biomass_p(ji,jj,jk) * c2p * 12.011e+6 * theta !* tmask(ji,jj,jk) 
              ! ENDIF
               chl_bling(ji,jj,jk) = MAX(chl_min, chl_dia)

               !--------------------------------------------------
               ! PHOSPHORUS CYCLE
               !--------------------------------------------------
               ! The uptake of nutrients is assumed to contribute to the growth of
               ! phytoplankton, which subsequently die and are consumed by heterotrophs.
               ! This can involve the transfer of nutrient elements between many
               ! organic pools, both particulate and dissolved, with complex histories.
               ! We take a simple approach here, partitioning the total uptake into two
               ! fractions - sinking and non-sinking - as a function of temperature, 
               ! following Dunne et al. (2005). 
               ! Then, the non-sinking fraction is further subdivided, such that the 
               ! majority is recycled instantaneously to the inorganic nutrient pool,
               ! representing the fast turnover of labile dissolved organic matter via
               ! the microbial loop, and the remainder is converted to semi-labile
               ! dissolved organic matter. Iron and phosphorus are treated identically 
               ! for the first step, but all iron is recycled instantaneously in the
               ! second step (i.e. there is no dissolved organic iron pool).
               !-----------------------------------------------------------------------

               !!! --- AGT: Add diazotroph contribution... --- !!!
           !    IF ( ln_nitro ) THEN
                       ! Phosphorous uptake flux [mol P/m3/s]
                       jp_uptake(ji,jj,jk)=biomass_p(ji,jj,jk)*mu(ji,jj,jk) + biomass_p_diaz(ji,jj,jk)*mu_diaz(ji,jj,jk)
                       !!!!!!!!! how to calculate export fraction with diazotrophs??????
            !   ELSE
               !!! ------ !!!
               ! Phosphorous uptake flux [mol P/m3/s]
             !  jp_uptake(ji,jj,jk)=biomass_p(ji,jj,jk)*mu(ji,jj,jk)
              ! ENDIF
               ! [no units]
               frac_pop=(phi_sm+phi_lg*(mulamb0expkT)**2) / (1+(mulamb0expkT)**2) * EXP(kappa_remin*ts(ji,jj,jk,jp_tem,Kmm))

               ! [mol P/m3/s]
               jp_pop(ji,jj,jk)=frac_pop*jp_uptake(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Then the remainder is divided between instantaneously recycled and
               ! long-lived dissolved organic matter,
               !-----------------------------------------------------------------------
               !
               ! [mol P/m3/s]
               jp_dop=phi_dop*(jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk))

               ! [mol P/m3/s]
               jp_recycle(ji,jj,jk)=jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk)-jp_dop
               !!! --- AGT: Start by adding nitrogen uptake and recycling as redfield multiple of phosphorous tracer, minus fixation: --- !!!
         !      IF ( ln_nitro ) THEN
                       ! calculate nitrogen fixation
                       jn_fix(ji,jj,jk) = biomass_p_diaz(ji,jj,jk)*mu_diaz(ji,jj,jk)
                       ! calculate nitrate uptake
                       jn_uptake(ji,jj,jk)=n2p*jp_uptake(ji,jj,jk) - jn_fix(ji,jj,jk)
                       ! calculate particulate flux of nitrogen
                       jn_pon(ji,jj,jk)=frac_pop*(jn_uptake(ji,jj,jk)+jn_fix(ji,jj,jk))
                       ! calculate dissolved flux of nitrogen
                       jn_don=phi_dop*(jn_uptake(ji,jj,jk)+jn_fix(ji,jj,jk)-jn_pon(ji,jj,jk))
                       ! calculate nitrate recycling
                       jn_recycle(ji,jj,jk)=jn_uptake(ji,jj,jk)+jn_fix(ji,jj,jk)-jn_pon(ji,jj,jk)-jn_don
          !     ENDIF
               !!! ------ !!!
!DO_CARBON<
               !---------------------------------------------------------------------
               ! As a helpful diagnostic, the implied fraction of production by large 
               ! phytoplankton is calculated, also following Dunne et al. 2005. This
               ! could be done more simply, but is done here in a complicated
               ! way as a sanity check. Looks fine.
               ! Note the calculation is made in P units, rather than C.
               ! This is also used for the CaCO3 production.

               ! MC. Note: s_over_p = (mu/lambda)**2 for gamma_b=0. 
               ! To do the maths consider jp_uptake=mu*B and lambda=lambda0*expkT  
               ! [no units]

               !!! --- AGT: For compatibility with nitro version calculate biomass x growth rate explicitly for non-diazotrophs
               s_over_p = ( ( 1.d0 + 4.d0 * biomass_p(ji,jj,jk)*mu(ji,jj,jk) / (lambda0*expkT(ji,jj,jk)*p_star) )**0.5 - 1.d0 ) / 2.d0
               !!! ------ !!!
               ![no units]
               frac_lg(ji,jj,jk) = s_over_p / (1.d0 + s_over_p)

               !-----------------------------------------------------------------------
               ! Calcium carbonate production
               ! Alkalinity is consumed through the production of CaCO3. Here, this is
               ! simply a linear function of the implied growth rate of small
               ! phytoplankton, which gave a reasonably good fit to the global 
               ! observational synthesis of Dunne (in prep., 2009).
               ! This is consistent with the findings of Jin et al. (GBC,2006).
               !-----------------------------------------------------------------------

               ! [mol Ca/m3/s]
               !!! --- AGT: For compatibility with nitro version calculate biomass x growth rate explicitly for non-diazotrophs
               jca_uptake(ji,jj,jk) = (1.d0-frac_lg(ji,jj,jk))*biomass_p(ji,jj,jk)*mu(ji,jj,jk)*ca2p         
               !!! ------ !!!
!>DO_CARBON
               !---------------------------------------------------------------------
               ! IRON
               !---------------------------------------------------------------------
               ! Iron is then taken up as a function of PO4 uptake and iron limitation,
               ! with a maximum Fe:P uptake ratio of fe2p_max:
               !-----------------------------------------------------------------------

               ! Iron uptake/POFe/recycle fluxes (units of jp_uptake [mol Fe/m3/s])
               jfe_uptake(ji,jj,jk) =jp_uptake(ji,jj,jk)*fe2p_up
               jfe_pofe             =frac_pop*jfe_uptake(ji,jj,jk)
               jfe_recycle(ji,jj,jk)=jfe_uptake(ji,jj,jk)-jfe_pofe

               !-----------------------------------------------------------------------
               ! Calculate free and inorganically associated iron concentration for
               ! scavenging.
               ! We assume that there is a 
               ! spectrum of iron ligands present in seawater, with varying binding
               ! strengths and whose composition varies with light and iron 
               ! concentrations. For example, photodissocation of ligand complexes 
               ! occurs under bright light, weakening the binding strength 
               ! (e.g. Barbeau et al., Nature 2001), while at very low iron 
               ! concentrations (order kfe_eq_lig_femin), siderophores are thought
               ! to be produced as a response to extreme
               ! iron stress.
               ! In anoxic waters, iron should be reduced, and therefore mostly 
               ! immune to scavenging. Easiest way to do this is to skip the feprime
               ! calculation if oxygen is less than 0.
               !-----------------------------------------------------------------------

               ! Calculate Fe prime [mol Fe/m3]
               if (foxy > oxy_min ) then
                 dum5       = irr_inst(ji,jj,jk)**2/(kfe_eq_lig_irr**2+irr_inst(ji,jj,jk)**2)        ! [no units]
                 dum2       = max(  0.e0, min( 1.e0, 1.2e0*(ffed-kfe_eq_lig_femin)/(epsln+ffed) )  ) ! [no units]
                 kfe_eq_lig(ji,jj,jk) = kfe_eq_lig_max -(kfe_eq_lig_max-kfe_eq_lig_min)*dum5*dum2    ! [m3/mol Fe]
                 feprime(ji,jj,jk)= 1.e0 + kfe_eq_lig(ji,jj,jk)*(felig_bkg - ffed)                   ! [no units]
                 ! units of (kfe_eq_lig)-1[mol Fe/m3]
                 feprime(ji,jj,jk)= (-feprime(ji,jj,jk) + sqrt(feprime(ji,jj,jk)**2 + 4.e0*kfe_eq_lig(ji,jj,jk)*ffed) )/(2.e0*kfe_eq_lig(ji,jj,jk))
               else
                 feprime(ji,jj,jk) = 0.e0
               endif

               ! [mol Fe/m3/s]
               jfe_ads_inorg(ji,jj,jk) = min( 0.5d0, kfe_inorg*sqrt(feprime(ji,jj,jk)) )*feprime(ji,jj,jk)

               ! [mol Fe/m3/s]
               dum4(jk) = jfe_pofe + jfe_ads_inorg(ji,jj,jk)

               !--------------------------------------------------
               ! COMPUTE TRENDS w/o remineralization processes
               !--------------------------------------------------
               !
               ! [mol P/m3/s]
               jpo4(ji,jj,jk) =   jp_recycle(ji,jj,jk) + gamma_dop*fdop -jp_uptake(ji,jj,jk)
               jdop(ji,jj,jk) = - gamma_dop*fdop + phi_dop*(jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk))
               ! [mol Fe/m3/s]
               jfed(ji,jj,jk) =   jfe_recycle(ji,jj,jk)-jfe_uptake(ji,jj,jk)-jfe_ads_inorg(ji,jj,jk)
               !!! --- AGT: Incoroporate nitrogen cycle --- !!!
               jno3(ji,jj,jk) =   jn_recycle(ji,jj,jk) + gamma_dop*fdon -jn_uptake(ji,jj,jk)
               jdon(ji,jj,jk) = - gamma_dop*fdon + phi_dop*(jn_uptake(ji,jj,jk)-jn_pon(ji,jj,jk))
               !!! --- !!!
               !DO_CARBON
               ! [mol C/m3/s]
               jdic(ji,jj,jk) = - jca_uptake(ji,jj,jk) 
               jalk(ji,jj,jk) = - 2.d0*jca_uptake(ji,jj,jk)

            ENDDO


            !-----------------------------------------------------------------------
            ! SINKING AND REMINERALIZATION
            !-----------------------------------------------------------------------
            ! Calculate the remineralization lengthscale matrix, zremin, a function 
            ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
            ! then  increases linearly with depth.
            ! The remineralization rate is a function of oxygen concentrations,
            ! following a Holling type 2 dependence, decreasing to a minimum value
            ! of remin_min. This is ad hoc, following work by Bianchi, Sarmiento,
            ! Galbraith and Kwon (unpublished).
            !-----------------------------------------------------------------------
            ! In general, the flux at the bottom of a grid cell should equal
            ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
            ! where Ft is the flux at the top, and prod*dz is the integrated 
            ! production of new sinking particles within the layer.
            ! Since Ft=0 in the first layer,
            !---------------------------------------------------------------------
            ! Calculate co3_ion first. Used to compute caco3 remineralization lengthscale
            ! Also calculate co2 fluxes csurf and alpha for the next round of exchange
            ! Note a scaled value of the PO4, rather than SiOH3, is used for all 
            ! calculations since there is no prognostic silica cycle 
            ! GFDL subroutine used to compute H+ that includes the OCMIP2 protocol
            !---------------------------------------------------------------------
            
            !!! --- AGT --- !!!
            ! Here try to account for decaying biomass. If growth rate is <0 then 
            ! change in biomass = - adjustment time scale x biomass
            ! assume that this leads to a sinking flux of particulate organic matter
            !
            ! Hence after production ends, the biomass field effectively acts as 
            ! a detritus field.  Note that with this modification there is now 
            ! a non-conservative term in the N/P budgets!
            
            IF ( mulamb0expkT<0 .AND. mulamb0expkT_diaz<0 ) THEN
                    jp_pop(ji,jj,jk) = - (biomass_p(ji,jj,jk) + biomass_p_diaz(ji,jj,jk))*MIN(1.d0,gam_biomass*rfact)
                    jn_pon(ji,jj,jk) = n2p * jp_pop(ji,jj,jk)
            ENDIF
                        
            !!! ------ !!!

            ! k=1: surface layer
            jk=1
            ! [m]
            zzz  =e3t(ji,jj,jk,Kmm)

            ! Sinking rate [m/s]
            IF (zzz .lt. wsink0_z) THEN
               wsink=wsink0
            ELSE
               wsink=wsink0+wsink_acc*(zzz-wsink0_z)
            ENDIF

            ! Remineralization lengthscale (oxygen dependent process)
            ! [mol O2/m3/s]
            foxy = tr(ji,jj,jk,jpOxy_bling,Kmm)       
            IF (foxy>oxy_min) THEN
            ! [no units]
                oxy_up =foxy**2 / (koxy**2 + foxy**2) 
            ELSE
                oxy_up = 0.0d0
            ENDIF
            ! [m-1]
            zremin(ji,jj,jk) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)
            

            ! [mol P/m2/s]
            fpop(ji,jj,jk)    = jp_pop(ji,jj,jk)*e3t(ji,jj,jk,Kmm)/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk)) 
            
            ! [mol P/m3/s]
            jp_remin(ji,jj,jk)=(jp_pop(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fpop(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))

            !!! --- AGT: Add nitrogen remineralization --- !!!
       !     IF ( ln_nitro ) THEN
                    ! [mol N/m2/s]
                    fpon(ji,jj,jk)    = jn_pon(ji,jj,jk)*e3t(ji,jj,jk,Kmm)/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk))
                    ! [mol N/m3/s]
                    jn_remin(ji,jj,jk)=(jn_pon(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fpon(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))
        !    ENDIF
            !!! ------ !!!

            !-----------------------------------------------------------------------
!DO_CARBON<-
            ! Using Sayles for calcite solubility (will change to Mucci later) [mol/kg]
            ! NOTE: co3_solubility units cancel with co3_ion units also in mol/kg.
            ! This algorithm is ported exactly from blingv0 in GFDL, which units are in mol/kg.
            co3_solubility=max(  4.95e-7    * exp( 0.05021/(ts(ji,jj,jk,jp_tem,Kmm)+273.15)*zzz )  &
                                *3.42031e-3 * rho0_co3sol * rho0_co3sol / max(epsln, ts(ji,jj,jk,jp_sal,Kmm))    &
                               , epsln )

            ! CaCO3 dissolution lengthscale is a function of the saturation
            ! state, CO3 / CO3solubility, such that the lengthscale decreases
            ! from infinity at CO3 > CO3solubility to zero at CO3 = 0.

            ! [m-1]
            zremin_caco3(ji,jj,jk)=  1.d0 / ca_remin_depth  &
                                  *(1.d0 - min(1.d0, co3_ion(ji,jj,jk)/(co3_solubility+epsln) ))

            ! Generate CaCO3 sinking flux, and dissolve it through the water column. 
            ! Same as for other elements - see below for more detailed explanation.

            ! MC note: fcaco3 is computed as fpop, and jca_remin as jp_remin
            ! [mol C/m2/s]
            fcaco3    (ji,jj,1)= jca_uptake(ji,jj,jk)*e3t(ji,jj,jk,Kmm)/(1.d0+e3t(ji,jj,jk,Kmm)*zremin_caco3(ji,jj,jk))
            jca_reminp(ji,jj,1)=(jca_uptake(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fcaco3(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))

            ! Add remineralization terms to sources
            ! [mol C/m2/s]
            jdic(ji,jj,jk)=jdic(ji,jj,jk)+     jca_reminp(ji,jj,jk)
            ! [mol eq/m2/s]
            jalk(ji,jj,jk)=jalk(ji,jj,jk)+2.d0*jca_reminp(ji,jj,jk)
!DO_CARBON>

            !-----------------------------------------------------------------------
            ! Now, calculate the Fe adsorption using this fpop:
            ! The absolute first order rate constant is calculated from the 
            ! concentration of organic particles, after Parekh et al. (2005). Never
            !  allowed to be greater than 1/2dt for numerical stability.
            !-----------------------------------------------------------------------

            !Iron [gC/m3]
            dum3                 = (fpop(ji,jj,jk)*c2p*12.011d0/(epsln+wsink))**(0.58)
            ! [mol Fe/m3/s]
            jfe_ads_org(ji,jj,jk) = min (0.5d0, kfe_org*dum3)*feprime(ji,jj,jk)
            ! [mol Fe/m2/s]
            dum4(jk)              =( dum4(jk)+jfe_ads_org(ji,jj,jk) )*e3t(ji,jj,jk,Kmm)
            fpofe(ji,jj,jk)       =  dum4(jk)/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk))
            ! [mol Fe/m3/s]
            jfe_remin(ji,jj,jk)   =( dum4(jk)-fpofe(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm) )

            ! Add remineralization terms to trends
            ! [mol P/m3/s]
            jpo4(ji,jj,1)=jpo4(ji,jj,1)+(1.d0-phi_dop)*jp_remin(ji,jj,1)
            jdop(ji,jj,1)=jdop(ji,jj,1)+       phi_dop*jp_remin(ji,jj,1)
            ! [mol Fe/m3/s]
            jfed(ji,jj,1)=jfed(ji,jj,1)+jfe_remin(ji,jj,1)-jfe_ads_org(ji,jj,1)   

            ! k=2:NK: rest of the water column
            DO jk=2, jpk

               fpopkm1 = fpop(ji,jj,jk-1)
               fcaco3km1 = fcaco3(ji,jj,jk-1)
               fpofekm1=fpofe(ji,jj,jk-1)
    
               zzz=zzz+e3t(ji,jj,jk,Kmm)

               IF (zzz .lt. wsink0_z) THEN
                  wsink=wsink0
               ELSE
                  wsink=wsink0+wsink_acc*(zzz-wsink0_z)
               ENDIF

               ! Phosphorous
               foxy = tr(ji,jj,jk,jpOxy_bling,Kmm)
               
               IF (foxy>oxy_min) THEN
                  oxy_up=foxy**2 / (koxy**2 + foxy**2)
               ELSE
                  oxy_up=0.0d0
               ENDIF
               ! [m-1]
               zremin(ji,jj,jk) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)

               ! [mol P/m2/s]
               fpop(ji,jj,jk)    =(fpopkm1+jp_pop(ji,jj,jk)*e3t(ji,jj,jk,Kmm))/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk)) 
               ! [mol P/kg/s]
               jp_remin(ji,jj,jk)=(fpopkm1+jp_pop(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fpop(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))

               !!! --- AGT: Add nitrogen remineralization --- !!!
     !          IF ( ln_nitro ) THEN    
                       fponkm1 = fpon(ji,jj,jk-1)
                       ! [mol N/m2/s]
                       fpon(ji,jj,jk)    = (fponkm1+jn_pon(ji,jj,jk)*e3t(ji,jj,jk,Kmm))/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk))
                       ! [mol N/m3/s]
                       jn_remin(ji,jj,jk)= (fponkm1+jn_pon(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fpon(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))
      !         ENDIF
               !!! ------ !!!


!DO_CARBON<
               co3_solubility=max(  4.95e-7    * exp( 0.05021/(ts(ji,jj,jk,jp_tem,Kmm)+273.15)*zzz )  &
                                   *3.42031e-3 * rho0_co3sol * rho0_co3sol / max(epsln, ts(ji,jj,jk,jp_sal,Kmm))    &
                                  , epsln )

               zremin_caco3(ji,jj,jk)=  1.d0 / ca_remin_depth  &
                                      *(1.d0 - min(1.d0, co3_ion(ji,jj,jk)/(co3_solubility+epsln) ))

               fcaco3(ji,jj,jk)=(fcaco3km1+jca_uptake(ji,jj,jk)*e3t(ji,jj,jk,Kmm))   &
                                / ( 1.d0+e3t(ji,jj,jk,Kmm) * zremin_caco3(ji,jj,jk))

               jca_reminp(ji,jj,jk)=  (fcaco3km1+jca_uptake(ji,jj,jk)*e3t(ji,jj,jk,Kmm)-fcaco3(ji,jj,jk)) &
                                    / (epsln+e3t(ji,jj,jk,Kmm))

               ! Add remineralization terms to sources
               ! [mol C/m3/s]
               jdic(ji,jj,jk)=jdic(ji,jj,jk)+     jca_reminp(ji,jj,jk)
               ! [mol eq/m3/s]
               jalk(ji,jj,jk)=jalk(ji,jj,jk)+2.d0*jca_reminp(ji,jj,jk)
!DO_CARBON>
               ! Iron
               ! [gC/m3]
               dum3            = (fpop(ji,jj,jk)*c2p*12.011d0/(epsln+wsink))**(0.58)
               ! [mol Fe/m3/s]
               jfe_ads_org(ji,jj,jk) = min (0.5d0, kfe_org*dum3)*feprime(ji,jj,jk)
               ! [mol Fe/m2/s]
               dum4(jk)        =( dum4(jk)+jfe_ads_org(ji,jj,jk) )*e3t(ji,jj,jk,Kmm)
               fpofe(ji,jj,jk) = (fpofekm1 + dum4(jk))/(1.d0+e3t(ji,jj,jk,Kmm)*zremin(ji,jj,jk))
               ! [mol Fe/m3/s]
               jfe_remin(ji,jj,jk) = (fpofekm1+dum4(jk)-fpofe(ji,jj,jk))/(epsln+e3t(ji,jj,jk,Kmm))

               ! Save fPOP and fPOFe at the bottom grid cell to compute bottom fluxes
               ! [mol/m2/s]
               IF (jk==mbkt(ji,jj)) THEN
                  fpop_b (ji,jj) = fpop(ji,jj,jk)
                  fpofe_b(ji,jj) = fpofe(ji,jj,jk)
                  fcaco3_b(ji,jj) = fcaco3(ji,jj,jk)  !DO_CARBON
                  !!! --- AGT --- !!!
                  ! Save fPON at the bottom grid cell to compute bottom fluxes
                  ! [mol/m2/s]
                  IF ( ln_nitro ) THEN
                          fpon_b (ji,jj) = fpon(ji,jj,jk)
                  ENDIF
                  !!! ------ !!!
               ENDIF

               ! Add remineralization terms to trends
               ! [mol P/m3/s]
               jpo4(ji,jj,jk)=jpo4(ji,jj,jk)+(1.d0-phi_dop)*jp_remin(ji,jj,jk)
               jdop(ji,jj,jk)=jdop(ji,jj,jk)+      phi_dop *jp_remin(ji,jj,jk)
               ! [mol Fe/m3/s]
               jfed(ji,jj,jk)=jfed(ji,jj,jk)+jfe_remin(ji,jj,jk)-jfe_ads_org(ji,jj,jk)
   !            IF ( ln_nitro ) THEN 
                       ! Add remineralization terms to trends
                       ! [mol N/m3/s]
                       jno3(ji,jj,jk)=jno3(ji,jj,jk)+(1.d0-phi_dop)*jn_remin(ji,jj,jk)
                       jdon(ji,jj,jk)=jdon(ji,jj,jk)+      phi_dop *jn_remin(ji,jj,jk)
    !           ENDIF

            ENDDO 

            !-----------------------------------------------------------------------
            ! OXYGEN
            !-----------------------------------------------------------------------
            ! Assuming constant P:O ratio.
            ! Optional prevention of negative oxygen (does not conserve ocean 
            ! redox potential) or alternatively it can be allowed to go negative, 
            ! keeping track of an implicit nitrate deficit 
            ! plus sulfate reduction.
            !-----------------------------------------------------------------------

            DO jk=1, jpk
               ! [mol O2/m3]
               foxy = tr(ji,jj,jk,jpOxy_bling,Kmm)
               IF ( (ln_prev_o2lt0) .and. (foxy<oxy_min) ) then
                   ! O2 below hypoxic threshold, no O2 consumption
                   joxy(ji,jj,jk)=0.d0
               ELSE
                   ! [mol O2/m3/s]
                   !!! --- AGT: Include oxygen correction from nitrogen fixation --- !!!
!                   IF ( ln_nitro ) THEN
                           joxy(ji,jj,jk)=-oxy2p*jno3(ji,jj,jk)/n2p - 1.25d0*jn_fix(ji,jj,jk) + oxy2p*jn_fix(ji,jj,jk)/n2p
 !                  ELSE
  !                         joxy(ji,jj,jk)=-oxy2p*jpo4(ji,jj,jk)
 !                  ENDIF
               ENDIF
                             ! Add dissolved inorganic carbon terms from po4 final fluxes
                jdic(ji,jj,jk)=jdic(ji,jj,jk)+jpo4(ji,jj,jk)*c2p
                jalk(ji,jj,jk)=jalk(ji,jj,jk)-jpo4(ji,jj,jk)*n2p

                ! Check for negative values

                tr(ji,jj,jk,jpPO4_bling,Krhs) = MAX(-tr(ji,jj,jk,jpPO4_bling,Kmm),jpo4(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpDOP_bling,Krhs) = MAX(-tr(ji,jj,jk,jpDOP_bling,Kmm),jdop(ji,jj,jk)*rfact)
                !!! AGT - try allowing negative oxygen
                tr(ji,jj,jk,jpOxy_bling,Krhs) = joxy(ji,jj,jk)*rfact
                !tr(ji,jj,jk,jpOxy_bling,Krhs) = MAX(-tr(ji,jj,jk,jpOxy_bling,Kmm),joxy(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpFed_bling,Krhs) = MAX(-tr(ji,jj,jk,jpFed_bling,Kmm),jfed(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpDIC_bling,Krhs) = MAX(-tr(ji,jj,jk,jpDIC_bling,Kmm),jdic(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpalk_bling,Krhs) = MAX(-tr(ji,jj,jk,jpalk_bling,Kmm),jalk(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpNO3_bling,Krhs) = MAX(-tr(ji,jj,jk,jpNO3_bling,Kmm),jno3(ji,jj,jk)*rfact)
                tr(ji,jj,jk,jpDON_bling,Krhs) = MAX(-tr(ji,jj,jk,jpDON_bling,Kmm),jdon(ji,jj,jk)*rfact)

            ENDDO
         ENDDO
      ENDDO

!WRITE(numout,*) '      max carbon trend        = ', maxval(jdic)
!WRITE(numout,*) '      max alkalinity trend    = ', maxval(jalk)
!WRITE(numout,*) '      max phosphate trend     = ', maxval(jpo4)
!WRITE(numout,*) '      max recycling           = ', maxval(jp_recycle)
!WRITE(numout,*) '      max remin               = ', maxval(jp_remin)
!WRITE(numout,*) '      max uptake              = ', maxval(jp_uptake)
!WRITE(numout,*) '      max particulate flux    = ', maxval(jp_pop)
!WRITE(numout,*) '      max dop trend           = ', maxval(jdop)
!WRITE(numout,*) '      max iron trend          = ', maxval(jfed)
!WRITE(numout,*) '      max recycling           = ', maxval(jfe_recycle)
!WRITE(numout,*) '      max remin               = ', maxval(jfe_remin)
!WRITE(numout,*) '      max uptake              = ', maxval(jfe_uptake)
!WRITE(numout,*) '      max organic ads         = ', maxval(jfe_ads_org)
!WRITE(numout,*) '      max inorganic ads       = ', maxval(jfe_ads_inorg)
!WRITE(numout,*) '      max oxygen trend        = ', maxval(joxy)
!WRITE(numout,*) '      max nitrate trend       = ', maxval(jno3)
!WRITE(numout,*) '      max n recycling         = ', maxval(jn_recycle)
!WRITE(numout,*) '      max n remin             = ', maxval(jn_remin)
!WRITE(numout,*) '      max n uptake            = ', maxval(jn_uptake)
!WRITE(numout,*) '      max n fixation          = ', maxval(jn_fix)
!WRITE(numout,*) '      max particulate n flux  = ', maxval(jn_pon)
!WRITE(numout,*) '      max don trend           = ', maxval(jdon)

!WRITE(numout,*) '      min carbon trend        = ', minval(jdic)
!WRITE(numout,*) '      min alkalinity trend    = ', minval(jalk)
!WRITE(numout,*) '      min phosphate trend     = ', minval(jpo4)
!WRITE(numout,*) '      min recycling           = ', minval(jp_recycle)
!WRITE(numout,*) '      min remin               = ', minval(jp_remin)
!WRITE(numout,*) '      min uptake              = ', minval(jp_uptake)
!WRITE(numout,*) '      min particulate flux    = ', minval(jp_pop)
!WRITE(numout,*) '      min dop trend           = ', minval(jdop)
!WRITE(numout,*) '      min iron trend          = ', minval(jfed)
!WRITE(numout,*) '      min recycling           = ', minval(jfe_recycle)
!WRITE(numout,*) '      min remin               = ', minval(jfe_remin)
!WRITE(numout,*) '      min uptake              = ', minval(jfe_uptake)
!WRITE(numout,*) '      min organic ads         = ', minval(jfe_ads_org)
!WRITE(numout,*) '      min inorganic ads       = ', minval(jfe_ads_inorg)
!WRITE(numout,*) '      min oxygen trend        = ', minval(joxy)
!WRITE(numout,*) '      min nitrate trend       = ', minval(jno3)
!WRITE(numout,*) '      min n recycling         = ', minval(jn_recycle)
!WRITE(numout,*) '      min n remin             = ', minval(jn_remin)
!WRITE(numout,*) '      min n uptake            = ', minval(jn_uptake)
!WRITE(numout,*) '      min n fixation          = ', minval(jn_fix)
!WRITE(numout,*) '      min particulate n flux  = ', minval(jn_pon)
!WRITE(numout,*) '      min don trend           = ', minval(jdon)

!WRITE(numout,*) '      max fpop_b              = ', maxval(fpop_b)
!WRITE(numout,*) '      max fpon_b              = ', maxval(fpon_b)

!WRITE(numout,*) '      min fpop_b              = ', minval(fpop_b)
!WRITE(numout,*) '      min fpon_b              = ', minval(fpon_b)

      ! Add dissolved inorganic carbon terms from po4 final fluxes
!      jdic(:,:,:)=jdic(:,:,:)+jpo4(:,:,:)*c2p
!      jalk(:,:,:)=jalk(:,:,:)-jpo4(:,:,:)*n2p

 !     tr(:,:,:,jpPO4_bling,Krhs) = tr(:,:,:,jpPO4_bling,Krhs) + jpo4(:,:,:)*rfact
 !     tr(:,:,:,jpDOP_bling,Krhs) = tr(:,:,:,jpDOP_bling,Krhs) + jdop(:,:,:)*rfact
 !     tr(:,:,:,jpFed_bling,Krhs) = tr(:,:,:,jpFed_bling,Krhs) + jfed(:,:,:)*rfact
 !     tr(:,:,:,jpOxy_bling,Krhs) = tr(:,:,:,jpOxy_bling,Krhs) + joxy(:,:,:)*rfact

      !!! --- AGT: Add nitrogen terms --- !!!
!      IF ( ln_nitro ) THEN
!              tr(:,:,:,jpNO3_bling,Krhs) = tr(:,:,:,jpNO3_bling,Krhs) + jno3(:,:,:)*rfact
!              tr(:,:,:,jpDON_bling,Krhs) = tr(:,:,:,jpDON_bling,Krhs) + jdon(:,:,:)*rfact
!      ENDIF
      !!! ------ !!!

      !jdic(:,:,:)=0.e0
!      tr(:,:,:,jpDIC_bling,Krhs) = tr(:,:,:,jpDIC_bling,Krhs) + jdic(:,:,:)*rfact

      !jalk(:,:,:)=0.e0
!      tr(:,:,:,jpalk_bling,Krhs) = tr(:,:,:,jpalk_bling,Krhs) + jalk(:,:,:)*rfact

      ! Feb 10, 2017, xianmin force last level to be zero
!      tr(:,:,jpk,jp_blg0:jp_blg1,Krhs)=0.0_wp

      !test if concentrations fall below 0
!      xnegtr(:,:,:) = 0.e0
!      DO jn = jp_blg0, jp_blg1
!         DO jk = 1, jpk
!            DO jj = 1, jpj
!               DO ji = 1, jpi
!                  IF( ( tr(ji,jj,jk,jn,Kmm) + tr(ji,jj,jk,jn,Krhs) ) < 0.e0 ) THEN 
!                     ztra             = ABS(  ( tr(ji,jj,jk,jn,Kmm) - rtrn ) &
!                                            / ( tr(ji,jj,jk,jn,Krhs) + rtrn ) )
!                     xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
!                      xnegtr(ji,jj,jk) = - tr(ji,jj,jk,jn,Kmm) -tr(ji,jj,jk,jn,Krhs)
!                  ENDIF
!              END DO
!            END DO
!         END DO
!         tr(:,:,:,jn,Krhs) = tr(:,:,:,jn,Krhs) + xnegtr(:,:,:)
!      END DO

      ! Prgonostic tracer fields
      DO jn = jp_blg0, jp_blg1
   !      WRITE(numout,*) '   ==>>>   AGT: updating tracer conc. '
         !write(*,'(I3,2(1X,E11.4))') jp_bling0, trn(ji,jj,jk,jn),tra(ji,jj,jk,jn)
!         tr(:,:,:,jn,Kmm) = tr(:,:,:,jn,Kmm) + xnegtr(:,:,:) * tr(:,:,:,jn,Krhs)
!         tr(:,:,:,jn,Kmm) = tr(:,:,:,jn,Kmm) + xnegtr(:,:,:) * tr(:,:,:,jn,Krhs) &
!                                 + (1 - xnegtr(:,:,:))*(1.e-11 - tr(:,:,:,jn,Kmm))
         tr(:,:,:,jn,Kmm) = tr(:,:,:,jn,Kmm) + tr(:,:,:,jn,Krhs)
      END DO

      ! Feb 10, 2017, xianmin force last level to be zero
!      tr(:,:,jpk,jp_blg0:jp_blg1,Kmm)=0.0_wp

      !DO jn = jp_bling0, jp_bling1
      !   trn(:,:,:,jn) = trn(:,:,:,jn) + tra(:,:,:,jn)
      !   DO jk = 1, jpk
      !      DO jj = 1, jpj
      !         DO ji = 1, jpi
      !           trn(ji,jj,jk,jn)=MAX( 0.e0, trn(ji,jj,jk,jn) )
      !         END DO
      !      END DO
      !   END DO
      !END DO

      ! Copy new arrays to trb (tracer fields before) and set tra to zero
      ! to compute tracer gradients with tracer fields after ecological forcing
!      tr(:,:,:,jp_blg0:jp_blg1,Krhs) = 0.e0

      ! add external fluxes
      CALL trc_ext_bling (kt,Kbb,Kmm,Krhs)

      DO jn=jp_blg0, jp_blg1
         CALL lbc_lnk( 'trcsms_blingv0', tr(:,:,:,jn,Kmm), 'T', 1.0_wp )
         !CALL lbc_lnk( trb(:,:,:,jn), 'T', 1.0_wp )  ! Feb 17, 2017, xianmin
         !CALL lbc_lnk( tra(:,:,:,jn), 'T', 1.0_wp )  ! Feb 17, 2017, xianmin
      ENDDO
      ! Copy to trb to use BLING updated tracer values to compute transport
      ! trends
      DO jn=jp_blg0, jp_blg1
!         tr(:,:,:,jn,Krhs)=tr(:,:,:,jn,Kmm)
         tr(:,:,:,jn,Kbb)=tr(:,:,:,jn,Kmm)
      ENDDO


      IF ( ln_bling_mass) THEN

        IF( kt == nittrc000 ) THEN 

          if (lwp) then ! xhu to enable text output
             CALL ctl_opn( numsms ,  'sms.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
             CALL ctl_opn( numsms2, 'sms2.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
          endif
        ENDIF

        ! xhu
        sumrecycle=0.0_wp
        sumremin=0.0_wp
        sumuptake=0.0_wp
        sumuptake2=0.0_wp
        sumremin2=0.0_wp
        sumpo4=0.0_wp
        sumdop=0.0_wp
        sumdic=0.0_wp
        DO jk=1,jpk
           sumrecycle= sumrecycle+ glob_sum ('trcsms_blingv0',jp_recycle(:,:,jk)*cvol(:,:,jk))
           sumremin  = sumremin  + glob_sum ('trcsms_blingv0',  jp_remin(:,:,jk)*cvol(:,:,jk))
           sumuptake = sumuptake + glob_sum ('trcsms_blingv0',-jp_uptake(:,:,jk)*cvol(:,:,jk))
           sumuptake2 = sumuptake + glob_sum ('trcsms_blingv0',(jp_uptake(:,:,jk)-jp_pop(:,:,jk))*cvol(:,:,jk))
           sumremin2 = sumremin2 + glob_sum ('trcsms_blingv0',jp_remin(:,:,jk)*cvol(:,:,jk))
           sumpo4    = sumpo4    + glob_sum ('trcsms_blingv0',jpo4(:,:,jk)*cvol(:,:,jk))
           sumdop    = sumdop    + glob_sum ('trcsms_blingv0',jdop(:,:,jk)*cvol(:,:,jk))
           sumdic    = sumdic    + glob_sum ('trcsms_blingv0',jdic(:,:,jk)*cvol(:,:,jk))
        ENDDO
        sumrecycle=sumrecycle*rfact
        sumremin  =sumremin  *(1.d0-phi_dop)*rfact
        sumuptake =sumuptake * rfact
        sumuptake2 =sumuptake2 *phi_dop * rfact
        sumremin2 =sumremin2 * phi_dop * rfact
        sumpo4    =sumpo4*rfact
        sumdop    =sumdop*rfact
        sumdic    =sumdic*rfact
        if (lwp) then
           WRITE(UNIT=numsms,FMT='(i10,5(3x,e18.10))')  kt, sumrecycle, sumremin, &
                                                        sumuptake, sumuptake2,sumremin2
           WRITE(UNIT=numsms2,FMT='(i10,3(3x,e18.10))')  kt, sumpo4, sumdop, sumdic
        endif

      ENDIF

      
      !IF (lk_iomput) THEN
         !CALL iom_put(      "expkT"  ,        expkT(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(   "irr_inst"  ,     irr_inst(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(    "irr_mix"  ,      irr_mix(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(       "irrk"  ,         irrk(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(       "pc_m"  ,         pc_m(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(         "mu"  ,           mu(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(  "biomass_p"  ,    biomass_p(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(       "fpop"  ,         fpop(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(     "zremin"  ,       zremin(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(     "def_fe"  ,       def_fe(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(    "feprime"  ,      feprime(:,:,:)*tmask(:,:,:) )
         !CALL iom_put( "kfe_eq_lig"  ,   kfe_eq_lig(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(      "fpofe"  ,        fpofe(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(   "frac_pop"  ,    frac_pop(:,:,:)*tmask(:,:,:) )

         !CALL iom_put(     "jp_pop"  ,       jp_pop(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(  "jp_uptake"  ,    jp_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(   "jp_remin"  ,     jp_remin(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put( "jp_recycle"  ,   jp_recycle(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                           
         !CALL iom_put( "jfe_uptake"  ,   jfe_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                          
         !CALL iom_put(  "jfe_remin"  ,    jfe_remin(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                          
         !CALL iom_put( "jfe_recycle" ,  jfe_recycle(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                          
         !CALL iom_put( "jfe_ads_org" ,  jfe_ads_org(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                           
         !CALL iom_put("jfe_ads_inorg",jfe_ads_inorg(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )                          
         !CALL iom_put(        "jpo4" ,         jpo4(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(        "jdop" ,         jdop(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(        "jfed" ,         jfed(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(        "joxy" ,         joxy(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )

!DO CARBON<-
         !CALL iom_put( "jca_uptake" ,  jca_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put( "jca_reminp" ,  jca_reminp(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(      "jdic"  ,        jdic(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(      "jalk"  ,        jalk(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:) )

         !CALL iom_put(    "frac_lg"  ,      frac_lg(:,:,:)*tmask(:,:,:) )
         !CALL iom_put("zremin_caco3" , zremin_caco3(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(     "fcaco3"  ,       fcaco3(:,:,:)*tmask(:,:,:) )
         !CALL iom_put(  "pco2_surf"  ,    pco2_surf(:,:)  *tmask(:,:,1) )
!DO CARBON>

      !ENDIF

      IF ( .NOT. lk_iomput )  THEN
        trc2d( :,:,jp_bling0_2d ) = pco2_surf(:,:)*tmask(:,:,1)

        trc3d( :,:,:,jp_bling0_3d    ) =    chl_bling(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+1  ) =         jpo4(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+2  ) =         jdop(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+3  ) =       jp_pop(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+4  ) =    jp_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+5  ) =   jp_recycle(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+6  ) =     jp_remin(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+7  ) =         jfed(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+8  ) =   jfe_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+9  ) =  jfe_recycle(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+10 ) =    jfe_remin(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+11 ) =  jfe_ads_org(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+12 ) =jfe_ads_inorg(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+13 ) =         joxy(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+14 ) =         jdic(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+15 ) =   jca_uptake(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+16 ) =   jca_reminp(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+17 ) =         jalk(:,:,:)*e3t_0(:,:,:)*tmask(:,:,:)

        trc3d( :,:,:,jp_bling0_3d+18 ) =        fpop(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+19 ) =       fpofe(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+20 ) =       expkT(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+21 ) =    irr_inst(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+22 ) =     irr_mix(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+23 ) =        irrk(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+24 ) =        pc_m(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+25 ) =          mu(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+26 ) =   biomass_p(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+27 ) =      zremin(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+28 ) =      def_fe(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+29 ) =     feprime(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+30 ) =  kfe_eq_lig(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+31 ) =     frac_lg(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+32 ) =zremin_caco3(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+33 ) =      fcaco3(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+34 ) =     co3_ion(:,:,:)*tmask(:,:,:)
        trc3d( :,:,:,jp_bling0_3d+35 ) =      htotal(:,:,:)*tmask(:,:,:)
      ENDIF

      DEALLOCATE ( rho )
      DEALLOCATE ( expkT )
      DEALLOCATE ( irr_inst )
      DEALLOCATE ( irr_mix )
      DEALLOCATE ( irrk )
      DEALLOCATE ( pc_m )
      DEALLOCATE ( mu )
      DEALLOCATE ( def_fe )
      DEALLOCATE ( feprime )
      DEALLOCATE ( kfe_eq_lig )
      DEALLOCATE ( fpofe )
      DEALLOCATE ( jp_pop )
      DEALLOCATE ( fpop )
      DEALLOCATE ( zremin )
      DEALLOCATE ( jp_uptake )
      DEALLOCATE ( jp_remin )
      DEALLOCATE ( jp_recycle )
      DEALLOCATE ( jfe_uptake )
      DEALLOCATE ( jfe_remin )
      DEALLOCATE ( jfe_recycle )
      DEALLOCATE ( jfe_ads_inorg )
      DEALLOCATE ( jfe_ads_org )
      DEALLOCATE ( jpo4 )
      DEALLOCATE ( jdop )
      DEALLOCATE ( jfed )
      DEALLOCATE ( joxy )
      DEALLOCATE ( jdic )
      DEALLOCATE ( jalk )
      DEALLOCATE ( jca_uptake )
      DEALLOCATE ( jca_reminp )
      DEALLOCATE ( fcaco3 )
      DEALLOCATE ( zremin_caco3 )
      DEALLOCATE ( frac_lg )
      DEALLOCATE ( pco2_surf )
      DEALLOCATE ( dum4 )
      DEALLOCATE ( xnegtr )

      !!! --- AGT --- !!!
    !  IF ( ln_nitro ) THEN
              DEALLOCATE ( jno3 )
              DEALLOCATE ( jdon )
              DEALLOCATE ( jn_pon )
              DEALLOCATE ( fpon )
              DEALLOCATE ( jn_uptake )
              DEALLOCATE ( jn_fix )
              DEALLOCATE ( jn_remin )
              DEALLOCATE ( jn_recycle )
              DEALLOCATE ( pc_m_diaz )
              DEALLOCATE ( mu_diaz )
     ! ENDIF
      !!! ------ !!!

      IF( ln_timing == 1 )  CALL timing_stop('trc_sms_bling')
      !
   END SUBROUTINE trc_sms_bling

   SUBROUTINE trc_sms_bling_mass_conserv (kt, Kbb, Kmm, Krhs)

      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level index
      REAL(wp) :: sum_phosp, sum_fed, sum_oxy, sum_carbon
      INTEGER  :: jk

      
      !!-------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN 
        if (lwp) then ! to enable text output
           CALL ctl_opn( numphp, 'php.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE.)
           CALL ctl_opn( numfed, 'fed.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE.)
           CALL ctl_opn( numoxy, 'oxy.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE.)
           CALL ctl_opn( numcar, 'car.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE.)
           WRITE(numphp,9500) kt,  areatot
           WRITE(numfed,9500) kt,  areatot
           WRITE(numoxy,9500) kt,  areatot
           WRITE(numcar,9500) kt,  areatot
        endif
      ENDIF

         ! total mass of phosphate
         ! xhu
         sum_phosp=0.0_wp
         sum_fed=0.0_wp
         sum_oxy=0.0_wp
         sum_carbon=0.0_wp
         DO jk=1,jpk
            ! PO4+DOP is conserved over the whole domain
            sum_phosp = sum_phosp + glob_sum ( 'trcsms_blingv0',( tr(:,:,jk,jpPO4_bling,Kmm) + tr(:,:,jk,jpDOP_bling,Kmm) )*cvol(:,:,jk)  )

            ! DIC+DOP*106 is conserved over the whole domain w/o air-sea fluxes
            sum_carbon = sum_carbon + glob_sum ( 'trcsms_blingv0',( tr(:,:,jk,jpDIC_bling,Kmm) + tr(:,:,jk,jpDOP_bling,Kmm)*106 )*cvol(:,:,jk)  )

            ! Non-conservative tracers
            sum_fed   = sum_fed   + glob_sum ( 'trcsms_blingv0',  tr(:,:,jk,jpFed_bling,Kmm)*cvol(:,:,jk)  )
            sum_oxy   = sum_oxy   + glob_sum ( 'trcsms_blingv0',  tr(:,:,jk,jpOxy_bling,Kmm)*cvol(:,:,jk)  )
         ENDDO
         
      IF ( lwp ) THEN
         WRITE(numphp,9500) kt,  sum_phosp, areatot
         WRITE(numfed,9500) kt,  sum_fed, areatot
         WRITE(numoxy,9500) kt,  sum_oxy, areatot
         WRITE(numcar,9500) kt,  sum_carbon, areatot
      ENDIF

9500  FORMAT(i10,2(3x,e18.10))

   END SUBROUTINE trc_sms_bling_mass_conserv
   
   SUBROUTINE trc_sms_init_bling

        INTEGER  :: numbio, kt 

         !!! --- AGT: Add option to read in initial biomass, assumed in mol/kg...  --- !!!
        ALLOCATE( sf_biomass_init(1) )
        CALL fld_fill( sf_biomass_init, (/ sn_biomass_init /), cn_dir_biomass_init, 'trc_sms_init_bling', 'Initial biomass ', 'namblingprod' )
        ALLOCATE( sf_biomass_init(1)%fnow(jpi,jpj,jpk)   )
        ! Doesn't work for time interpolation !
        !CALL fld_read( kt, 1, sf_biomass_init ) 
        !biomass_p(:,:,:) = sf_biomass_init(1)%fnow(:,:,:)     
        CALL iom_open (  TRIM( sn_biomass_init%clname ) , numbio ) 
        CALL iom_get( numbio, jpdom_global, TRIM( sn_biomass_init%clvar ), biomass_p(:,:,:) )
   
   END SUBROUTINE trc_sms_init_bling

END MODULE trcsms_blingv0
