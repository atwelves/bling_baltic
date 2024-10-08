MODULE vars_bling
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE par_oce    ! ocean parameters
   USE par_trc    ! passive tracer parameters
   USE lib_mpp
   USE fldread    ! defines FLD_N structure var
   USE iom
   IMPLICIT NONE
   PUBLIC

   PUBLIC bling_alloc

   !! Diagnostic variables
   !! ----------------------

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: chl_bling, irr_mem, biomass_p
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: htotal, co3_ion  !DO_CARBON
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fpop_b,fpofe_b,fcaco3_b
   !!! --- AGT --- !!!
   !IF ( ln_nitro ) THEN
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: biomass_p_diaz
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fpon_b
   !ENDIF
   !!! ------ !!!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: co2_csurf, co2_alpha

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: coast_bling
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: dust_bling     ! dust fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: patm_bling    ! atmospheric pressure at kt [N/m2]

#if defined key_rnf_nutrients
   ! add nutrients in the river runoff 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: rnfpo4, rnfdop, rnffed, rnfdic, rnfalk
#endif

   !! Numerical parameter
   !! ----------------------
   REAL(wp), PUBLIC :: rfact
   REAL(wp), PUBLIC, PARAMETER :: epsln=1.0e-30
   ! Minimum chl value allowed for numerical stability
   REAL(wp), PUBLIC, PARAMETER :: chl_min=1.e-8 !mg chl/kg]

   !! Stochiometric ratios
   !! ----------------------
   REAL(wp) :: c2n     !: redfield ratio c:n                       (NAMELIST)
   REAL(wp) :: c2p
   REAL(wp) :: oxy2p
   REAL(wp) :: n2p  !DO_CARBON
   REAL(wp) :: ca2p !DO_CARBON
   REAL(wp) :: rho0_co3sol !DO_CARBON
   REAL(wp) :: ca_remin_depth !DO_CARBON
   REAL(wp) :: htotal_scale_lo, htotal_scale_hi

   !! Production parameters
   !! ----------------------
   REAL(wp) :: pc_0
   REAL(wp) :: kappa_eppley
   REAL(wp) :: kpo4
   REAL(wp) :: kfe          
   REAL(wp) :: fe2p_max     
   REAL(wp) :: kfe2p_up     
   REAL(wp) :: def_fe_min   
   REAL(wp) :: thetamax_lo
   REAL(wp) :: thetamax_hi
   REAL(wp) :: alpha_max
   REAL(wp) :: alpha_min
   REAL(wp) :: resp_frac
   REAL(wp) :: p_star
   REAL(wp) :: lambda0
   REAL(wp) :: gam_biomass
   !REAL(wp) ::
   !!! --- AGT --- !!!
   !IF ( ln_nitro ) THEN
   LOGICAL  :: ln_nitro
   REAL(wp) :: kno3
   !ENDIF 
   ! Initial biomass
   TYPE(FLD_N) :: sn_biomass_init
   CHARACTER(len=100) :: cn_dir_biomass_init 

   !! Optical parameters                                
   !! ------------------                                
   REAL(wp) ::   xkr0     !: water coefficient absorption in red      (NAMELIST)
   REAL(wp) ::   xkb0     !: water coefficient absorption in green    (NAMELIST)
   REAL(wp) ::   xkrp     !: pigment coefficient absorption in red    (NAMELIST)
   REAL(wp) ::   xkbp     !: pigment coefficient absorption in green  (NAMELIST)
   REAL(wp) ::   xlr      !: exposant for pigment absorption in red   (NAMELIST)
   REAL(wp) ::   xlb      !: exposant for pigment absorption in green (NAMELIST)
   REAL(wp) ::   rpig     !: chla/chla+phea ratio                     (NAMELIST)
   REAL(wp) ::   rcchl    !: ???                                              
   REAL(wp) ::gam_irr_mem !: photoadaptation time constant            (NAMELIST) 

   !! Remineralization parameters
   !! ---------------------------
   REAL(wp) :: wsink0_z  
   REAL(wp) :: wsink0 
   REAL(wp) :: wsink_acc
   REAL(wp) :: koxy
   REAL(wp) :: remin_min
   REAL(wp) :: phi_dop
   REAL(wp) :: phi_sm
   REAL(wp) :: phi_lg
   REAL(wp) :: kappa_remin
   REAL(wp) :: gamma_dop
   REAL(wp) :: gamma_pop
   !!! --- AGT: Add variables for nitrogen in bottom level --- !!!
   REAL(wp) :: fden
   REAL(wp) :: fbur
   !...and benthic oxygen consumption....
   REAL(wp) :: fben
   ! Also, add variables for organic matter in rivers...
   REAL(wp) :: river_dop
   REAL(wp) :: river_don
   !!! ------ !!!

   !! Air-sea interaction parameters
   !! ------------------------------
   REAL(wp)    :: a_0, a_1, a_2, a_3, a_4, a_5
   REAL(wp)    :: b_0, b_1, b_2, b_3, b_4
   REAL(wp)    :: c_0
   REAL(wp)    :: a_1_o2, a_2_o2, a_3_o2, a_4_o2, a_5_o2
   REAL(wp)    :: a_1_co2, a_2_co2, a_3_co2, a_4_co2, a_5_co2
   REAL(wp)    :: pco2_cte=278._wp    ! units ppm
   REAL(wp)    :: pco2(142)            ! length of vector as defined in trcini_blingv0
   REAL(wp)    :: po2_atm=0.20946_wp  ! units atm
   LOGICAL     :: ln_patm_bling
   LOGICAL     :: ln_Pa2atm=.FALSE.   ! Apr 25, 2016 (xhu): convert Pa to atm, xh
   TYPE(FLD_N) :: sn_patm_bling
   CHARACTER(len=100) :: cn_dir_patm_bling
                                                        
   !! Iron parameters                  
   !! ------------------                                
   REAL(wp) :: kfe_eq_lig_irr
   REAL(wp) :: kfe_eq_lig_femin
   REAL(wp) :: kfe_eq_lig_max
   REAL(wp) :: kfe_eq_lig_min
   REAL(wp) :: felig_bkg
   REAL(wp) :: kfe_inorg
   REAL(wp) :: kfe_org
   REAL(wp) :: oxy_min
   LOGICAL  :: ln_prev_o2lt0
   LOGICAL  :: ln_dust_bling
   TYPE(FLD_N) ::   sn_dust_bling
   CHARACTER(len=100) ::  cn_dir_bling

   !! Budget parameters
   !! -----------------
   LOGICAL :: ln_bling_mass, ln_bling_ext

   !! Cumulative sums
   !! ---------------
   REAL(wp) :: sumtffed, sumtfoxy, sumbfpo4, sumbffed, sumbfoxy

   !REAL(wp) :: 
CONTAINS

   INTEGER FUNCTION bling_alloc()

#if defined key_rnf_nutrients
     INTEGER :: ierr(7)
#else
     INTEGER :: ierr(6)
#endif

     ierr(:)=0

     ! optical model
     ALLOCATE( chl_bling(jpi,jpj,jpk), irr_mem(jpi,jpj,jpk) &
               ,biomass_p(jpi,jpj,jpk), STAT=ierr(1) )

     ! bottom fluxes
     ALLOCATE(  fpop_b(jpi,jpj),  fpofe_b(jpi,jpj), fcaco3_b(jpi,jpj) , STAT=ierr(3) )
     
     !!! --- AGT --- !!!
     IF ( ln_nitro ) THEN
             ALLOCATE(  biomass_p_diaz(jpi,jpj,jpk) )
             ALLOCATE(  fpon_b(jpi,jpj) )
     ENDIF
     !!! ------ !!!

     ! dust fluxes
     ALLOCATE( dust_bling(jpi,jpj), coast_bling(jpi,jpj,jpk), STAT=ierr(4) )

     ! carbon cycle DO_CARBON
     ALLOCATE(  co3_ion  (jpi,jpj,jpk), htotal   (jpi,jpj,jpk), STAT=ierr(5) )
     ALLOCATE(  co2_csurf(jpi,jpj)    , co2_alpha(jpi,jpj)    , patm_bling(jpi,jpj), STAT=ierr(6) )

#if defined key_rnf_nutrients
      ALLOCATE( rnfpo4(jpi,jpj), rnfdop(jpi, jpj), rnfdic(jpi,jpj), rnfalk(jpi,jpj), rnffed(jpi,jpj), STAT=ierr(7))   
#endif
    ! WRITE(numout,*) '   ==>>>   AGT: irradiance memory allocated '
     bling_alloc=MAXVAL(ierr)

     IF( lk_mpp )   CALL mpp_sum ( 'vars_bling',bling_alloc )

   END FUNCTION bling_alloc

   !!======================================================================
END MODULE vars_bling
