MODULE trcnam_bling
   !!======================================================================
   !!                      ***  MODULE trcnam_bling  ***
   !! TOP :   initialisation of some run parameters for LOBSTER bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_bling'   :                                       BLINGv0 model
   !!----------------------------------------------------------------------
   !! trc_nam_bling      : BLINGv0 model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE iom
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_bling   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_bling.F90 -1   $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_bling
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_bling  ***  
      !!
      !! ** Purpose :   read BLINGv0 namelist
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER :: numnatg, jl, jn

      TYPE(DIAG), DIMENSION(jp_bling_2d )  :: blidia2d
      TYPE(DIAG), DIMENSION(jp_bling_3d )  :: blidia3d
      TYPE(DIAG), DIMENSION(jp_bling_trd ) :: blidiabio

      !NAMELIST/namblingpatm/   ln_patm_bling, sn_patm_bling, cn_dir_patm_bling
      ! Apr 25, 2016 (xhu): add ln_Pa2atm to convert Pa to atm
      NAMELIST/namblingpatm/   ln_patm_bling, sn_patm_bling, cn_dir_patm_bling, ln_Pa2atm

      NAMELIST/namblingrat/    c2n, c2p, oxy2p, n2p, ca2p, rho0_co3sol, ca_remin_depth, htotal_scale_lo, htotal_scale_hi
      NAMELIST/namblingopt/    xkr0, xkb0, xkrp, xkbp, xlr, xlb, rpig, rcchl, gam_irr_mem
      NAMELIST/namblingprod/   pc_0, kappa_eppley, kpo4, kfe, fe2p_max, kfe2p_up, def_fe_min, thetamax_lo, thetamax_hi
      NAMELIST/namblingprod/   alpha_min, alpha_max, resp_frac, p_star, lambda0, gam_biomass, ln_nitro, kno3
      NAMELIST/namblingprod/   sn_biomass_init, cn_dir_biomass_init
      NAMELIST/namblingremin/  wsink0_z, wsink0, wsink_acc, koxy, remin_min, phi_dop, phi_sm
      NAMELIST/namblingremin/  phi_lg, kappa_remin, gamma_dop, gamma_pop
      !!! --- AGT --- !!!
      NAMELIST/namblingremin/  fden, fbur, fben
      NAMELIST/namblingremin/  river_dop, river_don
      !!! ------ !!!
      NAMELIST/namblingairsea/ a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
      NAMELIST/namblingairsea/ a_1_o2, a_2_o2, a_3_o2, a_4_o2, a_5_o2
      NAMELIST/namblingairsea/ a_1_co2, a_2_co2, a_3_co2, a_4_co2, a_5_co2
      NAMELIST/namblingiron/   oxy_min, kfe_eq_lig_max, kfe_eq_lig_min, kfe_eq_lig_irr, kfe_eq_lig_femin
      NAMELIST/namblingiron/   felig_bkg, kfe_inorg, kfe_org, ln_prev_o2lt0
      NAMELIST/namblingiron/   cn_dir_bling, sn_dust_bling, ln_dust_bling
      NAMELIST/namblingbudget/ ln_bling_mass, ln_bling_ext               ! additional diagnostics
      NAMELIST/namblingdia/    blidia3d, blidia2d                        ! additional diagnostics
      NAMELIST/namblingdbi/    blidiabio                                 ! additional diagnostics
      !
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_bling : read BLINGv0 namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

      CALL ctl_opn( numnatg, 'namelist_blingv0', 'OLD', 'FORMATTED','SEQUENTIAL', -1, numout, .FALSE. )

      REWIND (numnatg)
      READ   (numnatg, namblingpatm)

      ! namblingrat  : Stochiometric ratios
      REWIND (numnatg)
      READ   (numnatg, namblingrat)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingrat:'
         WRITE(numout,*) '                           c2n     = ', c2n
         WRITE(numout,*) '                           c2p     = ', c2p
         WRITE(numout,*) '                           oxy2p   = ', oxy2p
         WRITE(numout,*) '    CaCO3    to Phosphorus ratio (of small phytoplankton)   = ', ca2p
         WRITE(numout,*) '    Nitrogen to Phosphorus ratio                            = ', n2p
         WRITE(numout,*) '    Constant rho value used to compute co3 solubility       = ', rho0_co3sol
         WRITE(numout,*) '   = ', htotal_scale_lo
         WRITE(numout,*) '   = ', htotal_scale_hi
         WRITE(numout,*) '   = ', ca_remin_depth
      ENDIF

      ! namblingprod : Production parameters
      REWIND (numnatg)
      READ   (numnatg, namblingprod)

      pc_0       =    pc_0/86400.d0
      lambda0    =    lambda0/86400.d0
      gam_biomass=gam_biomass/86400.d0

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingprod:'
         WRITE(numout,*) '    max growth rate at 0 C                          pc_0        = ', pc_0
         WRITE(numout,*) '    T dependence of growth                          kappa_epply = ', kappa_eppley 
         WRITE(numout,*) '    PO4 uptake half-saturation cte                  kpo4        = ', kpo4
         WRITE(numout,*) '    dissolved Fe uptake half-saturation cte         kfe         = ', kfe
         WRITE(numout,*) '    minimum value for iron deficiency term          def_fe_min  = ', def_fe_min
         WRITE(numout,*) '    maximum Fe:P uptake ratio                       fe2p_up_max = ', fe2p_max 
         WRITE(numout,*) '    half-saturation cellular Fe:P                   kfe2p_up    = ', kfe2p_up
         WRITE(numout,*) '    max Chl:C ratio, extreme iron limitation        thetamax_lo = ', thetamax_lo
         WRITE(numout,*) '    max Chl:C ratio, abundant iron                  thetamax_hi = ', thetamax_hi
         WRITE(numout,*) '    quantum yield under low light, iron limited     alpha_min   = ', alpha_min 
         WRITE(numout,*) '    quantum yield under low light, abundant iron    alpha_max   = ', alpha_max
         WRITE(numout,*) '    fraction of gross production respirated         resp_frac   = ', resp_frac
         WRITE(numout,*) '    pivotal phytoplankton biomass                   p_star      = ', p_star
         WRITE(numout,*) '    carbon-specific phytoplankton mortality rate    lambda0     = ', lambda0
         WRITE(numout,*) '    biomass adjustment time constant                gam_biomass = ', gam_biomass
      ENDIF

      ! namblingopt : Optical parameters
      REWIND (numnatg)
      READ   (numnatg, namblingopt)

      gam_irr_mem=gam_irr_mem/86400.d0

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingopt:'
         WRITE(numout,*) '    green   water absorption coeff                       xkg0  = ', xkb0
         WRITE(numout,*) '    red water absorption coeff                           xkr0  = ', xkr0
         WRITE(numout,*) '    pigment red absorption coeff                         xkrp  = ', xkrp
         WRITE(numout,*) '    pigment green absorption coeff                       xkgp  = ', xkbp
         WRITE(numout,*) '    green chl exposant                                   xlg   = ', xlb
         WRITE(numout,*) '    red   chl exposant                                   xlr   = ', xlr
         WRITE(numout,*) '    chla/chla+phea ratio                                 rpig  = ', rpig
         WRITE(numout,*) '    Photoadaptation time constant                 gam_irr_mem  = ', gam_irr_mem
      ENDIF

      ! namblingremin : Remineralization parameters
      REWIND (numnatg)
      READ   (numnatg, namblingremin)

      wsink0   =wsink0/86400.d0
      wsink_acc=wsink_acc/86400.d0
      gamma_dop=gamma_dop/365.25/86400.d0
      gamma_pop=gamma_pop/86400.d0
      !!! --- AGT --- !!!
      fden     =fden/86400.d0
      fbur     =fbur/86400.d0
      fben     =fben/86400.d0
      !!! ------ !!!

      IF (lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingremin:'
         WRITE(numout,*) '    depth at which sinking rate starts increasing    wsink0_z    =', wsink0_z 
         WRITE(numout,*) '    initial sinking rate                             wsink0      =', wsink0
         WRITE(numout,*) '    accerelation rate of sinking with depth          wsink_acc   =', wsink_acc
         WRITE(numout,*) '    half saturation const for aerobic respiration    koxy        =', koxy
         WRITE(numout,*) '    minimum anaerobic respiration rate               remin_min   =', remin_min 
         WRITE(numout,*) '    fraction of non-particulate uptake to DOM        phi_dop     =', phi_dop
         WRITE(numout,*) '    detritus production by small phyto               phi_sm      =', phi_sm
         WRITE(numout,*) '    detritus production by large phyto               phi_lg      =', phi_lg
         WRITE(numout,*) '    T dependence of particulate production           kappa_remin =', kappa_remin
         WRITE(numout,*) '    decay timescale of DOM                           gamma_dop   =', gamma_dop
         WRITE(numout,*) '    remineralization rate of sinking POM             gamma_pop   =', gamma_pop
         !!! --- AGT: parameters for nitrogen in sediment --- !!!
         WRITE(numout,*) '    rate of benthic denitrification                  fden        =', fden
         WRITE(numout,*) '    rate of sediment burial                          fbur        =', fbur
         ! Oxygen consumption in benthos
         WRITE(numout,*) '    rate of benthic oxygen consumption               fben        =', fben
         ! Also, global concentrations for DOM in river inflows...
         WRITE(numout,*) '    concentration of DOP in river runoff             river_dop   =', river_dop
         WRITE(numout,*) '    concentration of DON in river runoff             river_don   =', river_don
         !!! ------ !!!
         !WRITE(numout,*) ' '
      ENDIF

      ! namblingairsea : Air-sea interaction parameters
      REWIND (numnatg)
      READ   (numnatg, namblingairsea)

      IF (lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingairsea:'
         WRITE(numout,*) '    saturation coefficient for O2    a_0 =', a_0
         WRITE(numout,*) '    saturation coefficient for O2    a_1 =', a_1
         WRITE(numout,*) '    saturation coefficient for O2    a_2 =', a_2
         WRITE(numout,*) '    saturation coefficient for O2    a_3 =', a_3
         WRITE(numout,*) '    saturation coefficient for O2    a_4 =', a_4
         WRITE(numout,*) '    saturation coefficient for O2    a_5 =', a_5
         WRITE(numout,*) '    saturation coefficient for O2    b_0 =', b_0
         WRITE(numout,*) '    saturation coefficient for O2    b_1 =', b_1
         WRITE(numout,*) '    saturation coefficient for O2    b_2 =', b_2
         WRITE(numout,*) '    saturation coefficient for O2    b_3 =', b_3
         WRITE(numout,*) '    saturation coefficient for O2    c_0 =', c_0
         WRITE(numout,*) '       Schmidt coefficient for O2 a_1_o2 =', a_1_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_2_o2 =', a_2_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_3_o2 =', a_3_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_4_o2 =', a_4_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_5_o2 =', a_5_o2
         WRITE(numout,*) '     Schmidt coefficient for CO2 a_1_co2 =', a_1_co2
         WRITE(numout,*) '     Schmidt coefficient for CO2 a_2_co2 =', a_2_co2
         WRITE(numout,*) '     Schmidt coefficient for CO2 a_3_co2 =', a_3_co2
         WRITE(numout,*) '     Schmidt coefficient for CO2 a_4_co2 =', a_4_co2
         WRITE(numout,*) '     Schmidt coefficient for CO2 a_5_co2 =', a_5_co2
      ENDIF

      ! namblingiron : Iron cycle parameters
      REWIND (numnatg)
      READ   (numnatg, namblingiron)

      kfe_inorg = kfe_inorg/86400.d0
      kfe_org   = kfe_org/86400.d0

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingiron'
         WRITE(numout,*) '    Constant for iron binding with ligands  kfe_eq_lig_max   = ', kfe_eq_lig_max
         WRITE(numout,*) '    Minimum ligand strength                 kfe_eq_lig_min   = ', kfe_eq_lig_min
         WRITE(numout,*) '    Irradiance scaling cte for Kfe          kfe_eq_lig_irr   = ', kfe_eq_lig_irr
         WRITE(numout,*) '    Low-iron threshold for lig stability    kfe_eq_lig_femin = ', kfe_eq_lig_femin
         WRITE(numout,*) '    Global uniform iron ligand concentration     felig_bkg   = ', felig_bkg
         WRITE(numout,*) '    1.5-order iron scavenging                    kfe_inorg   = ', kfe_inorg
         WRITE(numout,*) '    Adsorption rate coefficient for detritus       kfe_org   = ', kfe_org
         WRITE(numout,*) '    Minimum [o2] for oxic remineralization         oxy_min   = ', oxy_min        
         WRITE(numout,*) '    Prevent oxygen from becoming negative    ln_prev_o2lt0   = ', ln_prev_o2lt0
      ENDIF

      ! namblingbudget : Spatially integrated values for tracer budget
      REWIND (numnatg)
      READ   (numnatg, namblingbudget)

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist parameter for mass conservation checking BLING'
         WRITE(numout,*) '    Flag to check mass conservation of phosphate = ', ln_bling_mass
         WRITE(numout,*) '    Flag to check external fluxes of tracers     = ', ln_bling_ext
      ENDIF

      !                
   END  SUBROUTINE  trc_nam_bling
   
   !!======================================================================
END MODULE trcnam_bling
