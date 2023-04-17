MODULE trcext_blingv0

   !!======================================================================
   !!                     ***  MODULE trcext_bling  ***
   !! TOP :  Iron external inputs
   !!======================================================================
   !! History :  MClaret@McGill@04/2016. CO2 external fluxes added
   !! History :  MClaret@McGill@04-07/2014
   !!======================================================================

   USE oce_trc
   USE par_trc
   USE trc
   USE fldread 
   USE iom
   USE dom_oce

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_ext_bling
   PUBLIC trc_ext_init_bling

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_patm_bling ! structure of input fields (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_dust_bling

#include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_ext_bling ( kt, Kbb, Kmm, Krhs ) 

      !-----------------------------------------------------------------------
      ! The prefixes "b" refers to a bottom flux and "s" to a sea surface flux
      !-----------------------------------------------------------------------

      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level index

      INTEGER  :: ji, jj, jk, ikb, nkk

      REAL(wp) :: zrfact
      REAL(wp) :: sss, sst, o2fact
      REAL(wp) :: tt, tk, ts1, ts2, ts3, ts4, ts5

      REAL(wp) :: ztc, ztc2, ztc3, ztc4, zws, zkgwan, xconv
      REAL(wp) :: zsch_co2, zsch_o2
      REAL(wp) :: co2_sat, co2_sur, tmp1, tmp2
      REAL(wp) :: o2_alpha, o2_alpha2, o2_sat, o2_sur, pco2_atm

      REAL(wp) :: foxy, fe_2_p_sed
      REAL(wp) :: sumbpo4_glob, sumbfed_glob, sumdust_glob, sumboxy_glob,sumo2flx_glob
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: tmask_ikb
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: sch_no_term_co2, co2_flx
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: sch_no_term_o2 , o2_flx
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: bpo4, bfed, boxy, bdic, balk
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: bno3

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ext_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ALLOCATE ( tmask_ikb(jpi,jpj) )
      ALLOCATE ( sch_no_term_co2(jpi,jpj) )
      ALLOCATE ( co2_flx(jpi,jpj) )
      ALLOCATE ( sch_no_term_o2(jpi,jpj) )
      ALLOCATE ( o2_flx(jpi,jpj) )
      ALLOCATE ( bpo4(jpi,jpj) )
      ALLOCATE ( bfed(jpi,jpj) )
      ALLOCATE ( boxy(jpi,jpj) )
      ALLOCATE ( bdic(jpi,jpj) )
      ALLOCATE ( balk(jpi,jpj) )
      IF ( ln_nitro ) THEN
              ALLOCATE ( bno3(jpi,jpj) )
      ENDIF

      !---------------------------------------------------------------------
      ! Calculate air-sea exhange for O2/CO2
      !---------------------------------------------------------------------

      ! Get atmospheric pressure if read from file
      IF( ln_patm_bling ) THEN
         CALL fld_read( kt, 1, sf_patm_bling )
         if (ln_Pa2atm) then  ! Apr 25, 2016 (xhu): convert Pa to atm, xhu
            patm_bling(:,:) = 9.86923e-6_wp*sf_patm_bling(1)%fnow(:,:,1)
         else
            patm_bling(:,:) = sf_patm_bling(1)%fnow(:,:,1)
         endif
      ENDIF

      ! Get pco2 atmospheric value
      IF ( (nyear .GT. 1958) .AND. (nyear .LT. 2101) ) THEN
        pco2_atm=pco2(nyear-1958)
      ELSE
        pco2_atm=pco2(1)
      ENDIF
#if defined key_xhu
      if (narea==0) write(*,*) 'trcext_bling: pco2_atm=',pco2_atm,' ln_Pa2atm=',ln_Pa2atm
#endif

      xconv  = 0.01_wp / 3600._wp ! convert from cm/h to m/s (piston velocity)
      o2fact = 1.d0 /22.3916d0    ! convert from ml/L to mol/m3, old version BLING: o2fact = 1027.d0 /22391.6d0
      ! Compute gas exchange coefficients (wind and ice effects considered)
      ! Ported exactly from PISCES scheme p4zflx
!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztc  = MIN( 35., ts(ji,jj,1,jp_tem,Kmm) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ztc4 = ztc * ztc3 ! add term for Shimdtt number calculatio following Wanninkhof 2014 (LCastro 2018)

            ! Compute the schmidt Number both O2 and CO2 
!            zsch_co2 = 2073.1 - 125.62 * ztc + 3.6276 * ztc2 - 0.043126 * ztc3
!            zsch_o2  = 1953.4 - 128.0  * ztc + 3.9918 * ztc2 - 0.050091 * ztc3

            ! Compute the schmidt Number both O2 and CO2 using Wanninkhof 2014 and variables in namelist (LCastro 2018) 
             zsch_co2 = a_1_co2 + a_2_co2 * ztc + a_3_co2 * ztc2 + a_4_co2 * ztc3 + a_5_co2 * ztc4
             zsch_o2  = a_1_o2  + a_2_o2  * ztc + a_3_o2  * ztc2 + a_4_o2  * ztc3 + a_5_o2  * ztc4



            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
           ! zkgwan = 0.3 * zws  + 2.5 * ( 0.5246 + 0.016256 * ztc + 0.00049946  * ztc2 ) ! coeficient = 0.3 (cm/h)(m/s)^-2  ! Wanninkhof 1992
           !zkgwan = zkgwan * xconv * ( 1.d0 - fr_i(ji,jj) ) * tmask(ji,jj,1)
           
           ! adjust coeficient to 0.251(cm/h)(m/s)^-2 follwing Wanninkhof 2014 (LCastro 2018) 
           !"air/sea gas exchange follows OMIP6 protocols (Orr et al 2016) based on Waninkhof 2014"
           zkgwan = 0.251 * zws * xconv * ( 1.d0 - fr_i(ji,jj) ) * tmask(ji,jj,1)

            ! air-sea gas transfer velocity (piston velocity). Units are m/s
            sch_no_term_co2(ji,jj) = zkgwan * SQRT( 660.d0 / zsch_co2 )
            sch_no_term_o2 (ji,jj) = zkgwan * SQRT( 660.d0 / zsch_o2 )

         END DO
      END DO

#if defined key_xhu
      write(*,*) 'narea=',narea,' max sch_no_term_co2:',maxval(sch_no_term_co2(:,:)),' sch_no_term_o2=',maxval(sch_no_term_o2), &
                                      ' co2_csurf=',maxval(co2_csurf),' patm_bling=',maxval(patm_bling), ' co2_alpha=',maxval(co2_alpha)
#endif
      ! Gas air-sea fluxes
      DO jj=1,jpj
        DO ji=1,jpi

          zrfact = rfact / e3t(ji,jj,1,Kmm) ! [s/m]

          !-----------------
          ! CO2 air-sea flux
          ! saturation concentration (solubility from co2calc subroutine) * piston v [mol/m3 m/s]
          co2_sat=patm_bling(ji,jj)*pco2_atm*1.e-6*co2_alpha(ji,jj)*sch_no_term_co2(ji,jj) 

          ! surface concentration (computed in co2calc subroutine) * piston v [mol/m3 m/s]
          co2_sur=co2_csurf(ji,jj)*sch_no_term_co2(ji,jj)

          ! air-sea flux [mol/m3 m/s]
          co2_flx(ji,jj)=co2_sat-co2_sur

          ! add trend
          tr(ji,jj,1,jpDIC_bling,Krhs)=tr(ji,jj,1,jpDIC_bling,Krhs)+co2_flx(ji,jj)*zrfact

          !----------------
          ! O2 air-sea flux
          ! Compute O2 solubility (GARCIA and GORDON 1992, eq.8)
          !sst=MIN( 35.,tsn(ji,jj,1,jp_tem) )
          sst=ts(ji,jj,1,jp_tem,Kmm)
          sss=ts(ji,jj,1,jp_sal,Kmm)

          tt=298.15d0-sst
          tk=273.15d0+sst
          ts1=LOG(tt/tk)
          ts2=ts1 *ts1
          ts3=ts2*ts1
          ts4=ts3*ts1
          ts5=ts4*ts1

          ! [mol/m3]
          tmp1= a_0 + a_1*ts1 + a_2*ts2 + a_3*ts3 + a_4*ts4 + a_5*ts5
          tmp2= sss * (b_0 + b_1*ts1 + b_2*ts2 + b_3*ts3) + c_0*sss*sss
          o2_alpha=o2fact*EXP(tmp1+tmp2)

          ! [mol/m3/atm]
          o2_alpha2=o2_alpha/po2_atm

          ! saturation concentration * piston v [mol/m2/s]
          !o2_sat=patm_bling(ji,jj)*po2_atm*o2_alpha2*sch_no_term_o2(ji,jj)
          o2_sat=po2_atm*o2_alpha2*sch_no_term_o2(ji,jj)

          ! surface concentration * piston v [mol/m2/s]
          o2_sur=tr(ji,jj,1,jpOxy_bling,Kbb)*sch_no_term_o2(ji,jj)


          ! Flux of oxygen (mol/m2/s)
          o2_flx(ji,jj)=o2_sat-o2_sur

          tr(ji,jj,1,jpOxy_bling,Krhs)=tr(ji,jj,1,jpOxy_bling,Krhs)+o2_flx(ji,jj)*zrfact

          !sumtfoxy=sumtfoxy+o2flx(ji,jj)*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)

        ENDDO
      ENDDO
#if defined key_xhu
      write(*,*) 'narea=',narea,' max o2flx=',maxval(o2_flx),' co2_flx=',maxval(co2_flx)
#endif

      !---------------------------------------------------------------------

#if defined key_rnf_nutrients
     ! Nutrient in river inflow
     ! Add nutrients in the river inflow (Laura Castro Oct 2016)
     ! -------------------------------------
     !rdttrc(1) model time step in seconds
     !nkk=nkrnf ! ! nkrnf is the k level for enhanced mixing @ the river mouth
     zrfact=rfact/1027.0_wp !density 1027kg/m3 consistant with other uses of density by bling
     !  deltaC = rnf_c * rnf/h_rnf * rdt / rho
     !         = mol/m^3 * kg/m^2/s * 1/m * s * m^3/kg = mol/m^3
     tr(:,:,1,jpPO4_bling,Krhs) = tr(:,:,1,jpPO4_bling,Krhs) + rnfpo4(:,:)*rnf(:,:)/e3t(:,:,1,Kmm)*zrfact
     tr(:,:,1,jpFed_bling,Krhs) = tr(:,:,1,jpFed_bling,Krhs) + rnffed(:,:)*rnf(:,:)/e3t(:,:,1,Kmm)*zrfact
     tr(:,:,1,jpDOP_bling,Krhs) = tr(:,:,1,jpDOP_bling,Krhs) + rnfdop(:,:)*rnf(:,:)/e3t(:,:,1,Kmm)*zrfact
     tr(:,:,1,jpDIC_bling,Krhs) = tr(:,:,1,jpDIC_bling,Krhs) + rnfdic(:,:)*rnf(:,:)/e3t(:,:,1,Kmm)*zrfact
     tr(:,:,1,jpalk_bling,Krhs) = tr(:,:,1,jpalk_bling,Krhs) + rnfalk(:,:)*rnf(:,:)/e3t(:,:,1,Kmm)*zrfact     
    !nkk=1
    !DO jk = 1, nkk  
    !   DO jj = 1, jpj
    !      DO ji = 1, jpi
    !          trn(ji,jj,jk,jpPO4_bling) = trn(ji,jj,jk,jpPO4_bling) + rnfpo4(ji,jj)*rnf(ji,jj)/fse3t(ji,jj,jk)*zrfact
    !          trn(ji,jj,jk,jpFed_bling) = trn(ji,jj,jk,jpFed_bling) + rnffed(ji,jj)*rnf(ji,jj)/fse3t(ji,jj,jk)*zrfact
    !          trn(ji,jj,jk,jpDOP_bling) = trn(ji,jj,jk,jpDOP_bling) + rnfdop(ji,jj)*rnf(ji,jj)/fse3t(ji,jj,jk)*zrfact
    !      END DO
    !   END DO
    !END DO
#endif
      !---------------------------------------------------------------------
      ! Calculate external fluxes for iron. 
      !---------------------------------------------------------------------

      !Get dust field
      IF( ln_dust_bling ) THEN
         CALL fld_read( kt, 1, sf_dust_bling )
         dust_bling(:,:) = sf_dust_bling(1)%fnow(:,:,1)
      ENDIF

      fe_2_p_sed=106.0e-4

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! [s/m] -> dust deposition in mol Fe/m3
            zrfact  = rfact / e3t(ji,jj,1,Kmm)

            ! [mol Fe/m3]
            tr(ji,jj,1,jpFed_bling,Krhs) =  tr(ji,jj,1,jpFed_bling,Krhs) + dust_bling(ji,jj)*zrfact

            !sumtffed=sumtffed+dust_bling(ji,jj)*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)

         ENDDO
      ENDDO

      !---------------------------------------------------------------------

      ! Exchange with sediments
      ! -------------------------------------
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor. For P, the bottom flux  
      ! puts the sinking flux reaching the bottom cell into the water column 
      ! through diffusion. For iron, the sinking flux disappears into the 
      ! sediments if bottom waters are oxic (assumed adsorbed as oxides),
      ! while an efflux of dissolved iron occurs dependent on the supply of
      ! reducing organic matter (scaled by the org-P sedimentation rate).
      ! If bottom waters are anoxic, the sinking flux of Fe is returned to
      ! the water column. Note this is not appropriate for very long runs
      ! with an anoxic ocean (iron will keep accumulating forever).
      ! For oxygen, the consumption of oxidant required to respire  
      ! the settling flux of organic matter (in support of the
      ! PO4 bottom flux) diffuses from the bottom water into the sediment.
      ! Do not bury any C - all goes back to water column (for DIC).
      ! Return all CaCO3 to the water column to preserve the alkalinity inventory
      ! Note that the flux of NO3 out of the sediment, inferred from 
      ! the PO4 flux, causes a negative flux of alkalinity.
      ! -------------------------------------

      DO jj = 1, jpj
        DO ji = 1, jpi

            ! mbkt is a matrix containing the vertical index of the
            ! bottom layer at each horizontal point
            ikb     = mbkt(ji,jj)
            tmask_ikb(ji,jj)=tmask(ji,jj,ikb)

            foxy   = tr(ji,jj,ikb,jpOxy_bling,Krhs)
            zrfact  = rfact / e3t(ji,jj,ikb,Kmm)

            ! Phosphate [mol P/m2/s]
            bpo4(ji,jj)=fpop_b(ji,jj)

            ! Oxygen [mol O2/m2/s]
            boxy(ji,jj) = -oxy2p*fpop_b(ji,jj)

            ! Iron [mol Fe/m2/s]
            IF (foxy>oxy_min) THEN
               bfed(ji,jj)= fe_2_p_sed*fpop_b(ji,jj)
            ELSE
               bfed(ji,jj)=(fe_2_p_sed*fpop_b(ji,jj)+fpofe_b(ji,jj))
            ENDIF

            ! DIC [mol C/m2/s]
            bdic(ji,jj)=fpop_b(ji,jj)*c2p+fcaco3_b(ji,jj)

            ! ALK [mol eq/m2/s]
            balk(ji,jj)=2.d0*fcaco3_b(ji,jj)-fpop_b(ji,jj)*n2p

            ! Add the bottom flux trend [mol/m3]
            tr(ji,jj,ikb,jpPO4_bling,Krhs) = tr(ji,jj,ikb,jpPO4_bling,Krhs) + bpo4(ji,jj)*zrfact
            tr(ji,jj,ikb,jpFed_bling,Krhs) = tr(ji,jj,ikb,jpFed_bling,Krhs) + bfed(ji,jj)*zrfact
            tr(ji,jj,ikb,jpOxy_bling,Krhs) = tr(ji,jj,ikb,jpOxy_bling,Krhs) + boxy(ji,jj)*zrfact
            tr(ji,jj,ikb,jpDIC_bling,Krhs) = tr(ji,jj,ikb,jpDIC_bling,Krhs) + bdic(ji,jj)*zrfact
            tr(ji,jj,ikb,jpalk_bling,Krhs) = tr(ji,jj,ikb,jpalk_bling,Krhs) + balk(ji,jj)*zrfact

            bpo4(ji,jj)=bpo4(ji,jj)*tmask(ji,jj,ikb)
            bfed(ji,jj)=bfed(ji,jj)*tmask(ji,jj,ikb)
            boxy(ji,jj)=boxy(ji,jj)*tmask(ji,jj,ikb)

            !!! --- AGT: Nitrogen fluxes in bottom cell --- !!!
            ! Choose here not to calculate denitrification rates (as in MITgcm), or burial rates
            ! Instead impose constant denitrification and sediment burial rates
            ! based on review paper "Nitrogen in the Baltic Sea..." (Lonborg & Markager, 2021)
            ! i) If the flux of particulate nitrogen is less than the combined loss rate 
            ! due to burial and denitrification, there is zero nitrate flux from sediments.
            ! ii) If the flux of particulate nitrogen excees the combined loss rate then
            ! the excess particulate flux forms the nitrate flux out of sediments.
            
            IF ( ln_nitro ) THEN
                    ! Nitrate [mol N/m2/s]
                    bno3(ji,jj) = fpon_b(ji,jj) - fbur - fden 
                    bno3(ji,jj) = MAX(bno3(ji,jj),0.e0)
                    ! Also need to reduce oxygen consumption according to denitrification
                    boxy(ji,jj) = boxy(ji,jj) + 1.25d0*fden
                    boxy(ji,jj) = MIN(boxy(ji,jj),0.e0)

                    tr(ji,jj,ikb,jpNO3_bling,Krhs) = tr(ji,jj,ikb,jpNO3_bling,Krhs) + bno3(ji,jj)*zrfact
                    bno3(ji,jj)=bno3(ji,jj)*tmask(ji,jj,ikb)
            ENDIF

         END DO
      END DO

      IF (ln_bling_ext) THEN

          IF( kt == nittrc000 ) THEN
            if (lwp) then 
               ! xhu
               CALL ctl_opn( numphpext, 'php.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
               CALL ctl_opn( numfedext, 'fed.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
               CALL ctl_opn( numoxyext, 'oxy.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
           endif
          ENDIF ! nittrc000

          ! xhu
          sumbpo4_glob= glob_sum( 'trcext_blingv0',bpo4(:,:)*e1e2t(:,:) )
          sumbfed_glob= glob_sum( 'trcext_blingv0',bfed(:,:)*e1e2t(:,:) )
          sumdust_glob= glob_sum( 'trcext_blingv0',dust_bling(:,:)*e1e2t(:,:)*tmask(:,:,1) )
          sumboxy_glob= glob_sum( 'trcext_blingv0',boxy(:,:)*e1e2t(:,:) )
          sumo2flx_glob= glob_sum( 'trcext_blingv0',o2_flx(:,:)*e1e2t(:,:)*tmask(:,:,1) )
          if (lwp) then
             WRITE(numphpext,9500) kt, sumbpo4_glob 
             WRITE(numfedext,9501) kt, sumbfed_glob, sumdust_glob
             WRITE(numoxyext,9501) kt, sumboxy_glob, sumo2flx_glob 
          endif ! lwp

9500  FORMAT(i10,e18.10)
9501  FORMAT(i10,2(e18.10))

      ENDIF ! ext

!      IF (lk_iomput) THEN
!        CALL iom_put( "po4_btf", bpo4      (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "oxy_stf", o2_flx    (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "oxy_btf", boxy      (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "fed_stf", dust_bling(:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "fed_btf", bfed      (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "fed_bur", fpofe_b   (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "dic_stf", co2_flx   (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "dic_btf", bdic      (:,:)*tmask_ikb(:,:) )
!        CALL iom_put( "alk_btf", balk      (:,:)*tmask_ikb(:,:) )
!      ENDIF

      IF (  .NOT. lk_iomput ) THEN
        trc2d( :,:,jp_bling0_2d+1 ) = bpo4      (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+2 ) = o2_flx    (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+3 ) = boxy      (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+4 ) = dust_bling(:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+5 ) = bfed      (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+6 ) = fpofe_b   (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+7 ) = co2_flx   (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+8 ) = bdic      (:,:)*tmask_ikb(:,:)
        trc2d( :,:,jp_bling0_2d+9 ) = balk      (:,:)*tmask_ikb(:,:)
      ENDIF

      DEALLOCATE ( tmask_ikb )
      DEALLOCATE ( sch_no_term_co2 )
      DEALLOCATE ( co2_flx )
      DEALLOCATE ( sch_no_term_o2 )
      DEALLOCATE ( o2_flx )
      DEALLOCATE ( bpo4 )
      DEALLOCATE ( bfed )
      DEALLOCATE ( boxy )
      DEALLOCATE ( bdic )
      DEALLOCATE ( balk )
      IF ( ln_nitro ) THEN
              DEALLOCATE ( bno3 )
      ENDIF

   END SUBROUTINE trc_ext_bling

   SUBROUTINE trc_ext_init_bling

      INTEGER  :: jm, ierr
      INTEGER  :: numdust, ntimes_dust
      INTEGER  :: numpatm, ntimes_patm
#if defined key_rnf_nutrients
      INTEGER  :: numrnf_nutrient
#endif 

      REAL(wp), DIMENSION(12) :: zsteps                 ! times records

      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) :: zdust, zpatm
      !!----------------------------------------------------------------------

      ALLOCATE( sf_dust_bling(1), STAT=ierr )    
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_ext_init_bling: unable to allocate sf_apr structure' )

      CALL fld_fill( sf_dust_bling, (/ sn_dust_bling /), cn_dir_bling, 'trc_ext_init_bling', 'Iron from sediment ', 'namblingiron' )

      ALLOCATE( sf_dust_bling(1)%fnow(jpi,jpj,1) )
      IF( sn_dust_bling%ln_tint ) ALLOCATE( sf_dust_bling(1)%fdta(jpi,jpj,1,2) )

#if defined key_rnf_nutrients
      CALL iom_open ( 'bling_rnf_nutrients.nc',numrnf_nutrient)         ! open file
      CALL iom_get  ( numrnf_nutrient, jpdom_data,'rnf_po4',rnfpo4)    ! read the river mouth array
      CALL iom_get  ( numrnf_nutrient, jpdom_data,'rnf_dop',rnfdop)    ! read the river mouth array
      CALL iom_get  ( numrnf_nutrient, jpdom_data,'rnf_fed',rnffed)    ! read the river mouth array
      CALL iom_get  ( numrnf_nutrient, jpdom_data,'rnf_dic',rnfdic)    ! read the river mouth array
      CALL iom_get  ( numrnf_nutrient, jpdom_data,'rnf_alk',rnfalk)    ! read the river mouth array
      CALL iom_close( numrnf_nutrient )                                   ! close file
#endif

     !CALL iom_open (  TRIM( sn_dust_bling%clname ) , numdust )

     !CALL iom_gettime( numdust, zsteps, kntime=ntimes_dust)  ! get number of record in file

     !ALLOCATE( zdust(jpi,jpj,ntimes_dust) )
     !DO jm = 1, ntimes_dust
     !   CALL iom_get( numdust, jpdom_data, TRIM( sn_dust_bling%clvar ), zdust(:,:,jm), jm )
     !ENDDO
     ! 
     !CALL iom_close( numdust )
     !DEALLOCATE( zdust)

      !!-------------
      
      IF ( ln_patm_bling ) THEN

         ALLOCATE( sf_patm_bling(1), STAT=ierr )       
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_patm structure' )
         
         CALL fld_fill( sf_patm_bling, (/ sn_patm_bling /), cn_dir_patm_bling, 'trc_ext_init_bling', 'Atmospheric pressure ', 'namblingpatm' )
         
         ALLOCATE( sf_patm_bling(1)%fnow(jpi,jpj,1)   )
         IF( sn_patm_bling%ln_tint )  ALLOCATE( sf_patm_bling(1)%fdta(jpi,jpj,1,2) )
         
         !CALL iom_open (  TRIM( sn_patm_bling%clname ) , numpatm )
         
         !CALL iom_gettime( numpatm, zsteps, kntime=ntimes_patm)  ! get number of record in file
         
        !ALLOCATE( zpatm(jpi,jpj,ntimes_patm) )
        !
        !DO jm = 1, ntimes_patm
        !  CALL iom_get( numpatm, jpdom_data, TRIM( sn_patm_bling%clvar ), zpatm(:,:,jm), jm )
        !ENDDO
        
        !CALL iom_close( numpatm )
        !DEALLOCATE( zpatm)

      ELSE

        patm_bling(:,:)=1.e0;

      ENDIF

   END SUBROUTINE trc_ext_init_bling

END MODULE trcext_blingv0
numpatm )
        !DEALLOCATE( zpatm)

      ELSE

        patm_bling(:,:)=1.e0;

      ENDIF

   END SUBROUTINE trc_ext_init_bling

END MODULE trcext_blingv0
