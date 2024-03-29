MODULE par_bling
   !!======================================================================
   !!                        ***  par_bling  ***
   !! TOP :   set the BLINGv0 parameters
   !!======================================================================
   !! History : 2014-2016 BLING tracers added
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_bling.F90 -1   $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   USE par_trc         ! TOP parameters

!   USE par_lobster, ONLY : jp_lobster      !: number of tracers in LOBSTER
!   USE par_lobster, ONLY : jp_lobster_2d   !: number of 2D diag in LOBSTER
!   USE par_lobster, ONLY : jp_lobster_3d   !: number of 3D diag in LOBSTER
!   USE par_lobster, ONLY : jp_lobster_trd  !: number of biological diag in LOBSTER

!   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
!   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
!   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
!   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

!   USE par_cfc    , ONLY : jp_cfc          !: number of tracers in CFC
!   USE par_cfc    , ONLY : jp_cfc_2d       !: number of tracers in CFC
!   USE par_cfc    , ONLY : jp_cfc_3d       !: number of tracers in CFC
!   USE par_cfc    , ONLY : jp_cfc_trd      !: number of tracers in CFC

!   USE par_c14b   , ONLY : jp_c14b         !: number of tracers in C14
!   USE par_c14b   , ONLY : jp_c14b_2d      !: number of tracers in C14
!   USE par_c14b   , ONLY : jp_c14b_3d      !: number of tracers in C14
!   USE par_c14b   , ONLY : jp_c14b_trd     !: number of tracers in C14

   IMPLICIT NONE

!   INTEGER, PARAMETER ::   jp_lm      = jp_lobster     + jp_pisces     + jp_cfc     + jp_c14b     !: 
!   INTEGER, PARAMETER ::   jp_lm_2d   = jp_lobster_2d  + jp_pisces_2d  + jp_cfc_2d  + jp_c14b_2d  !:
!   INTEGER, PARAMETER ::   jp_lm_3d   = jp_lobster_3d  + jp_pisces_3d  + jp_cfc_3d  + jp_c14b_3d  !:
   !INTEGER, PARAMETER ::   jp_lm_trd  = jp_lobster_trd + jp_pisces_trd + jp_cfc_trd + jp_c14b_trd !:
!   INTEGER, PARAMETER ::   jp_lm_trd  = jp_lobster_trd + jp_cfc_trd + jp_c14b_trd !:

! AGT: Remove dependencies on other bgc modules.
   INTEGER, PARAMETER :: jp_lm     = 0
   INTEGER, PARAMETER :: jp_lm_2d  = 0
   INTEGER, PARAMETER :: jp_lm_3d  = 0
   INTEGER, PARAMETER :: jp_lm_trd = 0
! -----------------------------------------

   !!---------------------------------------------------------------------
   !!   'key_bling'                     user defined tracers (BLINGv0)
   !!---------------------------------------------------------------------

   LOGICAL, PUBLIC, PARAMETER ::   lk_bling     = .TRUE.   !: PTS flag 
! AGT: jp_bling now defined in par_trc
!   INTEGER, PUBLIC, PARAMETER ::   jp_bling     =  6       !: number of PTS tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling_2d  =  10      !: additional 2d output arrays. If not 1
   INTEGER, PUBLIC, PARAMETER ::   jp_bling_3d  =  36      !: additional 3d output arrays. If not 1
   INTEGER, PUBLIC, PARAMETER ::   jp_bling_trd =  1       !: number of sms trends for BLINGv0. If not 1

   ! assign an index in trc arrays for each PTS prognostic variables
   ! Note by MClaret@McGill@March16: First three tracers are preserved as PO4+DOP=cte. and DOP*106+DIC=cte.
   INTEGER, PUBLIC, PARAMETER ::   jpPO4_bling = jp_lm + 1     !: Dissolved phosphate
   INTEGER, PUBLIC, PARAMETER ::   jpDOP_bling = jp_lm + 2     !: Dissolved organic phosphate
   INTEGER, PUBLIC, PARAMETER ::   jpFed_bling = jp_lm + 3     !: Dissolved iron
   INTEGER, PUBLIC, PARAMETER ::   jpOxy_bling = jp_lm + 4     !: Dissolved oxygen
   INTEGER, PUBLIC, PARAMETER ::   jpDIC_bling = jp_lm + 5     !: Dissolved inorganic carbon
   INTEGER, PUBLIC, PARAMETER ::   jpalk_bling = jp_lm + 6     !: Total alkalinity concentration
   INTEGER, PUBLIC, PARAMETER ::   jpNO3_bling = jp_lm + 7     !: Dissolved nitrate
   INTEGER, PUBLIC, PARAMETER ::   jpDON_bling = jp_lm + 8     !: Dissolved organic nitrogen
   !DO_CARBON

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_bling0     = jp_lm     + 1              !: First index of BLINGv0 passive tracers
   ! AGT : Need to hardcode jp_blg1 = 8 if using nitrogen... 
   INTEGER, PUBLIC, PARAMETER ::  jp_blg0 = 1  !: First index of BLING tracers
   INTEGER, PUBLIC, PARAMETER ::  jp_blg1 = 8 !: Last  index of BLING tracers
!  INTEGER, PUBLIC, PARAMETER ::   jp_bling1     = jp_lm     + 6              !: Last  index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling0_2d  = jp_lm_2d  + 1              !: First index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling1_2d  = jp_lm_2d  + jp_bling_2d    !: Last  index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling0_3d  = jp_lm_3d  + 1              !: First index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling1_3d  = jp_lm_3d  + jp_bling_3d    !: Last  index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling0_trd = jp_lm_trd + 1              !: First index of BLINGv0 passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bling1_trd = jp_lm_trd + jp_bling_trd   !: Last  index of BLINGv0 passive tracers

   !!======================================================================
END MODULE par_bling
