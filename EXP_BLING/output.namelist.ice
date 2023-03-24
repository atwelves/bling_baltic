 &NAMPAR
 JPL     =           5,
 NLAY_I  =           2,
 NLAY_S  =           1,
 LN_VIRTUAL_ITD  = F,
 LN_ICEDYN       = T,
 LN_ICETHD       = T,
 RN_AMAX_N       =  0.997000000000000     ,
 RN_AMAX_S       =  0.997000000000000     ,
 CN_ICERST_IN    = restart_ice                                                                                                       
                                                                                                                                     
           ,
 CN_ICERST_INDIR = .                                                                                                                 
                                                                                                                                     
           ,
 CN_ICERST_OUT   = restart_ice                                                                                                       
                                                                                                                                     
           ,
 CN_ICERST_OUTDIR        = .                                                                                                         
                                                                                                                                     
                   
 /
 &NAMITD
 LN_CAT_HFN      = T,
 RN_HIMEAN       =   2.00000000000000     ,
 LN_CAT_USR      = F,
 RN_CATBND       =  0.000000000000000E+000,  0.450000000000000     ,   1.10000000000000     ,   2.10000000000000     ,
    3.70000000000000     ,   6.00000000000000     , 95*0.000000000000000E+000  ,
 RN_HIMIN        =  0.100000000000000     ,
 RN_HIMAX        =   99.0000000000000     
 /
 &NAMTHD
 LN_ICEDH        = T,
 LN_ICEDA        = T,
 LN_ICEDO        = T,
 LN_ICEDS        = T,
 LN_LEADHFX      = T
 /
 &NAMTHD_ZDF
 LN_ZDF_BL99     = T,
 LN_CNDI_U64     = F,
 LN_CNDI_P07     = T,
 RN_CND_S        =  0.310000000000000     ,
 RN_KAPPA_I      =   1.00000000000000     ,
 RN_KAPPA_S      =   10.0000000000000     ,
 RN_KAPPA_SMLT   =   7.00000000000000     ,
 RN_KAPPA_SDRY   =   10.0000000000000     ,
 LN_ZDF_CHKCVG   = F
 /
 &NAMTHD_DA
 RN_BETA =   1.00000000000000     ,
 RN_DMIN =   8.00000000000000     
 /
 &NAMTHD_DO
 RN_HINEW        =  0.100000000000000     ,
 LN_FRAZIL       = F,
 RN_MAXFRAZ      =   1.00000000000000     ,
 RN_VFRAZ        =  0.417000000000000     ,
 RN_CFRAZ        =   5.00000000000000     
 /
 &NAMTHD_SAL
 NN_ICESAL       =           2,
 RN_ICESAL       =   4.00000000000000     ,
 RN_SAL_GD       =   5.00000000000000     ,
 RN_TIME_GD      =   1730000.00000000     ,
 RN_SAL_FL       =   2.00000000000000     ,
 RN_TIME_FL      =   864000.000000000     ,
 RN_SIMAX        =   20.0000000000000     ,
 RN_SIMIN        =  0.100000000000000     
 /
 &NAMTHD_PND
 LN_PND  = T,
 LN_PND_LEV      = T,
 RN_APND_MIN     =  0.150000000000000     ,
 RN_APND_MAX     =  0.850000000000000     ,
 RN_PND_FLUSH    =  0.100000000000000     ,
 LN_PND_CST      = F,
 RN_APND =  0.200000000000000     ,
 RN_HPND =  5.000000000000000E-002,
 LN_PND_TOPO     = F,
 LN_PND_LIDS     = T,
 LN_PND_ALB      = T
 /
 &NAMSBC
 RN_CIO  =  5.000000000000000E-003,
 NN_SNWFRA       =           2,
 RN_SNWBLOW      =  0.660000000000000     ,
 NN_FLXDIST      =          -1,
 LN_CNDFLX       = F,
 LN_CNDEMULATE   = F,
 NN_QTRICE       =           0
 /
 &NAMINI
 LN_ICEINI       = F,
 NN_ICEINI_FILE  =           0,
 RN_THRES_SST    =   2.00000000000000     ,
 RN_HTI_INI_N    =   3.00000000000000     ,
 RN_HTI_INI_S    =   1.00000000000000     ,
 RN_HTS_INI_N    =  0.300000000000000     ,
 RN_HTS_INI_S    =  0.300000000000000     ,
 RN_ATI_INI_N    =  0.900000000000000     ,
 RN_ATI_INI_S    =  0.900000000000000     ,
 RN_SMI_INI_N    =   6.30000000000000     ,
 RN_SMI_INI_S    =   6.30000000000000     ,
 RN_TMI_INI_N    =   270.000000000000     ,
 RN_TMI_INI_S    =   270.000000000000     ,
 RN_TSU_INI_N    =   270.000000000000     ,
 RN_TSU_INI_S    =   270.000000000000     ,
 RN_TMS_INI_N    =   270.000000000000     ,
 RN_TMS_INI_S    =   270.000000000000     ,
 RN_APD_INI_N    =  0.200000000000000     ,
 RN_APD_INI_S    =  0.200000000000000     ,
 RN_HPD_INI_N    =  5.000000000000000E-002,
 RN_HPD_INI_S    =  5.000000000000000E-002,
 RN_HLD_INI_N    =  0.000000000000000E+000,
 RN_HLD_INI_S    =  0.000000000000000E+000,
 SN_HTI%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_HTI%FREQH   =  -12.0000000000000     ,
 SN_HTI%CLVAR   = hti                               ,
 SN_HTI%LN_TINT = F,
 SN_HTI%LN_CLIM = T,
 SN_HTI%CLFTYP  = yearly  ,
 SN_HTI%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_HTI%VCOMP   =                                   ,
 SN_HTI%LNAME   =                                   ,
 SN_HTS%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_HTS%FREQH   =  -12.0000000000000     ,
 SN_HTS%CLVAR   = hts                               ,
 SN_HTS%LN_TINT = F,
 SN_HTS%LN_CLIM = T,
 SN_HTS%CLFTYP  = yearly  ,
 SN_HTS%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_HTS%VCOMP   =                                   ,
 SN_HTS%LNAME   =                                   ,
 SN_ATI%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_ATI%FREQH   =  -12.0000000000000     ,
 SN_ATI%CLVAR   = ati                               ,
 SN_ATI%LN_TINT = F,
 SN_ATI%LN_CLIM = T,
 SN_ATI%CLFTYP  = yearly  ,
 SN_ATI%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_ATI%VCOMP   =                                   ,
 SN_ATI%LNAME   =                                   ,
 SN_TSU%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_TSU%FREQH   =  -12.0000000000000     ,
 SN_TSU%CLVAR   = tsu                               ,
 SN_TSU%LN_TINT = F,
 SN_TSU%LN_CLIM = T,
 SN_TSU%CLFTYP  = yearly  ,
 SN_TSU%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_TSU%VCOMP   =                                   ,
 SN_TSU%LNAME   =                                   ,
 SN_TMI%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_TMI%FREQH   =  -12.0000000000000     ,
 SN_TMI%CLVAR   = tmi                               ,
 SN_TMI%LN_TINT = F,
 SN_TMI%LN_CLIM = T,
 SN_TMI%CLFTYP  = yearly  ,
 SN_TMI%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_TMI%VCOMP   =                                   ,
 SN_TMI%LNAME   =                                   ,
 SN_SMI%CLNAME  = Ice_initialization                                                                                                 
                                                                                                                                     
          ,
 SN_SMI%FREQH   =  -12.0000000000000     ,
 SN_SMI%CLVAR   = smi                               ,
 SN_SMI%LN_TINT = F,
 SN_SMI%LN_CLIM = T,
 SN_SMI%CLFTYP  = yearly  ,
 SN_SMI%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_SMI%VCOMP   =                                   ,
 SN_SMI%LNAME   =                                   ,
 SN_TMS%CLNAME  = NOT USED                                                                                                           
                                                                                                                                     
          ,
 SN_TMS%FREQH   =  -12.0000000000000     ,
 SN_TMS%CLVAR   = tms                               ,
 SN_TMS%LN_TINT = F,
 SN_TMS%LN_CLIM = T,
 SN_TMS%CLFTYP  = yearly  ,
 SN_TMS%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_TMS%VCOMP   =                                   ,
 SN_TMS%LNAME   =                                   ,
 SN_APD%CLNAME  = NOT USED                                                                                                           
                                                                                                                                     
          ,
 SN_APD%FREQH   =  -12.0000000000000     ,
 SN_APD%CLVAR   = apd                               ,
 SN_APD%LN_TINT = F,
 SN_APD%LN_CLIM = T,
 SN_APD%CLFTYP  = yearly  ,
 SN_APD%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_APD%VCOMP   =                                   ,
 SN_APD%LNAME   =                                   ,
 SN_HPD%CLNAME  = NOT USED                                                                                                           
                                                                                                                                     
          ,
 SN_HPD%FREQH   =  -12.0000000000000     ,
 SN_HPD%CLVAR   = hpd                               ,
 SN_HPD%LN_TINT = F,
 SN_HPD%LN_CLIM = T,
 SN_HPD%CLFTYP  = yearly  ,
 SN_HPD%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_HPD%VCOMP   =                                   ,
 SN_HPD%LNAME   =                                   ,
 SN_HLD%CLNAME  = NOT USED                                                                                                           
                                                                                                                                     
          ,
 SN_HLD%FREQH   =  -12.0000000000000     ,
 SN_HLD%CLVAR   = hld                               ,
 SN_HLD%LN_TINT = F,
 SN_HLD%LN_CLIM = T,
 SN_HLD%CLFTYP  = yearly  ,
 SN_HLD%WNAME   =                                                                                                                    
                                                                                                                                     
          ,
 SN_HLD%VCOMP   =                                   ,
 SN_HLD%LNAME   =                                   ,
 CN_DIR  = ./                                                                                                                        
                                                                                                                                     
   
 /
 &NAMDYN
 LN_DYNALL       = T,
 LN_DYNRHGADV    = F,
 LN_DYNADV1D     = F,
 LN_DYNADV2D     = F,
 RN_UICE =  0.500000000000000     ,
 RN_VICE =  0.500000000000000     ,
 RN_ISHLAT       =   2.00000000000000     ,
 LN_LANDFAST_L16 = F,
 RN_LF_DEPFRA    =  0.125000000000000     ,
 RN_LF_BFR       =   15.0000000000000     ,
 RN_LF_RELAX     =  1.000000000000000E-005,
 RN_LF_TENSILE   =  5.000000000000000E-002,
 SN_ICBMSK%CLNAME  = NOT USED                                                                                                        
                                                                                                                                     
             ,
 SN_ICBMSK%FREQH   =  -12.0000000000000     ,
 SN_ICBMSK%CLVAR   = icb_mask                          ,
 SN_ICBMSK%LN_TINT = F,
 SN_ICBMSK%LN_CLIM = T,
 SN_ICBMSK%CLFTYP  = yearly  ,
 SN_ICBMSK%WNAME   =                                                                                                                 
                                                                                                                                     
             ,
 SN_ICBMSK%VCOMP   =                                   ,
 SN_ICBMSK%LNAME   =                                   ,
 CN_DIR  = ./                                                                                                                        
                                                                                                                                     
   
 /
 &NAMDYN_RDGRFT
 LN_STR_H79      = T,
 RN_PSTAR        =   20000.0000000000     ,
 RN_CRHG =   20.0000000000000     ,
 LN_STR_R75      = F,
 RN_PE_RDG       =   17.0000000000000     ,
 LN_STR_CST      = F,
 RN_STR  =  0.000000000000000E+000,
 LN_STR_SMOOTH   = T,
 LN_DISTF_LIN    = T,
 LN_DISTF_EXP    = F,
 RN_MURDG        =   3.00000000000000     ,
 RN_CSRDG        =  0.500000000000000     ,
 LN_PARTF_LIN    = F,
 RN_GSTAR        =  0.150000000000000     ,
 LN_PARTF_EXP    = T,
 RN_ASTAR        =  3.000000000000000E-002,
 LN_RIDGING      = T,
 RN_HSTAR        =   25.0000000000000     ,
 RN_PORORDG      =  0.300000000000000     ,
 RN_FSNWRDG      =  0.500000000000000     ,
 RN_FPNDRDG      =   1.00000000000000     ,
 LN_RAFTING      = T,
 RN_HRAFT        =  0.750000000000000     ,
 RN_CRAFT        =   5.00000000000000     ,
 RN_FSNWRFT      =  0.500000000000000     ,
 RN_FPNDRFT      =   1.00000000000000     
 /
 &NAMDYN_RHG
 LN_RHG_EVP      = T,
 LN_AEVP = T,
 LN_RHG_EAP      = F,
 RN_CREEPL       =  2.000000000000000E-009,
 RN_ECC  =   2.00000000000000     ,
 NN_NEVP =         100,
 RN_RELAST       =  0.333000000000000     ,
 NN_RHG_CHKCVG   =           0,
 LN_RHG_VP       = F,
 NN_VP_NOUT      =          10,
 NN_VP_NINN      =        1500,
 NN_VP_CHKCVG    =           5
 /
 &NAMDYN_ADV
 LN_ADV_PRA      = T,
 LN_ADV_UMX      = F,
 NN_UMX  =           5
 /
 &NAMALB
 RN_ALB_SDRY     =  0.850000000000000     ,
 RN_ALB_SMLT     =  0.750000000000000     ,
 RN_ALB_IDRY     =  0.600000000000000     ,
 RN_ALB_IMLT     =  0.500000000000000     ,
 RN_ALB_DPND     =  0.270000000000000     ,
 RN_ALB_HPIV     =   1.00000000000000     
 /
 &NAMDIA
 LN_ICEDIACHK    = F,
 RN_ICECHK_CEL   =   1.00000000000000     ,
 RN_ICECHK_GLO   =  1.000000000000000E-004,
 LN_ICEDIAHSB    = F,
 LN_ICECTL       = F,
 IICEPRT =          10,
 JICEPRT =          10
 /
