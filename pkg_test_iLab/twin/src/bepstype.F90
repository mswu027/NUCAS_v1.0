

module bepstype
  use shr_kind_mod,only:r8=>shr_kind_r8
  !--iLab::no need for use-statement
  ! use beps_par
  implicit none

  ! For restart
  !-- iLab::module field should have save attribute
  real(r8), allocatable,public,save     ::  v2last(:,:,:)


  !*********************Climate Forcing datasets**********************
  type,public:: forc
     real(r8),pointer::  Temp(:)
     real(r8),pointer::  Tempmx(:)
     real(r8),pointer::  Tempmn(:)
     real(r8),pointer::  Wind(:)
     real(r8),pointer::  Srad(:)   !downward solar raditation
     real(r8),pointer::  Rh(:)     !relative humidity(%)
     real(r8),pointer::  Rain(:)   ! liquid precipitation
     real(r8),pointer::  Snow(:)   ! Snow rate
     real(r8),pointer::  Swdr(:)   ! SW direct radiation
     real(r8),pointer::  Swdf(:)   ! SW diffuse radiation
     character(len=16) :: meteo_ref_yyyymmdd !--iLab::to consistenly handle tempooral settings
                                             !--      expected format yyyy-mm-dd
  end type forc

  type(forc),save,target,public:: clim

  !***********************for CPL datasets****************************

  !******************for boundary/yrdata/lai/cpools/
  type,public::surf
     integer,pointer::  lcno(:,:)     ! PFT types
     integer,pointer::  stext(:)      ! soil texture
     real(r8),pointer:: PCT_PFT(:,:) ! PFT fraction
     real(r8),pointer:: clumping(:)
     real(r8),pointer:: longitude(:)
     real(r8),pointer:: latitude(:)

     real(r8),pointer:: sdp(:)       ! snowdepth
     real(r8),pointer:: st(:)        !
     real(r8),pointer:: sw(:)        ! soil moisture

     real(r8),pointer:: laiyr(:,:)   ! for plant resp
     real(r8),pointer:: nppyr(:,:)   ! for soil resp

!!! soil carbon pools  units(g C/m2)
     real(r8),pointer:: ccd(:,:)
     real(r8),pointer:: cfmd(:,:)
     real(r8),pointer:: cfsd(:,:)
     real(r8),pointer:: cm(:,:)
     real(r8),pointer:: cp(:,:)
     real(r8),pointer:: cs(:,:)
     real(r8),pointer:: csm(:,:)
     real(r8),pointer:: csmd(:,:)
     real(r8),pointer:: cssd(:,:)
     !     real(r8),pointer:: p_Vcmax(:)
     !     real(r8),pointer:: p_q10(:)
     !     real(r8),pointer:: p_drainage(:)
     !     real(r8),pointer:: p_beta(:)
     !     real(r8),pointer:: p_Ksat(:)
     !     real(r8),pointer:: p_b(:)
     !     real(r8),pointer:: u_Vcmax(:)
     !     real(r8),pointer:: u_q10(:)
     !     real(r8),pointer:: u_drainage(:)
     !     real(r8),pointer:: u_beta(:)
     !     real(r8),pointer:: u_Ksat(:)
     !     real(r8),pointer:: u_b(:)

     real(r8),pointer:: lai(:,:)     ! for photosynthesis
     real(r8),pointer:: Vcmax(:,:)   ! for data assimilation
     !     character,pointer:: name(:)
  end type surf

  type(surf),save,target,public:: bound   ! boundary conditions

  !******************for assimilation and parameter optimization/
  type,public::para
     real(r8),pointer:: p_Vcmax(:)
     real(r8),pointer:: p_VJ_slope(:)
     real(r8),pointer:: p_q10(:)
     real(r8),pointer:: p_sif_alpha(:)
     real(r8),pointer:: p_sif_beta(:)
     real(r8),pointer:: p_taweff(:)
     real(r8),pointer:: p_D0(:)
     real(r8),pointer:: p_Ksat_scalar(:)
     real(r8),pointer:: p_b_scalar(:)
     real(r8),pointer:: p_f_leaf
     real(r8),pointer:: p_kc25
     real(r8),pointer:: p_ko25
     real(r8),pointer:: p_tau25
 !    real(r8),pointer:: p_f_lr
     real(r8),pointer:: p_agb2vod

     real(r8),pointer:: u_Vcmax(:)
     real(r8),pointer:: u_VJ_slope(:)
     real(r8),pointer:: u_q10(:)
     real(r8),pointer:: u_sif_alpha(:)
     real(r8),pointer:: u_sif_beta(:)
     real(r8),pointer:: u_taweff(:)
     real(r8),pointer:: u_D0(:)
     real(r8),pointer:: u_Ksat_scalar(:)
     real(r8),pointer:: u_b_scalar(:)
     real(r8),pointer:: u_f_leaf
     real(r8),pointer:: u_kc25
     real(r8),pointer:: u_ko25
     real(r8),pointer:: u_tau25
!     real(r8),pointer:: u_f_lr
     real(r8),pointer:: u_agb2vod

  end type para

  type(para),save,target,public:: assim   ! optimization of parameters


  !*******************************************************************
  !------------------Soil Status----------
  !*******************************************************************
  type,public  :: soils
     integer,pointer   ::  n_layer(:)
     real(r8),pointer  ::  Zp(:,:)
     real(r8),pointer  ::  Zsp(:,:)
     real(r8),pointer  ::  r_rain_g(:,:)
     real(r8),pointer  ::  r_drainage(:,:)
     real(r8),pointer  ::  r_root_decay(:,:)
     real(r8),pointer  ::  psi_min(:,:)
     real(r8),pointer  ::  alpha(:,:)
     real(r8),pointer  ::  f_soilwater(:,:)


     real(r8),pointer  ::  d_soil(:,:)
     real(r8),pointer  ::  f_root(:,:,:)
     real(r8),pointer  ::  dt(:,:,:)
     real(r8),pointer  ::  thermal_cond(:,:,:)
     real(r8),pointer  ::  theta_vfc(:,:,:)
     real(r8),pointer  ::  theta_vwp(:,:,:)
     real(r8),pointer  ::  fei(:,:,:)
     real(r8),pointer  ::  Ksat(:,:,:)
     real(r8),pointer  ::  psi_sat(:,:,:)
     real(r8),pointer  ::  b(:,:,:)
     real(r8),pointer  ::  density_soil(:,:)   !!(npoints,0:Max_LayerS-1)
     real(r8),pointer  ::  f_org(:,:,:)
     real(r8),pointer  ::  ice_ratio(:,:,:)
     real(r8),pointer  ::  thetam(:,:,:)
     real(r8),pointer  ::  thetam_prev(:,:,:)
     real(r8),pointer  ::  temp_soil_p(:,:,:)
     real(r8),pointer  ::  temp_soil_c(:,:,:)
     real(r8),pointer  ::  f_ice(:,:,:)
     real(r8),pointer  ::  psim(:,:,:)
     real(r8),pointer  ::  thetab(:,:,:)
     real(r8),pointer  ::  psib(:,:,:)
     real(r8),pointer  ::  r_waterflow(:,:,:)
     real(r8),pointer  ::  km(:,:,:)
     real(r8),pointer  ::  kb(:,:,:)
     real(r8),pointer  ::  KK(:,:,:)
     real(r8),pointer  ::  Cs(:,:,:)
     real(r8),pointer  ::  lambda(:,:,:)
     real(r8),pointer  ::  Ett(:,:,:)
     real(r8),pointer  ::  G(:,:,:)
  end type soils

  type(soils),target,save,public :: soilstat


  !*********************************************************
  !-----------Interest Variables For output-----------------
  !*********************************************************
  type,public :: res
     real(r8),pointer::  GPPpft(:,:)
     real(r8),pointer::  SIFpft(:,:)
     real(r8),pointer::  SIFpft_sat(:,:)   !accord with the satellite data
     real(r8),pointer::  NPPpft(:,:)
     real(r8),pointer::  NEPpft(:,:)
     real(r8),pointer::  SHpft(:,:)
     real(r8),pointer::  LHpft(:,:)
     real(r8),pointer::  Transpft(:,:)
     real(r8),pointer::  Evappft(:,:)
     real(r8),pointer::  Net_Radpft(:,:)
     real(r8),pointer::  GPP(:)
     real(r8),pointer::  SIF(:)
     real(r8),pointer::  SIF_sat(:)
     real(r8),pointer::  NPP(:)
     real(r8),pointer::  NEP(:)
     real(r8),pointer::  LAIpft(:,:)
     real(r8),pointer::  LAI(:)
     real(r8),pointer::  SH(:)
     real(r8),pointer::  LH(:)
     real(r8),pointer::  Trans(:)
     real(r8),pointer::  Evap(:)
     real(r8),pointer::  Net_Rad(:)
     real(r8),pointer::  Thetampft(:,:)
     real(r8),pointer::  Thetam(:)
     real(r8),pointer::  fAPARpft(:,:)
     real(r8),pointer::  fAPAR(:)
     real(r8),pointer::  VODpft(:,:)
     real(r8),pointer::  VOD(:)
     real(r8),pointer::  COS_fluxpft(:,:)
     real(r8),pointer::  COS_flux(:)
     real(r8),pointer::  NPP_yr_acc(:,:)
  end type res

  type(res),save,target,public:: output

end module bepstype
