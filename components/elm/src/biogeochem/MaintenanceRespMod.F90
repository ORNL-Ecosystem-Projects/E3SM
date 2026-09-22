module MaintenanceRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding maintenance respiration routines for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varpar          , only : nlevgrnd
  use shr_const_mod       , only : SHR_CONST_TKFRZ
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use pftvarcon           , only : iscft
  use SharedParamsMod   , only : ParamsShareInst
  use VegetationPropertiesType      , only : veg_vp
  use SoilStateType       , only : soilstate_type
  use CanopyStateType     , only : canopystate_type
  use CNStateType         , only : cnstate_type
  use TemperatureType     , only : temperature_type
  use PhotosynthesisType  , only : photosyns_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnDataType      , only : col_es
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_es, veg_cs, veg_cf, veg_ns
  use TopounitType        , only : top_pp
  use elm_varctl          , only : use_humhol
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MaintenanceResp
  public :: readMaintenanceRespParams

   type, private :: MaintenanceRespParamsType
      real(r8):: br_mr        !base rate for maintenance respiration(gC/gN/s)
      real(r8):: mr_acclim_warming_frac ! fraction of warming acclimated by woody respiration
      integer :: mr_acclim_spinup_years ! spinup years used to define the temperature baseline
   end type MaintenanceRespParamsType

  !type(MaintenanceRespParamsType),private ::  MaintenanceRespParamsInst
  real(r8), public :: br_mr_Inst
  real(r8), public :: mr_acclim_warming_frac_Inst
  integer, public :: mr_acclim_spinup_years_Inst
  real(r8), public :: cpool_target_wood_frac_Inst
  real(r8), public :: cpool_target_leafroot_frac_Inst
  real(r8), public :: cpool_xr_scale_max_Inst
  !$acc declare create(br_mr_Inst)
  !$acc declare create(mr_acclim_warming_frac_Inst)
  !$acc declare create(mr_acclim_spinup_years_Inst)
  !$acc declare create(cpool_target_wood_frac_Inst)
  !$acc declare create(cpool_target_leafroot_frac_Inst)
  !$acc declare create(cpool_xr_scale_max_Inst)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine readMaintenanceRespParams ( ncid )
     !
     ! !DESCRIPTION:
     ! Read parameters
     !
     ! !USES:
     use ncdio_pio , only : file_desc_t,ncd_io
     !
     ! !ARGUMENTS:
     implicit none
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     !
     ! !LOCAL VARIABLES:
     character(len=32)  :: subname = 'MaintenanceRespParamsType'
     character(len=100) :: errCode = '-Error reading in parameters file:'
     logical            :: readv ! has variable been read in or not
     real(r8)           :: tempr ! temporary to read in constant
     character(len=100) :: tString ! temp. var for reading
     !-----------------------------------------------------------------------

     tString='br_mr'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     br_mr_Inst = tempr

     tString='mr_acclim_warming_frac'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if (.not. readv) then
        mr_acclim_warming_frac_Inst = 0._r8
     else
        mr_acclim_warming_frac_Inst = tempr
     end if

     tString='mr_acclim_spinup_years'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if (.not. readv) then
        mr_acclim_spinup_years_Inst = 0
     else
        mr_acclim_spinup_years_Inst = nint(tempr)
     end if

     tString='cpool_target_wood_frac'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if (.not. readv) then
        cpool_target_wood_frac_Inst = 0.03_r8
     else
        cpool_target_wood_frac_Inst = tempr
     end if
     if (cpool_target_wood_frac_Inst < 0._r8) then
        call endrun(msg='ERROR: cpool_target_wood_frac must be nonnegative'//errMsg(__FILE__, __LINE__))
     end if

     tString='cpool_target_leafroot_frac'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if (.not. readv) then
        cpool_target_leafroot_frac_Inst = 0.10_r8
     else
        cpool_target_leafroot_frac_Inst = tempr
     end if
     if (cpool_target_leafroot_frac_Inst < 0._r8) then
        call endrun(msg='ERROR: cpool_target_leafroot_frac must be nonnegative'//errMsg(__FILE__, __LINE__))
     end if

     tString='cpool_xr_scale_max'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if (.not. readv) then
        cpool_xr_scale_max_Inst = 10._r8
     else
        cpool_xr_scale_max_Inst = tempr
     end if
     if (cpool_xr_scale_max_Inst <= 0._r8) then
        call endrun(msg='ERROR: cpool_xr_scale_max must be positive'//errMsg(__FILE__, __LINE__))
     end if

   end subroutine readMaintenanceRespParams

  !-----------------------------------------------------------------------
  ! FIX(SPM,032414) this shouldn't even be called with ED on.
  !
  subroutine MaintenanceResp(bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       canopystate_vars, soilstate_vars, photosyns_vars, cnstate_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_soilp       ! number of soil points in patch filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j ! indices
    integer :: fp    ! soil filter patch index
    integer :: fc    ! soil filter column index
    real(r8):: br_mr ! base rate (gC/gN/s)
    real(r8):: br_mr_woody ! acclimated woody base rate (gC/gN/s)
    real(r8):: spinup_temp_offset ! annual temperature departure from spinup baseline (K)
    real(r8):: q10   ! temperature dependence
    real(r8):: tc    ! temperature correction, 2m air temp (unitless)
    real(r8):: tc_root ! PFT-specific root temperature correction (unitless)
    real(r8):: cpool_target ! target nonstructural C pool (gC m-2)
    real(r8):: cpool_relative ! actual cpool relative to its target (unitless)
    real(r8):: tcsoi(bounds%begc:bounds%endc,nlevgrnd) ! temperature correction by soil layer (unitless)
    !-----------------------------------------------------------------------

    associate(                                                        &
         ivt            =>    veg_pp%itype                             , & ! Input:  [integer  (:)   ]  patch vegetation type
         woody          =>    veg_vp%woody                      , & ! Input:  [real(r8) (:)   ]  woody lifeform flag (0 = non-woody, 1 = tree, 2 = shrub)
         br_xr          =>    veg_vp%br_xr                      , & ! Input:  [real(r8) (:)   ]  base rate for excess respiration
         br_mr_pft      =>    veg_vp%br_mr_pft                  , & ! Input:  [real(r8) (:)   ]  PFT-specific base rate
         q10_mr_pft     =>    veg_vp%q10_mr_pft                 , & ! Input:  [real(r8) (:)   ]  PFT-specific Q10
         frac_veg_nosno =>    canopystate_vars%frac_veg_nosno_patch , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         laisun         =>    canopystate_vars%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha         =>    canopystate_vars%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index

         rootfr         =>    soilstate_vars%rootfr_patch           , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)
         annavg_t2m     =>    cnstate_vars%annavg_t2m_patch         , & ! Input:  [real(r8) (:)   ]  annual average air temperature (K)
         spinup_t       =>    cnstate_vars%spinup_t_patch           , & ! Input:  [real(r8) (:)   ]  spinup reference temperature (K)
         spinup_t_nyears =>   cnstate_vars%spinup_t_nyears_patch    , & ! Input:  [integer  (:)   ]  years in reference

         t_soisno       =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         t_ref2m        =>    veg_es%t_ref2m          , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)

         lmrsun         =>    photosyns_vars%lmrsun_patch           , & ! Input:  [real(r8) (:)   ]  sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
         lmrsha         =>    photosyns_vars%lmrsha_patch           , & ! Input:  [real(r8) (:)   ]  shaded leaf maintenance respiration rate (umol CO2/m**2/s)

         cpool          =>    veg_cs%cpool          , & ! Input: [real(r8) (:)   ]   plant carbon pool (gC m-2)
         leafc          =>    veg_cs%leafc          , & ! Input: [real(r8) (:)   ]   leaf carbon (gC m-2)
         frootc         =>    veg_cs%frootc         , & ! Input: [real(r8) (:)   ]   fine-root carbon (gC m-2)
         livestemc     =>    veg_cs%livestemc      , & ! Input: [real(r8) (:)   ]   live stem carbon (gC m-2)
         deadstemc     =>    veg_cs%deadstemc      , & ! Input: [real(r8) (:)   ]   structural stem carbon (gC m-2)
         livecrootc    =>    veg_cs%livecrootc     , & ! Input: [real(r8) (:)   ]   live coarse-root carbon (gC m-2)
         deadcrootc    =>    veg_cs%deadcrootc     , & ! Input: [real(r8) (:)   ]   structural coarse-root carbon (gC m-2)

         leaf_mr        =>    veg_cf%leaf_mr         , & ! Output: [real(r8) (:)   ]
         froot_mr       =>    veg_cf%froot_mr        , & ! Output: [real(r8) (:)   ]
         livestem_mr    =>    veg_cf%livestem_mr     , & ! Output: [real(r8) (:)   ]
         livecroot_mr   =>    veg_cf%livecroot_mr    , & ! Output: [real(r8) (:)   ]
         grain_mr       =>    veg_cf%grain_mr        , & ! Output: [real(r8) (:)   ]
         xr             =>    veg_cf%xr              , & ! Output: [real(r8) (:)   ]  (gC/m2) respiration of excess C

         frootn         =>    veg_ns%frootn       , & ! Input:  [real(r8) (:)   ]  (gN/m2) fine root N
         livestemn      =>    veg_ns%livestemn    , & ! Input:  [real(r8) (:)   ]  (gN/m2) live stem N
         livecrootn     =>    veg_ns%livecrootn   , & ! Input:  [real(r8) (:)   ]  (gN/m2) live coarse root N
         grainn         =>    veg_ns%grainn         & ! Output: [real(r8) (:)   ]  (kgN/m2) grain N
         )

      ! base rate for maintenance respiration is from:
      ! M. Ryan, 1991. Effects of climate change on plant respiration.
      ! Ecological Applications, 1(2), 157-167.
      ! Original expression is br = 0.0106 molC/(molN h)
      ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
      ! set constants
      br_mr = br_mr_Inst

      ! Peter Thornton: 3/13/09
      ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
      ! to improve seasonal cycle of atmospheric CO2 concentration in global
      ! simulatoins

      ! Set Q10 from SharedParamsMod
      Q10 = ParamsShareInst%Q10_mr

      ! column loop to calculate temperature factors in each soil layer
      do j=1,nlevgrnd
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! calculate temperature corrections for each soil layer, for use in
            ! estimating fine root maintenance respiration with depth
            tcsoi(c,j) = Q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)

         end do
      end do

      ! patch loop for leaves and live wood
      do fp = 1, num_soilp
         p = filter_soilp(fp)
         br_mr = br_mr_Inst
         br_mr_woody = br_mr

         ! calculate maintenance respiration fluxes in
         ! gC/m2/s for each of the live plant tissues.
         ! Leaf and live wood MR

         tc = Q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
         if (use_humhol .and. top_pp%peat_depth(veg_pp%topounit(p)) > 0._r8) then
            br_mr = br_mr_pft(ivt(p))
            tc = q10_mr_pft(ivt(p))**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
            br_mr_woody = br_mr
            if (mr_acclim_warming_frac_Inst > 0._r8 .and. spinup_t_nyears(p) > 0) then
               spinup_temp_offset = annavg_t2m(p) - spinup_t(p)
               br_mr_woody = br_mr * q10_mr_pft(ivt(p))** &
                    ((-spinup_temp_offset * mr_acclim_warming_frac_Inst) / 10._r8)
            end if
         end if
         if (frac_veg_nosno(p) == 1) then
            leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_r8 + &
                         lmrsha(p) * laisha(p) * 12.011e-6_r8

         else !nosno
             leaf_mr(p) = 0._r8

         end if

         if (woody(ivt(p)) >= 1.0_r8) then
            livestem_mr(p) = livestemn(p)*br_mr_woody*tc
            livecroot_mr(p) = livecrootn(p)*br_mr_woody*tc
         else if (iscft(ivt(p)) .and. livestemn(p) .gt. 0._r8) then
            livestem_mr(p) = livestemn(p)*br_mr*tc
            grain_mr(p) = grainn(p)*br_mr*tc
         end if
         if (br_xr(ivt(p)) .gt. 1e-9_r8) then
            if (top_pp%peat_depth(veg_pp%topounit(p)) > 0._r8) then
               ! Treat cpool as a biomass-scaled nonstructural-C reserve on
               ! peatland topounits. The effective whole-organ wood fraction
               ! includes both live and structural wood because br_xr acts on
               ! one accessible patch-level reserve rather than organ pools.
               cpool_target = cpool_target_wood_frac_Inst * &
                    (max(livestemc(p), 0._r8) + max(deadstemc(p), 0._r8) + &
                     max(livecrootc(p), 0._r8) + max(deadcrootc(p), 0._r8)) + &
                    cpool_target_leafroot_frac_Inst * &
                    (max(leafc(p), 0._r8) + max(frootc(p), 0._r8))

               ! br_xr remains a fractional turnover rate (s-1). At the
               ! target, the existing XR equation is recovered. The turnover
               ! coefficient scales linearly with cpool/target, so at twice
               ! the target the coefficient is 2*br_xr and total XR is four
               ! times its value at the target. Before structural biomass is
               ! established but cpool remains after senescence, use the same
               ! bounded maximum rather than divide by a vanishing target.
               if (cpool_target > tiny(1._r8)) then
                  cpool_relative = min(cpool_xr_scale_max_Inst, &
                       max(cpool(p), 0._r8) / cpool_target)
               else if (cpool(p) > 0._r8) then
                  cpool_relative = cpool_xr_scale_max_Inst
               else
                  cpool_relative = 1._r8
               end if
               xr(p) = cpool(p) * br_xr(ivt(p)) * cpool_relative * tc
            else
               xr(p) = cpool(p) * br_xr(ivt(p)) * tc
            end if
            !xr_above(p) = xr(p) * (leafn(p) + livestemn(p)) / &
            !          (leafn(p) + livestemn(p) + frootn(p))
            !xr_below(p) = xr(p) - xr_above(p)
         else
            xr(p) = 0._r8
            !xr_above(p) = 0._r8
            !xr_below(p) = 0._r8
         end if
      end do

      ! soil and patch loop for fine root

      do j = 1,nlevgrnd
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = veg_pp%column(p)
            br_mr = br_mr_Inst
            tc_root = tcsoi(c,j)

            ! Fine root MR
            ! rootfr(j) sums to 1.0 over all soil layers, and
            ! describes the fraction of root mass that is in each
            ! layer.  This is used with the layer temperature correction
            ! to estimate the total fine root maintenance respiration as a
            ! function of temperature and N content.
            if (use_humhol .and. top_pp%peat_depth(veg_pp%topounit(p)) > 0._r8) then
               br_mr = br_mr_pft(ivt(p))
               tc_root = q10_mr_pft(ivt(p))** &
                    ((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
            end if
            froot_mr(p) = froot_mr(p) + frootn(p)*br_mr*tc_root*rootfr(p,j)
         end do
      end do

    end associate

  end subroutine MaintenanceResp

end module MaintenanceRespMod
