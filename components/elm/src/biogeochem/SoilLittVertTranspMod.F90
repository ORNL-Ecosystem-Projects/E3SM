module SoilLittVertTranspMod

  !-----------------------------------------------------------------------
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varctl             , only : iulog, use_c13, use_c14, spinup_state, use_vertsoilc
  use elm_varctl             , only : use_microbe_methane
  use elm_varctl             , only : use_peatland_vertical_transport
  use elm_varctl             , only : use_peatland_compaction_profile
  use elm_varctl             , only : use_microbe_aqueous_transport
  use elm_varcon             , only : secspday
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use CNDecompCascadeConType , only : decomp_cascade_con
  use MicrobeDecompMod       , only : MicrobeDecompParamsInst
  use CanopyStateType        , only : canopystate_type
  use CNStateType            , only : cnstate_type
  use elm_varctl             , only : nu_com
  use ColumnDataType         , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType         , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType         , only : col_ns, col_nf, col_ps, col_pf
  use ColumnType             , only : col_pp
  use TopounitType           , only : top_pp
  use timeinfoMod
  !
  implicit none
  save
  !
  public :: SoilLittVertTransp
  public :: createLitterTransportList
  public :: readSoilLittVertTranspParams
  private :: calc_diffus_advflux

  type, public :: SoilLittVertTranspParamsType
     real(r8)  :: som_diffus                  ! Soil organic matter diffusion
     real(r8)  :: cryoturb_diffusion_k        ! The cryoturbation diffusive constant cryoturbation to the active layer thickness
     real(r8)  :: max_altdepth_cryoturbation  ! (m) maximum active layer thickness for cryoturbation to occur
  end type SoilLittVertTranspParamsType

  type(SoilLittVertTranspParamsType), public ::  SoilLittVertTranspParamsInst
  !$acc declare create(SoilLittVertTranspParamsInst)

  type, public :: ConcTransportType
     !! Type that points to decomposition pools for vertical transport calculations

     real(r8), pointer :: conc_ptr(:,:,:) => null()
     real(r8), pointer :: src_ptr(:,:,:)  => null()
     real(r8), pointer :: trcr_tend_ptr(:,:,:) => null()
  end type ConcTransportType
  type(ConcTransportType), public, allocatable :: transport_ptr_list(:)
  !$acc declare create(transport_ptr_list(:))

  real(r8), public :: som_adv_flux =  0._r8
  !$acc declare create(som_adv_flux)
  real(r8), public :: peat_som_adv_flux = 0.0004_r8 / (secspday * 365._r8)
  !$acc declare create(peat_som_adv_flux)
  real(r8), public :: peat_som_diffus = 0._r8
  !$acc declare create(peat_som_diffus)
  real(r8), public :: peat_adv_reference_depth = 3._r8
  !$acc declare create(peat_adv_reference_depth)
  real(r8), public :: peat_compaction_surface_density = 25._r8
  !$acc declare create(peat_compaction_surface_density)
  real(r8), public :: peat_compaction_deep_density = 100._r8
  !$acc declare create(peat_compaction_deep_density)
  real(r8), public :: peat_compaction_efolding_depth = 0.25_r8
  !$acc declare create(peat_compaction_efolding_depth)
  real(r8), public :: peat_compaction_timescale_years = 1._r8
  !$acc declare create(peat_compaction_timescale_years)
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  !$acc declare create(max_depth_cryoturb)
  !-----------------------------------------------------------------------

contains

   subroutine createLitterTransportList()
      ! This subroutine creates a list that will point to the
      ! litter/som fields needed for the vertical transport
      ! calculations.

      implicit none

      integer :: ntype

      ntype = 3
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif

      allocate(transport_ptr_list(ntype))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! C
      transport_ptr_list(1)%conc_ptr        => col_cs%decomp_cpools_vr
      transport_ptr_list(1)%src_ptr         => col_cf%decomp_cpools_sourcesink
      transport_ptr_list(1)%trcr_tend_ptr   => col_cf%decomp_cpools_transport_tendency
      ! N
      transport_ptr_list(2)%conc_ptr        => col_ns%decomp_npools_vr
      transport_ptr_list(2)%src_ptr         => col_nf%decomp_npools_sourcesink
      transport_ptr_list(2)%trcr_tend_ptr   => col_nf%decomp_npools_transport_tendency
      ! P
      transport_ptr_list(3)%conc_ptr        => col_ps%decomp_ppools_vr
      transport_ptr_list(3)%src_ptr         => col_pf%decomp_ppools_sourcesink
      transport_ptr_list(3)%trcr_tend_ptr   => col_pf%decomp_ppools_transport_tendency
      ! c13 and c14 if there
      if(use_c14 .and. use_c13) then
         !
         transport_ptr_list(4)%conc_ptr       => c13_col_cs%decomp_cpools_vr
         transport_ptr_list(4)%src_ptr        => c13_col_cf%decomp_cpools_sourcesink
         transport_ptr_list(4)%trcr_tend_ptr  => c13_col_cf%decomp_cpools_transport_tendency
         !
         transport_ptr_list(5)%conc_ptr       => c14_col_cs%decomp_cpools_vr
         transport_ptr_list(5)%src_ptr        => c14_col_cf%decomp_cpools_sourcesink
         transport_ptr_list(5)%trcr_tend_ptr  => c14_col_cf%decomp_cpools_transport_tendency
      else
         if(use_c13) then
            transport_ptr_list(4)%conc_ptr       => c13_col_cs%decomp_cpools_vr
            transport_ptr_list(4)%src_ptr        => c13_col_cf%decomp_cpools_sourcesink
            transport_ptr_list(4)%trcr_tend_ptr  => c13_col_cf%decomp_cpools_transport_tendency
         end if
         if (use_c14) then
            transport_ptr_list(4)%conc_ptr       => c14_col_cs%decomp_cpools_vr
            transport_ptr_list(4)%src_ptr        => c14_col_cf%decomp_cpools_sourcesink
            transport_ptr_list(4)%trcr_tend_ptr  => c14_col_cf%decomp_cpools_transport_tendency
         end if
      end if

   end subroutine createLitterTransportList

   subroutine cleanupLitterTransportList()
      !! Nullifies pointers and frees allocated memory
      integer :: i

      if (allocated(transport_ptr_list)) then
         ! Loop through each element and nullify pointers
         do i = 1, size(transport_ptr_list)
            if (associated(transport_ptr_list(i)%conc_ptr)) then
               nullify(transport_ptr_list(i)%conc_ptr)
            endif
            if (associated(transport_ptr_list(i)%src_ptr)) then
               nullify(transport_ptr_list(i)%src_ptr)
            endif
            if (associated(transport_ptr_list(i)%trcr_tend_ptr)) then
               nullify(transport_ptr_list(i)%trcr_tend_ptr)
            endif
         end do

         ! Deallocate the transport_ptr_list array itself
         deallocate(transport_ptr_list)
      endif

   end subroutine cleanupLitterTransportList


  !-----------------------------------------------------------------------
  subroutine readSoilLittVertTranspParams ( ncid )
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'SoilLittVertTranspType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in parameters
    !
    ! REMOVE THESE?
    tString='som_diffus'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    !SoilLittVertTranspParamsInst%som_diffus=tempr
    ! FIX(SPM,032414) - can't be pulled out since division makes things not bfb
    SoilLittVertTranspParamsInst%som_diffus = 1e-4_r8 / (secspday * 365._r8)

    tString='cryoturb_diffusion_k'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    !SoilLittVertTranspParamsInst%cryoturb_diffusion_k=tempr
    !FIX(SPM,032414) Todo.  This constant cannot be on file since the divide makes things
    !SPM Todo.  This constant cannot be on file since the divide makes things
    !not bfb

    SoilLittVertTranspParamsInst%cryoturb_diffusion_k = 5e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr

    tString='max_altdepth_cryoturbation'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    SoilLittVertTranspParamsInst%max_altdepth_cryoturbation=tempr

    !$acc enter data copyin(SoilLittVertTranspParamsInst)
  end subroutine readSoilLittVertTranspParams

  function aaa(pe) result(res)
     !$acc routine seq
     implicit none
     real(r8) :: res
     real(r8) :: pe
     res =  max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)
  end function

  !-----------------------------------------------------------------------
  subroutine SoilLittVertTransp(num_soilc, filter_soilc,   &
       canopystate_vars, cnstate_vars )
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools
    ! calculated in the CarbonStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
    use elm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use elm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: pe          ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1  ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1  ! Harmonic mean of diffusivity
    real(r8) :: d_p1_zp1    ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1    ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: pe_p1       ! Peclet # for next j
    real(r8) :: pe_m1       ! Peclet # for previous j
    real(r8) :: dz_node,dz_nodep1          ! difference between nodes
    real(r8) :: a_p_0
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l,t  ! indices
    integer  :: jtop(num_soilc)    ! top level at each column
    real(r8) :: spinup_term                  ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8), parameter :: epsilon=1.e-30     ! small number
    !!added to remove arrays:
    real(r8) :: diffus_j, diffus_jm1, diffus_jp1 ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux_j,adv_flux_jm1, adv_flux_jp1        ! advective flux (m/s)  (includes spinup correction, if any)
    real(r8) :: a_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "r" vector for tridiagonal solution
    real(r8) :: conc_trcr(num_soilc,0:nlevdecomp+1,ndecomp_pools)                  !
    real(r8) :: bet
    real(r8) :: gam(0:nlevdecomp+1)
    real(r8) :: peat_depth_factor
    real(r8) :: peat_physical_density
    real(r8) :: peat_target_density
    real(r8) :: peat_pool_factor
    real(r8) :: peat_excess_flux
    real(r8) :: peat_compaction_timescale_seconds
    real(r8) :: dom_diffusion_multiplier
    logical  :: peat_dynamic_column(num_soilc)
    logical  :: transport_pool(num_soilc,ndecomp_pools)
    real(r8) :: peat_c_density(num_soilc,nlevdecomp)
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------

    ! Set statement functions
    associate(                                                      &
         is_cwd           => decomp_cascade_con%is_cwd            , & ! Input:  [logical (:)    ]  TRUE => pool is a cwd pool
         is_dissolved     => decomp_cascade_con%is_dissolved      , &
         is_microbial     => decomp_cascade_con%is_microbial_biomass, &
         spinup_factor    => decomp_cascade_con%spinup_factor     , & ! Input:  [real(r8) (:)   ]  spinup accelerated decomposition factor, used to accelerate transport as well

         decomp_cpools_vr => col_cs%decomp_cpools_vr             , &
         decomp_csources  => col_cf%decomp_cpools_sourcesink     , &

         altmax           => canopystate_vars%altmax_col          , & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear  => canopystate_vars%altmax_lastyear_col , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw

         som_adv_coef     => cnstate_vars%som_adv_coef_col       , & ! Output: [real(r8) (:,:) ]  SOM advective flux (m/s)
         som_diffus_coef  => cnstate_vars%som_diffus_coef_col    ,  & ! Output: [real(r8) (:,:) ]  SOM diffusivity due to bio/cryo-turbation (m2/s)
         peat_c_density_diag => cnstate_vars%peat_c_density_col  , &
         peat_c_target_diag => cnstate_vars%peat_c_target_density_col, &
         peat_c_burial_flux_diag => cnstate_vars%peat_c_burial_flux_col, &
         ! !Set parameters of vertical mixing of SOM
          som_diffus                 => SoilLittVertTranspParamsInst%som_diffus   , &
          cryoturb_diffusion_k       => SoilLittVertTranspParamsInst%cryoturb_diffusion_k  , &
          max_altdepth_cryoturbation => SoilLittVertTranspParamsInst%max_altdepth_cryoturbation &
         )
      
      dom_diffusion_multiplier = 1._r8
      if (use_microbe_methane) then
         dom_diffusion_multiplier = MicrobeDecompParamsInst%dom_som_diffusion_multiplier
      end if

      !$acc enter data copyin(dom_diffusion_multiplier)
      !$acc enter data create(a_tri(:,:,:),b_tri(:,:,:),&
      !$acc     c_tri(:,:,:),r_tri(:,:,:), &
      !$acc     conc_trcr(:,:,:), gam(:), &
      !$acc     peat_dynamic_column(:), transport_pool(:,:), &
      !$acc     peat_c_density(:,:) )
      ntype = 3
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif
   
      !$acc enter data create(spinup_term, i_type) 
      spinup_term = 1._r8
      !$acc update device(spinup_term)

      peat_compaction_timescale_seconds = peat_compaction_timescale_years * &
           secspday * 365._r8
      !$acc enter data copyin(peat_compaction_timescale_seconds)

      ! Build the physical-equivalent solid-C density used by the peat
      ! capacity-overflow calculation. In AD spinup the prognostic slow-pool
      ! concentrations are deliberately reduced, so reconstruct the physical
      ! density with the same pool factors used by vertical transport.
      !$acc parallel loop independent gang vector default(present) private(c,t)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         t = col_pp%topounit(c)
         peat_dynamic_column(fc) = use_peatland_compaction_profile .and. &
              top_pp%peat_depth(t) > 0._r8
      end do

      !$acc parallel loop independent gang default(present)
      do j = 1,nlevdecomp
         !$acc loop vector independent &
         !$acc& private(c,peat_physical_density,peat_target_density, &
         !$acc& peat_pool_factor,s)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            peat_physical_density = 0._r8
            peat_target_density = 0._r8
            if (peat_dynamic_column(fc)) then
               !$acc loop seq
               do s = 1,ndecomp_pools
                  if (.not. is_dissolved(s)) then
                     peat_pool_factor = 1._r8
                     if (spinup_state == 1) then
                        peat_pool_factor = spinup_factor(s)
                        if (spinup_factor(s) > 1._r8 .and. year_curr >= 40) then
                           peat_pool_factor = peat_pool_factor / &
                                max(cnstate_vars%scalaravg_col(c,j), epsilon)
                        end if
                     end if
                     peat_physical_density = peat_physical_density + &
                          max(decomp_cpools_vr(c,j,s) + decomp_csources(c,j,s), 0._r8) * &
                          peat_pool_factor
                  end if
               end do
               peat_target_density = 1000._r8 * &
                    (peat_compaction_surface_density + &
                     (peat_compaction_deep_density - &
                      peat_compaction_surface_density) * &
                     (1._r8 - exp(-max(zsoi(j),0._r8) / &
                      peat_compaction_efolding_depth)))
            end if
            peat_c_density(fc,j) = peat_physical_density
            peat_c_density_diag(c,j) = peat_physical_density
            peat_c_target_diag(c,j) = peat_target_density
            peat_c_burial_flux_diag(c,j) = 0._r8
         end do
      end do

      !$acc parallel loop independent gang default(present)
      do s = 1,ndecomp_pools
         !$acc loop vector independent private(c)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (peat_dynamic_column(fc)) then
               ! Capacity overflow buries the complete solid peat matrix.
               transport_pool(fc,s) = .not. is_dissolved(s)
            else
               ! Preserve the standard ELM pool selection outside dynamic
               ! peat columns.
               transport_pool(fc,s) = .not. is_cwd(s) .and. &
                    .not. is_microbial(s) .and. &
                    .not. (use_microbe_aqueous_transport .and. is_dissolved(s))
            end if
         end do
      end do

      if (use_vertsoilc) then
         !------ first get diffusivity / advection terms -------!
         ! use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
         if (.not. use_peatland_vertical_transport) then
            ! Keep the standard path separate and unchanged for off-mode BFB.
            !$acc parallel loop independent gang default(present)
            do j = 1,nlevdecomp+1
               !$acc loop vector independent private(c)
               do fc = 1, num_soilc
                  c = filter_soilc (fc)
                  if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                     ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
                     ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
                     if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
                        som_diffus_coef(c,j) = cryoturb_diffusion_k
                        som_adv_coef(c,j) = 0._r8
                     else
                        som_diffus_coef(c,j) = max(cryoturb_diffusion_k * &
                             ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                             ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
                        som_adv_coef(c,j) = 0._r8
                     endif
                  elseif (  max(altmax(c), altmax_lastyear(c)) > 0._r8 ) then
                     ! constant advection, constant diffusion
                     som_adv_coef(c,j) = som_adv_flux
                     som_diffus_coef(c,j) = som_diffus
                  else
                     ! completely frozen soils--no mixing
                     som_adv_coef(c,j) = 0._r8
                     som_diffus_coef(c,j) = 0._r8
                  endif
               end do
            end do
         else
            !$acc parallel loop independent gang default(present)
            do j = 1,nlevdecomp+1
               !$acc loop vector independent &
               !$acc& private(c,t,l,peat_depth_factor,peat_physical_density, &
               !$acc& peat_target_density,peat_excess_flux)
               do fc = 1, num_soilc
                  c = filter_soilc (fc)
                  if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                     ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
                     if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
                        som_diffus_coef(c,j) = cryoturb_diffusion_k
                        som_adv_coef(c,j) = 0._r8
                     else
                        som_diffus_coef(c,j) = max(cryoturb_diffusion_k * &
                             ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                             ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)
                        som_adv_coef(c,j) = 0._r8
                     endif
                  elseif (max(altmax(c), altmax_lastyear(c)) > 0._r8) then
                     t = col_pp%topounit(c)
                     peat_depth_factor = max(top_pp%peat_depth(t), 0._r8) / &
                          peat_adv_reference_depth
                     if (peat_depth_factor > 0._r8) then
                        if (use_peatland_compaction_profile) then
                           ! The coefficient at j is the velocity through the
                           ! upper interface of layer j. There is no external
                           ! C supply at the surface and no loss through the
                           ! bottom of the decomposition column. Interior
                           ! velocities relax only donor-layer C above its
                           ! depth-dependent storage capacity.
                           som_adv_coef(c,j) = 0._r8
                           if (j > 1 .and. j <= nlevdecomp .and. &
                                zisoi(j-1) <= top_pp%peat_depth(t)) then
                              l = j - 1
                              peat_physical_density = peat_c_density(fc,l)
                              peat_target_density = 1000._r8 * &
                                   (peat_compaction_surface_density + &
                                    (peat_compaction_deep_density - &
                                     peat_compaction_surface_density) * &
                                    (1._r8 - exp(-max(zsoi(l),0._r8) / &
                                     peat_compaction_efolding_depth)))
                              if (peat_physical_density > peat_target_density) then
                                 peat_excess_flux = &
                                      (peat_physical_density - peat_target_density) * &
                                      dzsoi_decomp(l) / peat_compaction_timescale_seconds
                                 som_adv_coef(c,j) = peat_excess_flux / &
                                      max(peat_physical_density, epsilon)
                                 peat_c_burial_flux_diag(c,l) = peat_excess_flux
                              end if
                           end if
                        else
                           ! Compatibility path for the historical prescribed
                           ! velocity and total-peat-depth scaling.
                           som_adv_coef(c,j) = peat_som_adv_flux * peat_depth_factor
                        end if
                        som_diffus_coef(c,j) = peat_som_diffus
                     else
                        som_adv_coef(c,j) = som_adv_flux
                        som_diffus_coef(c,j) = som_diffus
                     end if
                  else
                     som_adv_coef(c,j) = 0._r8
                     som_diffus_coef(c,j) = 0._r8
                  endif
               end do
            end do
         end if
      endif
   
      !------ loop over litter/som types
      do i_type = 1, ntype

         !$acc update device(i_type)

         if (use_vertsoilc) then
            ! Set Pe (Peclet #) and D/dz throughout column

            !$acc parallel loop independent gang default(present)
            do s = 1, ndecomp_pools
               !$acc loop independent worker vector private(c)
               do fc = 1, num_soilc ! dummy terms here
                  if (transport_pool(fc,s)) then
                     c = filter_soilc (fc)
                     conc_trcr(fc,0,s) = 0._r8
                     conc_trcr(fc,nlevdecomp+1,s) = 0._r8

                     a_tri(fc,0,s) = 0._r8
                     b_tri(fc,0,s) = 1._r8
                     c_tri(fc,0,s) = -1._r8
                     r_tri(fc,0,s) = 0._r8

                     conc_trcr(fc,nlevdecomp+1,s) = transport_ptr_list(i_type)%conc_ptr(c,nlevdecomp+1,s)
                     a_tri(fc,nlevdecomp+1,s) = -1._r8
                     b_tri(fc,nlevdecomp+1,s) = 1._r8
                     c_tri(fc,nlevdecomp+1,s) = 0._r8
                     r_tri(fc,nlevdecomp+1,s) = 0._r8
                  end if
               end do
            end do

            !$acc parallel loop independent gang worker vector collapse(3) default(present) 
            do s = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     if (transport_pool(fc,s)) then

                        if ( spinup_state .eq. 1 ) then
                           ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
                           spinup_term = spinup_factor(s)
                        else
                           spinup_term = 1.
                        endif
                        conc_trcr(fc,j,s) = transport_ptr_list(i_type)%conc_ptr(c,j,s)
                        ! dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
                        ! dz_node_tracer is difference between cell centers
                        call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j), som_adv_coef(c,j), &
                                                 cnstate_vars%scalaravg_col(c,j),adv_flux_j, diffus_j)
                        if (is_dissolved(s)) diffus_j = diffus_j * dom_diffusion_multiplier

                        ! Calculate the D and F terms in the Patankar algorithm
                        if (j == 1) then
                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j+1), som_adv_coef(c,j+1), &
                                                   cnstate_vars%scalaravg_col(c,j+1),adv_flux_jp1, diffus_jp1)
                           if (is_dissolved(s)) diffus_jp1 = diffus_jp1 * dom_diffusion_multiplier
                           dz_nodep1 =  zsoi(j+1) - zsoi(j)
                           d_m1_zm1 = 0._r8
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_nodep1
                           if (diffus_jp1 > 0._r8 .and. diffus_j > 0._r8) then
                             d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus_j + w_p1 / diffus_jp1) ! Harmonic mean of diffus
                           else
                              d_p1 = 0._r8
                           endif

                           d_p1_zp1 = d_p1 / dz_nodep1
                           pe_m1 = 0._r8
                           pe_p1 = adv_flux_jp1 / d_p1_zp1 ! Peclet #

                           a_p_0 =  dzsoi_decomp(j) / dtime_mod
                           a_tri(fc,j,s) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j,s) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j,s) = -a_tri(fc,j,s) - c_tri(fc,j,s) + a_p_0
                           if (peat_dynamic_column(fc)) then
                              ! The historical operator assumed a constant
                              ! velocity. Include flux divergence for the
                              ! layer-varying capacity-overflow velocity.
                              b_tri(fc,j,s) = b_tri(fc,j,s) + &
                                   adv_flux_jp1 - adv_flux_j
                           end if
                           r_tri(fc,j,s) = transport_ptr_list(i_type)%src_ptr(c,j,s) * dzsoi_decomp(j) /dtime_mod + (a_p_0 - adv_flux_j) * conc_trcr(fc,j,s)
                        else
                          ! Use distance from j-1 node to interface with j divided by distance between nodes
                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j-1), som_adv_coef(c,j-1), &
                                                   cnstate_vars%scalaravg_col(c,j-1),adv_flux_jm1, diffus_jm1)
                          if (is_dissolved(s)) diffus_jm1 = diffus_jm1 * dom_diffusion_multiplier

                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j+1), som_adv_coef(c,j+1), &
                                                   cnstate_vars%scalaravg_col(c,j+1),adv_flux_jp1, diffus_jp1)
                          if (is_dissolved(s)) diffus_jp1 = diffus_jp1 * dom_diffusion_multiplier
                           ! Use distance from j-1 node to interface with j divided by distance between nodes
                           dz_node = zsoi(j) - zsoi(j-1)
                           w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node

                           if ( diffus_jm1 > 0._r8 .and. diffus_j > 0._r8) then
                              d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus_j + w_m1 / diffus_jm1) ! Harmonic mean of diffus
                           else
                              d_m1 = 0._r8
                           endif

                           dz_nodep1 = zsoi(j+1) - zsoi(j)
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_nodep1

                           if ( diffus_jp1 > 0._r8 .and. diffus_j > 0._r8) then
                              d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus_j + w_p1 / diffus_jp1) ! Harmonic mean of diffus
                           else
                              d_p1 = (1._r8 - w_m1) * diffus_j + w_p1 * diffus_jp1  ! Arithmetic mean of diffus
                           endif
                           d_m1_zm1 = d_m1 / dz_node
                           d_p1_zp1 = d_p1 / dz_nodep1
                           pe_m1 = adv_flux_j / d_m1_zm1 ! Peclet #
                           pe_p1 = adv_flux_jp1 / d_p1_zp1 ! Peclet #

                           a_p_0 =  dzsoi_decomp(j) / dtime_mod
                           a_tri(fc,j,s) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j,s) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j,s) = -a_tri(fc,j,s) - c_tri(fc,j,s) + a_p_0
                           if (peat_dynamic_column(fc)) then
                              b_tri(fc,j,s) = b_tri(fc,j,s) + &
                                   adv_flux_jp1 - adv_flux_j
                           end if
                           r_tri(fc,j,s) = transport_ptr_list(i_type)%src_ptr(c,j,s) * dzsoi_decomp(j) /dtime_mod + a_p_0 * conc_trcr(fc,j,s)
                        end if
                     end if
                  enddo ! fc
               enddo ! j; nlevdecomp
            end do ! s: ndecomp_pools

            ! subtract initial concentration and source terms for tendency calculation
            !$acc parallel loop independent collapse(3) gang vector default(present)
            do s = 1, ndecomp_pools
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     if (transport_pool(fc,s)) then
                        transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) = 0._r8 - (conc_trcr(fc,j,s) + transport_ptr_list(i_type)%src_ptr(c,j,s))
                     end if
                  end do
               end do
            end do

            ! Solve for the concentration profile for this time step

            !$acc parallel loop independent gang worker vector collapse(2) default(present) private(bet, gam(0:nlevdecomp+1))
            do s = 1, ndecomp_pools
               do fc = 1,num_soilc
                  if (transport_pool(fc,s)) then
                     bet = b_tri(fc,0,s)

                     !$acc loop seq
                     do j = 0, nlevdecomp+1
                        if (j == 0) then
                           conc_trcr(fc,j,s) = r_tri(fc,j,s) / bet
                        else
                           gam(j) = c_tri(fc,j-1,s) / bet
                           bet = b_tri(fc,j,s) - a_tri(fc,j,s) * gam(j)
                           conc_trcr(fc,j,s) = (r_tri(fc,j,s) - a_tri(fc,j,s)*conc_trcr(fc,j-1,s)) / bet
                        end if
                     end do

                     !$acc loop seq
                     do j = nlevdecomp,1,-1
                       conc_trcr(fc,j,s) = conc_trcr(fc,j,s) - gam(j+1) * conc_trcr(fc,j+1,s)
                     end do
                  end if
               end do
            end do

            ! add post-transport concentration to calculate tendency term
            !$acc parallel loop independent gang collapse(2) default(present)
            do s = 1, ndecomp_pools
               do j = 1, nlevdecomp
                  !$acc loop vector independent private(c)
                  do fc = 1, num_soilc
                     if (transport_pool(fc,s)) then
                        c = filter_soilc (fc)
                        transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) = &
                             (transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) + &
                              conc_trcr(fc,j,s)) / dtime_mod
                     end if
                  end do
                  !
               end do
            end do

            ! Pools not selected for matrix transport remain in their source
            ! layers. Dynamic peat columns bury all solid pools; dissolved
            ! matter remains available to the separate aqueous operator.
            !$acc parallel loop independent gang default(present)
            do s = 1, ndecomp_pools
               !$acc loop worker vector collapse(2) independent private(c)
               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     if (.not. transport_pool(fc,s)) then
                        c = filter_soilc (fc)
                        conc_trcr(fc,j,s) = &
                             transport_ptr_list(i_type)%conc_ptr(c,j,s) + &
                             transport_ptr_list(i_type)%src_ptr(c,j,s)
                        if (is_microbial(s) .or. &
                             (use_microbe_aqueous_transport .and. is_dissolved(s)) .or. &
                             (peat_dynamic_column(fc) .and. is_dissolved(s))) then
                           transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) = 0._r8
                        end if
                     end if
                  end do
               end do
            end do


            !$acc parallel loop independent gang collapse(2) default(present)
            do s = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  !$acc loop vector independent private(c)
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     transport_ptr_list(i_type)%conc_ptr(c,j,s) = conc_trcr(fc,j,s)
                  end do
               end do

            end do ! s (pool loop)

         else !use_vertsoilc?

            !! for single level case, no transport; just update the fluxes calculated in the StateUpdate1 subroutines
            !$acc parallel loop independent collapse(2) default(present)
            do l = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  !$acc loop vector independent private(c)
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     transport_ptr_list(i_type)%conc_ptr(c,j,l) = transport_ptr_list(i_type)%conc_ptr(c,j,l) &
                                                                  +transport_ptr_list(i_type)%src_ptr(c,j,l)
                     transport_ptr_list(i_type)%trcr_tend_ptr(c,j,l) = 0._r8
                  end do
               end do
            end do

         endif

      end do  ! i_type
   
      !$acc exit data delete(a_tri(:,:,:),b_tri(:,:,:),&
      !$acc     c_tri(:,:,:),r_tri(:,:,:), gam(:), &
      !$acc     conc_trcr(:,:,:), spinup_term, i_type, dom_diffusion_multiplier, &
      !$acc     peat_compaction_timescale_seconds, peat_dynamic_column(:), &
      !$acc     transport_pool(:,:), peat_c_density(:,:))
    end associate

  end subroutine SoilLittVertTransp

  subroutine calc_diffus_advflux(spinup_term, year, som_diffus_coef, som_adv_coef,&
                                   cnscalaravg_col, adv_fluxj, diffusj)
    !$acc routine seq
    real(r8), intent(in)  ::  spinup_term
    integer ,  intent(in) :: year
    real(r8), intent(in)  ::  som_diffus_coef, som_adv_coef, cnscalaravg_col
    real(r8), intent(out) ::  adv_fluxj, diffusj

    real(r8), parameter :: eps = 1d-30
    logical :: use_scaling

    ! Determine if we should scale by cnscalaravg_col
    use_scaling = (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1)

    ! Compute adv_fluxj
    if (abs(som_adv_coef) * spinup_term < eps) then
        adv_fluxj = eps
    else
        adv_fluxj = som_adv_coef * spinup_term
        if (use_scaling) adv_fluxj = adv_fluxj / cnscalaravg_col
    endif

    ! Compute diffusj
    if (abs(som_diffus_coef) * spinup_term < eps) then
        diffusj = eps
    else
        diffusj = som_diffus_coef * spinup_term
        if (use_scaling) diffusj = diffusj / cnscalaravg_col
    endif

 end subroutine calc_diffus_advflux


end module SoilLittVertTranspMod
