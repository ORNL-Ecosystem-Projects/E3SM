module MicrobeMethaneMod

  ! State and lifecycle owner for the revised microbial methane backend.
  ! Phase 3 step 1 adds only parameter, state, history, and restart support.
  ! CH4Mod remains the executing backend until later explicit dispatch work.
  ! DOM, bacteria, and fungi remain authoritative decomposition pools.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod, only : bounds_type
  use elm_varctl, only : use_microbe_methane
  use pio, only : file_desc_t

  implicit none
  private
  save

  type, public :: microbe_methane_type
     ! Carbon-bearing state is g C m-3 soil; gas state is mol gas m-3.
     ! Partition values are concentrations, not area-weighted bulk copies.
     real(r8), pointer :: acetate_c_unsat_col(:,:) => null()
     real(r8), pointer :: acetate_c_sat_col(:,:) => null()
     real(r8), pointer :: acetate_methanogen_c_unsat_col(:,:) => null()
     real(r8), pointer :: acetate_methanogen_c_sat_col(:,:) => null()
     real(r8), pointer :: h2_methanogen_c_unsat_col(:,:) => null()
     real(r8), pointer :: h2_methanogen_c_sat_col(:,:) => null()
     real(r8), pointer :: aerobic_methanotroph_c_unsat_col(:,:) => null()
     real(r8), pointer :: aerobic_methanotroph_c_sat_col(:,:) => null()
     real(r8), pointer :: anaerobic_methanotroph_c_unsat_col(:,:) => null()
     real(r8), pointer :: anaerobic_methanotroph_c_sat_col(:,:) => null()
     real(r8), pointer :: conc_ch4_unsat_col(:,:) => null()
     real(r8), pointer :: conc_ch4_sat_col(:,:) => null()
     real(r8), pointer :: conc_o2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_o2_sat_col(:,:) => null()
     real(r8), pointer :: conc_co2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_co2_sat_col(:,:) => null()
     real(r8), pointer :: conc_h2_unsat_col(:,:) => null()
     real(r8), pointer :: conc_h2_sat_col(:,:) => null()

     ! Future conservative partition remapping uses this prior-step fraction.
     real(r8), pointer :: sat_fraction_previous_col(:) => null()
   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitCold
     procedure, private :: InitHistory
     procedure, private :: ReadParams
     procedure, public  :: Restart
  end type microbe_methane_type

contains

  subroutine Init(this, bounds)
    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds

    if (.not. use_microbe_methane) return
    call this%ReadParams()
    call this%InitAllocate(bounds)
    call this%InitCold(bounds)
    call this%InitHistory(bounds)
  end subroutine Init

  subroutine InitAllocate(this, bounds)
    use elm_varpar, only : nlevdecomp_full
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer :: begc, endc

    begc = bounds%begc
    endc = bounds%endc
    allocate(this%acetate_c_unsat_col(begc:endc,1:nlevdecomp_full)); this%acetate_c_unsat_col = nan
    allocate(this%acetate_c_sat_col(begc:endc,1:nlevdecomp_full)); this%acetate_c_sat_col = nan
    allocate(this%acetate_methanogen_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%acetate_methanogen_c_unsat_col = nan
    allocate(this%acetate_methanogen_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%acetate_methanogen_c_sat_col = nan
    allocate(this%h2_methanogen_c_unsat_col(begc:endc,1:nlevdecomp_full)); this%h2_methanogen_c_unsat_col = nan
    allocate(this%h2_methanogen_c_sat_col(begc:endc,1:nlevdecomp_full)); this%h2_methanogen_c_sat_col = nan
    allocate(this%aerobic_methanotroph_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%aerobic_methanotroph_c_unsat_col = nan
    allocate(this%aerobic_methanotroph_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%aerobic_methanotroph_c_sat_col = nan
    allocate(this%anaerobic_methanotroph_c_unsat_col(begc:endc,1:nlevdecomp_full)); &
         this%anaerobic_methanotroph_c_unsat_col = nan
    allocate(this%anaerobic_methanotroph_c_sat_col(begc:endc,1:nlevdecomp_full)); &
         this%anaerobic_methanotroph_c_sat_col = nan
    allocate(this%conc_ch4_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_ch4_unsat_col = nan
    allocate(this%conc_ch4_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_ch4_sat_col = nan
    allocate(this%conc_o2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_o2_unsat_col = nan
    allocate(this%conc_o2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_o2_sat_col = nan
    allocate(this%conc_co2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_co2_unsat_col = nan
    allocate(this%conc_co2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_co2_sat_col = nan
    allocate(this%conc_h2_unsat_col(begc:endc,1:nlevdecomp_full)); this%conc_h2_unsat_col = nan
    allocate(this%conc_h2_sat_col(begc:endc,1:nlevdecomp_full)); this%conc_h2_sat_col = nan
    allocate(this%sat_fraction_previous_col(begc:endc)); this%sat_fraction_previous_col = nan
  end subroutine InitAllocate

  subroutine ReadParams(this)
    use MicrobeMethaneParamsMod, only : assertMicrobeMethaneParamsRead
    class(microbe_methane_type) :: this
    call assertMicrobeMethaneParamsRead()
  end subroutine ReadParams

  subroutine InitCold(this, bounds)
    use elm_varpar, only : nlevdecomp
    use shr_infnan_mod, only : spval => shr_infnan_nan, assignment(=)
    use ColumnType, only : col_pp
    use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsInst

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer :: c
    real(r8) :: biomass_seed

    biomass_seed = MicrobeMethaneParamsInst%mfg_biomass_min
    this%acetate_c_unsat_col = spval
    this%acetate_c_sat_col = spval
    this%acetate_methanogen_c_unsat_col = spval
    this%acetate_methanogen_c_sat_col = spval
    this%h2_methanogen_c_unsat_col = spval
    this%h2_methanogen_c_sat_col = spval
    this%aerobic_methanotroph_c_unsat_col = spval
    this%aerobic_methanotroph_c_sat_col = spval
    this%anaerobic_methanotroph_c_unsat_col = spval
    this%anaerobic_methanotroph_c_sat_col = spval
    this%conc_ch4_unsat_col = spval
    this%conc_ch4_sat_col = spval
    this%conc_o2_unsat_col = spval
    this%conc_o2_sat_col = spval
    this%conc_co2_unsat_col = spval
    this%conc_co2_sat_col = spval
    this%conc_h2_unsat_col = spval
    this%conc_h2_sat_col = spval
    this%sat_fraction_previous_col = spval

    do c = bounds%begc, bounds%endc
       if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
          this%acetate_c_unsat_col(c,1:nlevdecomp) = 0._r8
          this%acetate_c_sat_col(c,1:nlevdecomp) = 0._r8
          this%acetate_methanogen_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%acetate_methanogen_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%h2_methanogen_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%h2_methanogen_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%aerobic_methanotroph_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%aerobic_methanotroph_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%anaerobic_methanotroph_c_unsat_col(c,1:nlevdecomp) = biomass_seed
          this%anaerobic_methanotroph_c_sat_col(c,1:nlevdecomp) = biomass_seed
          this%conc_ch4_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_ch4_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_o2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_o2_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_co2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_co2_sat_col(c,1:nlevdecomp) = 0._r8
          this%conc_h2_unsat_col(c,1:nlevdecomp) = 0._r8
          this%conc_h2_sat_col(c,1:nlevdecomp) = 0._r8
          this%sat_fraction_previous_col(c) = 0._r8
       end if
    end do
  end subroutine InitCold

  subroutine InitHistory(this, bounds)
    use histFileMod, only : hist_addfld_decomp

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds

    call add_state('MM_ACETATE_C_UNSAT', 'gC/m^3', &
         'acetate carbon concentration in unsaturated subarea', this%acetate_c_unsat_col)
    call add_state('MM_ACETATE_C_SAT', 'gC/m^3', &
         'acetate carbon concentration in saturated subarea', this%acetate_c_sat_col)
    call add_state('MM_ACET_METH_C_UNSAT', 'gC/m^3', &
         'acetoclastic methanogen biomass carbon in unsaturated subarea', &
         this%acetate_methanogen_c_unsat_col)
    call add_state('MM_ACET_METH_C_SAT', 'gC/m^3', &
         'acetoclastic methanogen biomass carbon in saturated subarea', &
         this%acetate_methanogen_c_sat_col)
    call add_state('MM_H2_METH_C_UNSAT', 'gC/m^3', &
         'hydrogenotrophic methanogen biomass carbon in unsaturated subarea', &
         this%h2_methanogen_c_unsat_col)
    call add_state('MM_H2_METH_C_SAT', 'gC/m^3', &
         'hydrogenotrophic methanogen biomass carbon in saturated subarea', &
         this%h2_methanogen_c_sat_col)
    call add_state('MM_AER_METHANOTROPH_C_UNSAT', 'gC/m^3', &
         'aerobic methanotroph biomass carbon in unsaturated subarea', &
         this%aerobic_methanotroph_c_unsat_col)
    call add_state('MM_AER_METHANOTROPH_C_SAT', 'gC/m^3', &
         'aerobic methanotroph biomass carbon in saturated subarea', &
         this%aerobic_methanotroph_c_sat_col)
    call add_state('MM_ANAER_METHANOTROPH_C_UNSAT', 'gC/m^3', &
         'anaerobic methanotroph biomass carbon in unsaturated subarea', &
         this%anaerobic_methanotroph_c_unsat_col)
    call add_state('MM_ANAER_METHANOTROPH_C_SAT', 'gC/m^3', &
         'anaerobic methanotroph biomass carbon in saturated subarea', &
         this%anaerobic_methanotroph_c_sat_col)
    call add_state('MM_CONC_CH4_UNSAT', 'mol/m^3', &
         'dissolved methane concentration in unsaturated subarea', this%conc_ch4_unsat_col)
    call add_state('MM_CONC_CH4_SAT', 'mol/m^3', &
         'dissolved methane concentration in saturated subarea', this%conc_ch4_sat_col)
    call add_state('MM_CONC_O2_UNSAT', 'mol/m^3', &
         'dissolved oxygen concentration in unsaturated subarea', this%conc_o2_unsat_col)
    call add_state('MM_CONC_O2_SAT', 'mol/m^3', &
         'dissolved oxygen concentration in saturated subarea', this%conc_o2_sat_col)
    call add_state('MM_CONC_CO2_UNSAT', 'mol/m^3', &
         'dissolved carbon dioxide concentration in unsaturated subarea', this%conc_co2_unsat_col)
    call add_state('MM_CONC_CO2_SAT', 'mol/m^3', &
         'dissolved carbon dioxide concentration in saturated subarea', this%conc_co2_sat_col)
    call add_state('MM_CONC_H2_UNSAT', 'mol/m^3', &
         'dissolved hydrogen concentration in unsaturated subarea', this%conc_h2_unsat_col)
    call add_state('MM_CONC_H2_SAT', 'mol/m^3', &
         'dissolved hydrogen concentration in saturated subarea', this%conc_h2_sat_col)

  contains
    subroutine add_state(name, units, long_name, field)
      character(len=*), intent(in) :: name, units, long_name
      real(r8), pointer, intent(inout) :: field(:,:)
      call hist_addfld_decomp(fname=name, units=units, type2d='levdcmp', &
           avgflag='A', long_name=long_name, ptr_col=field, default='inactive')
    end subroutine add_state
  end subroutine InitHistory

  subroutine Restart(this, bounds, ncid, flag)
    use ncdio_pio, only : ncd_double
    use restUtilMod, only : restartvar
    use elm_varctl, only : nsrest, nsrStartup
    use abortutils, only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg

    class(microbe_methane_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in) :: flag
    logical :: readvar

    if (.not. use_microbe_methane) return
    call restart_state('MM_ACETATE_C_UNSAT', 'gC/m^3', this%acetate_c_unsat_col)
    call restart_state('MM_ACETATE_C_SAT', 'gC/m^3', this%acetate_c_sat_col)
    call restart_state('MM_ACET_METH_C_UNSAT', 'gC/m^3', this%acetate_methanogen_c_unsat_col)
    call restart_state('MM_ACET_METH_C_SAT', 'gC/m^3', this%acetate_methanogen_c_sat_col)
    call restart_state('MM_H2_METH_C_UNSAT', 'gC/m^3', this%h2_methanogen_c_unsat_col)
    call restart_state('MM_H2_METH_C_SAT', 'gC/m^3', this%h2_methanogen_c_sat_col)
    call restart_state('MM_AER_METHANOTROPH_C_UNSAT', 'gC/m^3', this%aerobic_methanotroph_c_unsat_col)
    call restart_state('MM_AER_METHANOTROPH_C_SAT', 'gC/m^3', this%aerobic_methanotroph_c_sat_col)
    call restart_state('MM_ANAER_METHANOTROPH_C_UNSAT', 'gC/m^3', this%anaerobic_methanotroph_c_unsat_col)
    call restart_state('MM_ANAER_METHANOTROPH_C_SAT', 'gC/m^3', this%anaerobic_methanotroph_c_sat_col)
    call restart_state('MM_CONC_CH4_UNSAT', 'mol/m^3', this%conc_ch4_unsat_col)
    call restart_state('MM_CONC_CH4_SAT', 'mol/m^3', this%conc_ch4_sat_col)
    call restart_state('MM_CONC_O2_UNSAT', 'mol/m^3', this%conc_o2_unsat_col)
    call restart_state('MM_CONC_O2_SAT', 'mol/m^3', this%conc_o2_sat_col)
    call restart_state('MM_CONC_CO2_UNSAT', 'mol/m^3', this%conc_co2_unsat_col)
    call restart_state('MM_CONC_CO2_SAT', 'mol/m^3', this%conc_co2_sat_col)
    call restart_state('MM_CONC_H2_UNSAT', 'mol/m^3', this%conc_h2_unsat_col)
    call restart_state('MM_CONC_H2_SAT', 'mol/m^3', this%conc_h2_sat_col)
    call restartvar(ncid=ncid, flag=flag, varname='MM_SAT_FRACTION_PREVIOUS', xtype=ncd_double, &
         dim1name='column', long_name='previous revised-methane saturated-area fraction', units='1', &
         readvar=readvar, interpinic_flag='interp', data=this%sat_fraction_previous_col)
    call require_state_on_restart('MM_SAT_FRACTION_PREVIOUS', readvar)

  contains
    subroutine restart_state(name, units, field)
      character(len=*), intent(in) :: name, units
      real(r8), pointer, intent(inout) :: field(:,:)
      call restartvar(ncid=ncid, flag=flag, varname=name, xtype=ncd_double, &
           dim1name='column', dim2name='levgrnd', switchdim=.true., &
           long_name='revised microbial methane prognostic state', units=units, &
           readvar=readvar, interpinic_flag='interp', data=field)
      call require_state_on_restart(name, readvar)
    end subroutine restart_state

    subroutine require_state_on_restart(name, was_read)
      character(len=*), intent(in) :: name
      logical, intent(in) :: was_read

      if (flag == 'read' .and. nsrest /= nsrStartup .and. .not. was_read) then
         call endrun(msg=' ERROR: revised methane restart is missing required state: '//trim(name)//&
              '. A legacy restart cannot be continued or branched as revised methane.'//&
              errMsg(__FILE__, __LINE__))
      end if
    end subroutine require_state_on_restart
  end subroutine Restart

end module MicrobeMethaneMod
