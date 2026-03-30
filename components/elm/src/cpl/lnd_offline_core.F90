module lnd_offline_core

  ! Reusable coupler-free ELM init/step/finalize hooks for offline runs.
  ! The first implementation is intentionally narrow: it assumes CPL_BYPASS
  ! forcing and allocates a dummy x2l buffer so the existing lnd_import path
  ! can still be reused without the MCT driver.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use decompMod   , only : bounds_type, get_proc_bounds
  use controlMod  , only : control_setNL
  use elm_cpl_indices, only : elm_cpl_indices_set_offline, nflds_x2l
  use elm_cpl_indices, only : index_x2l_Sa_co2prog, index_x2l_Sa_co2diag
  use elm_driver  , only : elm_drv
  use elm_finalizeMod, only : final
  use elm_initializeMod, only : initialize1, initialize2, initialize3
  use elm_instance , only : elm_instance_set_offline
  use elm_instMod , only : atm2lnd_vars, glc2lnd_vars, lnd2atm_vars
  use elm_instMod , only : iac2lnd_vars, ocn2lnd_vars
  use elm_time_manager, only : get_curr_calday, get_curr_date, get_nstep, get_step_size
  use elm_time_manager, only : set_nextsw_cday, set_timemgr_init, set_timemgr_quiet, update_rad_dtime
  use elm_time_manager, only : advance_timestep
  use elm_varctl  , only : inst_index, inst_name, inst_suffix
  use elm_varctl  , only : fatmlndfrc, iulog
  use elm_varctl  , only : noland, nsrStartup, use_lnd_rof_two_way, use_ocn_lnd_one_way
  use elm_varctl  , only : elm_varctl_set, elm_varctl_set_iac_flag, elm_varctl_set_offline_driver_mode
  use elm_varctl  , only : co2_ppmv
  use elm_varorb  , only : eccen, lambm0, mvelpp, obliqr
  use lnd_import_export, only : lnd_import
  use mpi         , only : mpi_comm_world, mpi_comm_rank, mpi_success
  use spmdMod     , only : masterproc, spmd_init
  use shr_file_mod, only : shr_file_getLogLevel, shr_file_getLogUnit
  use shr_file_mod, only : shr_file_getUnit, shr_file_setIO, shr_file_setLogLevel, shr_file_setLogUnit
  use shr_flds_mod , only : shr_flds_dom_coord, shr_flds_dom_other
  use shr_pio_mod , only : shr_pio_init1, shr_pio_init2, shr_pio_finalize
  use shr_orb_mod , only : shr_orb_decl

  implicit none
  private

  public :: lnd_offline_options_type
  public :: lnd_offline_runtime_type
  public :: lnd_offline_init
  public :: lnd_offline_step
  public :: lnd_offline_finalize
  public :: lnd_offline_is_last_step

  type :: lnd_offline_options_type
     character(len=256) :: nl_filename = 'lnd_in'
     character(len=256) :: pio_nl_filename = 'drv_in'
     character(len=256) :: caseid = 'offline_elm'
     character(len=256) :: ctitle = 'Offline ELM run'
     character(len=256) :: version = 'offline'
     character(len=256) :: hostname = 'offline'
     character(len=256) :: username = 'offline'
     character(len=64)  :: calendar = 'NO_LEAP'
     integer            :: start_ymd = 20000101
     integer            :: start_tod = 0
     integer            :: ref_ymd = 20000101
     integer            :: ref_tod = 0
     integer            :: stop_ymd = 20000102
     integer            :: stop_tod = 0
     integer            :: dtime = -1
     integer            :: sync_dtime = -1
     integer            :: nsrest = nsrStartup
     logical            :: single_column = .false.
     logical            :: scm_multcols = .false.
     real(r8)           :: scmlat = 0._r8
     real(r8)           :: scmlon = 0._r8
     integer            :: scm_nx = 1
     integer            :: scm_ny = 1
     logical            :: atm_present = .false.
     logical            :: iac_present = .false.
     logical            :: use_lnd_rof_two_way = .false.
     logical            :: use_ocn_lnd_one_way = .false.
     real(r8)           :: nextsw_cday = 1._r8
     real(r8)           :: orb_eccen = 0._r8
     real(r8)           :: orb_mvelpp = 0._r8
     real(r8)           :: orb_lambm0 = 0._r8
     real(r8)           :: orb_obliqr = 0._r8
  end type lnd_offline_options_type

  type :: lnd_offline_runtime_type
     type(bounds_type) :: bounds
     logical           :: initialized = .false.
     logical           :: atm_present = .false.
     logical           :: pio_initialized = .false.
     logical           :: using_logfile = .false.
     real(r8)          :: nextsw_cday = 1._r8
     integer           :: shrlogunit = 6
     integer           :: shrloglev = 1
     integer           :: import_nstep_interval = 1
     real(r8), allocatable :: x2l_dummy(:,:)
  end type lnd_offline_runtime_type

contains

  subroutine lnd_offline_init(options, runtime)

    type(lnd_offline_options_type), intent(in)    :: options
    type(lnd_offline_runtime_type), intent(inout) :: runtime

    integer :: thisng
    integer :: global_comm
    integer :: comp_comm(1)
    integer :: comp_comm_iam(1)
    integer :: ierr
    integer :: comp_id(1)
    logical :: comp_iamin(1)
    logical :: exists
    character(len=3) :: comp_name(1)

#ifndef CPL_BYPASS
    call endrun('lnd_offline_init requires CPL_BYPASS so forcing can be read without the driver')
#endif

    call control_setNL(trim(options%nl_filename))

    inst_name   = 'lnd'
    inst_index  = 1
    inst_suffix = ''
    call elm_instance_set_offline(in_lnd_id=1, in_inst_name=inst_name, in_inst_index=inst_index, &
         in_inst_suffix=inst_suffix)

    call elm_varctl_set_iac_flag(options%iac_present)
    call elm_varctl_set_offline_driver_mode(.true.)
    call elm_cpl_indices_set_offline()

    shr_flds_dom_coord = 'lat:lon:hgt'
    shr_flds_dom_other = 'area:aream:mask:frac'

    global_comm = mpi_comm_world
    call spmd_init(global_comm, 1)
    call shr_file_getLogUnit(runtime%shrlogunit)
    call shr_file_getLogLevel(runtime%shrloglev)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix), exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix), iulog)
          runtime%using_logfile = .true.
       else
          iulog = runtime%shrlogunit
       end if
    else
       iulog = runtime%shrlogunit
    end if
    call shr_file_setLogUnit(iulog)

    call shr_pio_init1(1, trim(options%pio_nl_filename), global_comm)
    comp_id(1) = 1
    comp_name(1) = 'LND'
    comp_iamin(1) = .true.
    comp_comm(1) = global_comm
    call mpi_comm_rank(global_comm, comp_comm_iam(1), ierr)
    if (ierr /= mpi_success) then
       call endrun('lnd_offline_init failed to query MPI rank for offline PIO setup')
    end if
    call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
    runtime%pio_initialized = .true.

    eccen   = options%orb_eccen
    mvelpp  = options%orb_mvelpp
    lambm0  = options%orb_lambm0
    obliqr  = options%orb_obliqr

    use_lnd_rof_two_way = options%use_lnd_rof_two_way
    use_ocn_lnd_one_way = options%use_ocn_lnd_one_way

    if (options%dtime > 0) then
       call set_timemgr_init(calendar_in=trim(options%calendar), start_ymd_in=options%start_ymd, &
            start_tod_in=options%start_tod, ref_ymd_in=options%ref_ymd, ref_tod_in=options%ref_tod, &
            stop_ymd_in=options%stop_ymd, stop_tod_in=options%stop_tod, dtime_in=options%dtime)
    else
       call set_timemgr_init(calendar_in=trim(options%calendar), start_ymd_in=options%start_ymd, &
            start_tod_in=options%start_tod, ref_ymd_in=options%ref_ymd, ref_tod_in=options%ref_tod, &
            stop_ymd_in=options%stop_ymd, stop_tod_in=options%stop_tod)
    end if
    call set_timemgr_quiet(.true.)

    call elm_varctl_set(caseid_in=options%caseid, ctitle_in=options%ctitle, &
         single_column_in=options%single_column, scm_multcols_in=options%scm_multcols, &
         scmlat_in=options%scmlat, scmlon_in=options%scmlon, scm_nx_in=options%scm_nx, &
         scm_ny_in=options%scm_ny, nsrest_in=options%nsrest, version_in=options%version, &
         hostname_in=options%hostname, username_in=options%username)

    call initialize1()

    if (noland) then
       call endrun('lnd_offline_init found noland=.true.; the offline core expects an active land domain')
    end if

    call initialize2()
    call initialize3()

    call get_proc_bounds(runtime%bounds)

    thisng = runtime%bounds%endg - runtime%bounds%begg + 1
    allocate(runtime%x2l_dummy(nflds_x2l, thisng))
    runtime%x2l_dummy(:,:) = 0._r8
    if (index_x2l_Sa_co2prog > 0) runtime%x2l_dummy(index_x2l_Sa_co2prog,:) = co2_ppmv
    if (index_x2l_Sa_co2diag > 0) runtime%x2l_dummy(index_x2l_Sa_co2diag,:) = co2_ppmv

    runtime%atm_present = options%atm_present
    runtime%nextsw_cday = options%nextsw_cday
    call set_nextsw_cday(runtime%nextsw_cday)

    if (.not. runtime%atm_present) then
       runtime%nextsw_cday = calc_nextsw_cday_from_nstep()
       call set_nextsw_cday(runtime%nextsw_cday)
    end if

    if (options%sync_dtime > 0 .and. get_step_size() > 0) then
       runtime%import_nstep_interval = max(1, options%sync_dtime / get_step_size())
    else
       runtime%import_nstep_interval = 1
    end if

    runtime%initialized = .true.

  end subroutine lnd_offline_init

  subroutine lnd_offline_step(runtime, rstwr, nlend, rdate)

    type(lnd_offline_runtime_type), intent(inout) :: runtime
    logical,                        intent(in)    :: rstwr
    logical,                        intent(in)    :: nlend
    character(len=*), optional,     intent(in)    :: rdate

    integer  :: nstep
    integer  :: dtime
    integer  :: yr
    integer  :: mon
    integer  :: day
    integer  :: tod
    logical  :: doalb
    real(r8) :: calday
    real(r8) :: caldayp1
    real(r8) :: declin
    real(r8) :: declinp1
    real(r8) :: eccf
    character(len=32) :: local_rdate

    if (.not. runtime%initialized) then
       call endrun('lnd_offline_step called before lnd_offline_init')
    end if

    nstep = get_nstep()
    if (nstep /= 1 .and. mod(nstep, runtime%import_nstep_interval) == 0) then
       call lnd_import(runtime%bounds, runtime%x2l_dummy, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, &
            lnd2atm_vars, iac2lnd_vars)
    end if

    if (.not. runtime%atm_present) then
       runtime%nextsw_cday = calc_nextsw_cday_from_nstep()
       call set_nextsw_cday(runtime%nextsw_cday)
    end if

    dtime = get_step_size()
    caldayp1 = get_curr_calday(offset=dtime)

    if (nstep == 0) then
       doalb = .false.
    else if (nstep == 1) then
       doalb = (abs(runtime%nextsw_cday - caldayp1) < 1.e-10_r8)
    else
       doalb = (runtime%nextsw_cday >= -0.5_r8)
    end if
    call update_rad_dtime(doalb)

    calday = get_curr_calday()
    call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
    call shr_orb_decl(runtime%nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf)

    if (present(rdate)) then
       local_rdate = rdate
    else
       call get_curr_date(yr, mon, day, tod)
       write(local_rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr, mon, day, tod
    end if

    call elm_drv(doalb, runtime%nextsw_cday, declinp1, declin, rstwr, nlend, local_rdate)
    call advance_timestep()

  end subroutine lnd_offline_step

  subroutine lnd_offline_finalize(runtime)

    type(lnd_offline_runtime_type), intent(inout) :: runtime

    if (allocated(runtime%x2l_dummy)) then
       deallocate(runtime%x2l_dummy)
    end if

    runtime%initialized = .false.

    if (runtime%pio_initialized) then
       call shr_pio_finalize()
       runtime%pio_initialized = .false.
    end if

    call final()

    call shr_file_setLogUnit(runtime%shrlogunit)
    call shr_file_setLogLevel(runtime%shrloglev)
    if (runtime%using_logfile) then
       close(iulog)
       iulog = runtime%shrlogunit
       runtime%using_logfile = .false.
    end if

  end subroutine lnd_offline_finalize

  logical function lnd_offline_is_last_step()

    use elm_time_manager, only : is_last_step

    lnd_offline_is_last_step = is_last_step()

  end function lnd_offline_is_last_step

  real(r8) function calc_nextsw_cday_from_nstep()

    integer  :: dtime
    integer  :: nstep

    dtime = get_step_size()
    nstep = get_nstep()
    calc_nextsw_cday_from_nstep = mod((nstep / (86400._r8 / dtime)) * 1.0_r8, 365._r8) + 1._r8

  end function calc_nextsw_cday_from_nstep

end module lnd_offline_core
