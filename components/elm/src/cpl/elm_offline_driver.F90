program elm_offline_driver

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use abortutils       , only : endrun
  use perf_mod         , only : t_adj_detailf, t_finalizef, t_initf, t_prf, t_startf, t_stopf
  use shr_orb_mod      , only : shr_orb_params
  use ESMF             , only : ESMF_Initialize, ESMF_Finalize
  use ESMF             , only : ESMF_LogKind_Flag, ESMF_LOGKIND_NONE
  use mpi
  use lnd_offline_core , only : lnd_offline_options_type, lnd_offline_runtime_type
  use lnd_offline_core , only : lnd_offline_finalize, lnd_offline_init, lnd_offline_is_last_step
  use lnd_offline_core , only : lnd_offline_step

  implicit none

  type(lnd_offline_options_type) :: options
  type(lnd_offline_runtime_type) :: runtime
  character(len=256)             :: driver_nml
  logical                        :: nlend
  integer                        :: ios
  integer                        :: ierr
  logical                        :: mpi_is_initialized
  logical                        :: mpi_is_finalized
  logical                        :: am_timer_master
  integer                        :: esmf_rc
  type(ESMF_LogKind_Flag)        :: esmf_logfile_kind
  integer                        :: drv_timemgr_find_status
  integer                        :: drv_timemgr_read_status
  character(len=256)             :: drv_timing_dir
  character(len=64)              :: drv_stop_option
  integer                        :: drv_stop_n
  integer                        :: drv_dtime
  integer                        :: drv_orb_iyear
  integer                        :: drv_orb_iyear_align
  character(len=256)             :: nl_filename
  character(len=256)             :: pio_nl_filename
  character(len=256)             :: caseid
  character(len=256)             :: ctitle
  character(len=256)             :: version
  character(len=256)             :: hostname
  character(len=256)             :: username
  character(len=64)              :: calendar
  character(len=64)              :: drv_orb_mode
  integer                        :: start_ymd
  integer                        :: start_tod
  integer                        :: ref_ymd
  integer                        :: ref_tod
  integer                        :: stop_ymd
  integer                        :: stop_tod
  integer                        :: dtime
  integer                        :: sync_dtime
  integer                        :: nsrest
  logical                        :: single_column
  logical                        :: scm_multcols
  real(r8)                       :: scmlat
  real(r8)                       :: scmlon
  integer                        :: scm_nx
  integer                        :: scm_ny
  logical                        :: atm_present
  logical                        :: iac_present
  logical                        :: use_lnd_rof_two_way
  logical                        :: use_ocn_lnd_one_way
  real(r8)                       :: nextsw_cday
  real(r8)                       :: orb_eccen
  real(r8)                       :: orb_obliq
  real(r8)                       :: orb_mvelp
  real(r8)                       :: orb_mvelpp
  real(r8)                       :: orb_lambm0
  real(r8)                       :: orb_obliqr
  logical                        :: have_timer_root

  namelist /offline_elm_nml/ nl_filename, pio_nl_filename, caseid, ctitle, version, hostname, username, calendar, &
       start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, dtime, sync_dtime, nsrest, single_column, &
       scm_multcols, scmlat, scmlon, scm_nx, scm_ny, atm_present, iac_present, &
       use_lnd_rof_two_way, use_ocn_lnd_one_way, nextsw_cday, orb_eccen, orb_mvelpp, &
       orb_lambm0, orb_obliqr

  call MPI_Initialized(mpi_is_initialized, ierr)
  if (ierr /= MPI_SUCCESS) then
     call endrun('elm_offline_driver failed querying MPI initialization state')
  end if
  if (.not. mpi_is_initialized) then
     call MPI_Init(ierr)
     if (ierr /= MPI_SUCCESS) then
        call endrun('elm_offline_driver failed to initialize MPI')
     end if
  end if

  drv_timing_dir = './timing'
  have_timer_root = .false.
  am_timer_master = .true.

  esmf_logfile_kind = ESMF_LOGKIND_NONE
  call ESMF_Initialize(logkindflag=esmf_logfile_kind, rc=esmf_rc)
  if (esmf_rc /= 0) then
     call endrun('elm_offline_driver failed to initialize ESMF')
  end if

  driver_nml = 'elm_offline_driver.nml'
  nl_filename = options%nl_filename
  pio_nl_filename = options%pio_nl_filename
  caseid = options%caseid
  ctitle = options%ctitle
  version = options%version
  hostname = options%hostname
  username = options%username
  calendar = options%calendar
  start_ymd = options%start_ymd
  start_tod = options%start_tod
  ref_ymd = options%ref_ymd
  ref_tod = options%ref_tod
  stop_ymd = options%stop_ymd
  stop_tod = options%stop_tod
  dtime = options%dtime
  sync_dtime = options%sync_dtime
  nsrest = options%nsrest
  single_column = options%single_column
  scm_multcols = options%scm_multcols
  scmlat = options%scmlat
  scmlon = options%scmlon
  scm_nx = options%scm_nx
  scm_ny = options%scm_ny
  atm_present = options%atm_present
  iac_present = options%iac_present
  use_lnd_rof_two_way = options%use_lnd_rof_two_way
  use_ocn_lnd_one_way = options%use_ocn_lnd_one_way
  nextsw_cday = options%nextsw_cday
  orb_eccen = options%orb_eccen
  orb_mvelpp = options%orb_mvelpp
  orb_lambm0 = options%orb_lambm0
  orb_obliqr = options%orb_obliqr
  if (command_argument_count() >= 1) then
     call get_command_argument(1, driver_nml)
  end if

  open(unit=10, file=trim(driver_nml), status='old', action='read', iostat=ios)
  if (ios == 0) then
     read(10, nml=offline_elm_nml, iostat=ios)
     close(10)
     if (ios /= 0) then
        call endrun('elm_offline_driver failed reading offline_elm_nml from '//trim(driver_nml))
     end if
  end if

  if (trim(caseid) == 'offline_elm') then
     call infer_caseid_from_pwd(caseid)
  end if
  if (trim(ctitle) == 'Offline ELM run') then
     ctitle = caseid
  end if

  call read_drv_timing_dir(trim(pio_nl_filename), drv_timing_dir)
  call ensure_directory_exists(trim(drv_timing_dir))
  call t_initf(trim(pio_nl_filename), LogPrint=.true., mpicom=MPI_COMM_WORLD, MasterTask=am_timer_master, MaxThreads=1)
  call t_startf('LND:OFFLINE')
  call t_adj_detailf(+1)
  have_timer_root = .true.

  call read_drv_timemgr(trim(pio_nl_filename), calendar, start_ymd, start_tod, stop_ymd, stop_tod, drv_stop_option, &
       drv_stop_n, drv_dtime, drv_timemgr_find_status, drv_timemgr_read_status)
  call validate_drv_timemgr(trim(pio_nl_filename), calendar, start_ymd, start_tod, drv_timemgr_find_status, drv_timemgr_read_status)
  ref_ymd = start_ymd
  ref_tod = start_tod
  if (stop_ymd <= 0) then
     call derive_stop_from_option(start_ymd, start_tod, calendar, trim(drv_stop_option), drv_stop_n, stop_ymd, stop_tod)
  end if
  if (drv_dtime > 0) then
     sync_dtime = drv_dtime
  end if
  drv_orb_mode = 'fixed_year'
  drv_orb_iyear = start_ymd_to_year(start_ymd)
  drv_orb_iyear_align = drv_orb_iyear
  call read_drv_orbit(trim(pio_nl_filename), drv_orb_mode, drv_orb_iyear, drv_orb_iyear_align, orb_eccen, orb_obliq, orb_mvelp)
  call initialize_orbit_from_drv(start_ymd, drv_orb_mode, drv_orb_iyear, drv_orb_iyear_align, &
       orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp)

  options%nl_filename = nl_filename
  options%pio_nl_filename = pio_nl_filename
  options%caseid = caseid
  options%ctitle = ctitle
  options%version = version
  options%hostname = hostname
  options%username = username
  options%calendar = calendar
  options%start_ymd = start_ymd
  options%start_tod = start_tod
  options%ref_ymd = ref_ymd
  options%ref_tod = ref_tod
  options%stop_ymd = stop_ymd
  options%stop_tod = stop_tod
  options%dtime = dtime
  options%sync_dtime = sync_dtime
  options%nsrest = nsrest
  options%single_column = single_column
  options%scm_multcols = scm_multcols
  options%scmlat = scmlat
  options%scmlon = scmlon
  options%scm_nx = scm_nx
  options%scm_ny = scm_ny
  options%atm_present = atm_present
  options%iac_present = iac_present
  options%use_lnd_rof_two_way = use_lnd_rof_two_way
  options%use_ocn_lnd_one_way = use_ocn_lnd_one_way
  options%nextsw_cday = nextsw_cday
  options%orb_eccen = orb_eccen
  options%orb_mvelpp = orb_mvelpp
  options%orb_lambm0 = orb_lambm0
  options%orb_obliqr = orb_obliqr

  call lnd_offline_init(options, runtime)

  do
     nlend = lnd_offline_is_last_step()
     call lnd_offline_step(runtime, rstwr=.false., nlend=nlend)
     if (nlend) exit
  end do

  call lnd_offline_finalize(runtime)

  if (have_timer_root) then
     call t_adj_detailf(-1)
     call t_stopf('LND:OFFLINE')
     call t_prf(trim(drv_timing_dir)//'/model_timing.lnd', mpicom=MPI_COMM_WORLD)
     call t_finalizef()
  end if

  call MPI_Finalized(mpi_is_finalized, ierr)
  if (ierr /= MPI_SUCCESS) then
     call endrun('elm_offline_driver failed querying MPI finalized state')
  end if
  if (.not. mpi_is_finalized) then
     call ESMF_Finalize(rc=esmf_rc)
     if (esmf_rc /= 0) then
        call endrun('elm_offline_driver failed to finalize ESMF')
     end if
     call MPI_Finalize(ierr)
     if (ierr /= MPI_SUCCESS) then
        call endrun('elm_offline_driver failed to finalize MPI')
     end if
   end if

contains

  subroutine read_drv_timemgr(nlfile, calendar, start_ymd, start_tod, stop_ymd, stop_tod, stop_option, stop_n, dtime, &
       find_status, read_status)

    character(len=*), intent(in)    :: nlfile
    character(len=*), intent(inout) :: calendar
    integer         , intent(inout) :: start_ymd
    integer         , intent(inout) :: start_tod
    integer         , intent(inout) :: stop_ymd
    integer         , intent(inout) :: stop_tod
    character(len=*), intent(inout) :: stop_option
    integer         , intent(inout) :: stop_n
    integer         , intent(inout) :: dtime
    integer         , intent(out)   :: find_status
    integer         , intent(out)   :: read_status

    integer :: unitn
    integer :: ios_local
    integer :: parse_status
    integer :: p1, p2
    logical :: exists
    logical :: in_group
    character(len=1024) :: line
    character(len=1024) :: line_lower
    character(len=1024) :: payload

    inquire(file=trim(nlfile), exist=exists)
    find_status = -999999
    read_status = -999999
    if (.not. exists) return
    stop_option = ' '
    stop_n = -1
    dtime = -1

    unitn = 11
    open(unit=unitn, file=trim(nlfile), status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) return
    in_group = .false.
    find_status = 1
    read_status = 0
    do
       read(unitn,'(A)',iostat=ios_local) line
       if (ios_local /= 0) exit
       line_lower = to_lower_ascii(line)
       if (.not. in_group) then
          p1 = index(line_lower, '&seq_timemgr_inparm')
          if (p1 > 0) then
             in_group = .true.
             find_status = 0
             payload = line(p1 + len('&seq_timemgr_inparm'):)
             call parse_drv_timemgr_assignments(payload, calendar, start_ymd, start_tod, stop_ymd, stop_tod, &
                  stop_option, stop_n, dtime, parse_status)
             if (parse_status /= 0 .and. read_status == 0) read_status = parse_status
          end if
       else
          p2 = index(line_lower, '/')
          if (p2 > 0) then
             payload = line(:p2-1)
             call parse_drv_timemgr_assignments(payload, calendar, start_ymd, start_tod, stop_ymd, stop_tod, &
                  stop_option, stop_n, dtime, parse_status)
             if (parse_status /= 0 .and. read_status == 0) read_status = parse_status
             exit
          else
             call parse_drv_timemgr_assignments(line, calendar, start_ymd, start_tod, stop_ymd, stop_tod, &
                  stop_option, stop_n, dtime, parse_status)
             if (parse_status /= 0 .and. read_status == 0) read_status = parse_status
          end if
       end if
    end do
    if (find_status /= 0) read_status = -999999
    close(unitn)

  end subroutine read_drv_timemgr

  subroutine read_drv_timing_dir(nlfile, timing_dir)

    character(len=*), intent(in)    :: nlfile
    character(len=*), intent(inout) :: timing_dir

    integer :: unitn
    integer :: ios_local
    integer :: eq_pos
    logical :: exists
    character(len=1024) :: line
    character(len=1024) :: line_clean
    character(len=1024) :: line_lower

    inquire(file=trim(nlfile), exist=exists)
    if (.not. exists) return

    unitn = 12
    open(unit=unitn, file=trim(nlfile), status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) return

    do
       read(unitn,'(A)',iostat=ios_local) line
       if (ios_local /= 0) exit
       line_clean = strip_comment(line)
       line_lower = to_lower_ascii(line_clean)
       eq_pos = index(line_lower, '=')
       if (eq_pos <= 0) cycle
       if (trim(adjustl(line_lower(:eq_pos-1))) == 'timing_dir') then
          timing_dir = trim(adjustl(strip_quotes(line_clean(eq_pos+1:))))
          exit
       end if
    end do

    close(unitn)

  end subroutine read_drv_timing_dir

  subroutine ensure_directory_exists(path)

    character(len=*), intent(in) :: path

    logical :: exists
    integer :: cmdstat
    character(len=2048) :: cmd

    inquire(file=trim(path), exist=exists)
    if (exists) return

    cmd = 'mkdir -p "' // trim(path) // '"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=cmdstat)
    if (cmdstat /= 0) then
       call endrun('elm_offline_driver failed to create timing_dir='//trim(path))
    end if

  end subroutine ensure_directory_exists

  subroutine parse_drv_timemgr_assignments(text, calendar, start_ymd, start_tod, stop_ymd, stop_tod, stop_option, &
       stop_n, dtime, parse_status)

    character(len=*), intent(in)    :: text
    character(len=*), intent(inout) :: calendar
    integer         , intent(inout) :: start_ymd, start_tod
    integer         , intent(inout) :: stop_ymd, stop_tod
    character(len=*), intent(inout) :: stop_option
    integer         , intent(inout) :: stop_n
    integer         , intent(inout) :: dtime
    integer         , intent(out)   :: parse_status

    integer :: comma_pos
    integer :: eq_pos
    integer :: ios_local
    character(len=1024) :: rest
    character(len=1024) :: token
    character(len=256) :: key
    character(len=256) :: value

    parse_status = 0
    rest = strip_comment(text)
    do while (len_trim(rest) > 0)
       comma_pos = index(rest, ',')
       if (comma_pos > 0) then
          token = rest(:comma_pos-1)
          rest = rest(comma_pos+1:)
       else
          token = rest
          rest = ' '
       end if
       eq_pos = index(token, '=')
       if (eq_pos <= 0) cycle
       key = trim(adjustl(to_lower_ascii(token(:eq_pos-1))))
       value = trim(adjustl(strip_quotes(token(eq_pos+1:))))
       select case (trim(key))
       case ('calendar')
          calendar = trim(value)
       case ('start_ymd')
          read(value,*,iostat=ios_local) start_ymd
          if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
       case ('start_tod')
          read(value,*,iostat=ios_local) start_tod
          if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
       case ('stop_ymd')
          read(value,*,iostat=ios_local) stop_ymd
          if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
       case ('stop_tod')
          read(value,*,iostat=ios_local) stop_tod
          if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
       case ('stop_option')
          stop_option = trim(value)
      case ('stop_n')
         read(value,*,iostat=ios_local) stop_n
         if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
       case ('dtime')
          read(value,*,iostat=ios_local) dtime
          if (ios_local /= 0 .and. parse_status == 0) parse_status = ios_local
      end select
    end do

  end subroutine parse_drv_timemgr_assignments

  subroutine read_drv_orbit(nlfile, orb_mode, orb_iyear, orb_iyear_align, orb_eccen, orb_obliq, orb_mvelp)

    character(len=*), intent(in)    :: nlfile
    character(len=*), intent(inout) :: orb_mode
    integer         , intent(inout) :: orb_iyear
    integer         , intent(inout) :: orb_iyear_align
    real(r8)        , intent(inout) :: orb_eccen
    real(r8)        , intent(inout) :: orb_obliq
    real(r8)        , intent(inout) :: orb_mvelp

    integer :: unitn
    integer :: ios_local
    integer :: eq_pos
    logical :: exists
    character(len=1024) :: line
    character(len=1024) :: line_clean
    character(len=1024) :: key
    character(len=1024) :: value

    inquire(file=trim(nlfile), exist=exists)
    if (.not. exists) return

    unitn = 13
    open(unit=unitn, file=trim(nlfile), status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) return

    do
       read(unitn,'(A)',iostat=ios_local) line
       if (ios_local /= 0) exit
       line_clean = strip_comment(line)
       eq_pos = index(line_clean, '=')
       if (eq_pos <= 0) cycle
       key = trim(adjustl(to_lower_ascii(line_clean(:eq_pos-1))))
       value = trim(adjustl(strip_quotes(line_clean(eq_pos+1:))))
       select case (trim(key))
       case ('orb_mode')
          orb_mode = trim(to_lower_ascii(value))
       case ('orb_iyear')
          read(value,*,iostat=ios_local) orb_iyear
       case ('orb_iyear_align')
          read(value,*,iostat=ios_local) orb_iyear_align
       case ('orb_eccen')
          read(value,*,iostat=ios_local) orb_eccen
       case ('orb_obliq')
          read(value,*,iostat=ios_local) orb_obliq
       case ('orb_mvelp')
          read(value,*,iostat=ios_local) orb_mvelp
       end select
    end do

    close(unitn)

  end subroutine read_drv_orbit

  subroutine initialize_orbit_from_drv(start_ymd, orb_mode, orb_iyear, orb_iyear_align, orb_eccen, orb_obliq, orb_mvelp, &
       orb_obliqr, orb_lambm0, orb_mvelpp)

    integer         , intent(in)    :: start_ymd
    character(len=*), intent(in)    :: orb_mode
    integer         , intent(in)    :: orb_iyear
    integer         , intent(in)    :: orb_iyear_align
    real(r8)        , intent(inout) :: orb_eccen
    real(r8)        , intent(inout) :: orb_obliq
    real(r8)        , intent(inout) :: orb_mvelp
    real(r8)        , intent(out)   :: orb_obliqr
    real(r8)        , intent(out)   :: orb_lambm0
    real(r8)        , intent(out)   :: orb_mvelpp

    integer :: year
    integer :: orb_cyear

    year = start_ymd_to_year(start_ymd)
    orb_cyear = orb_iyear
    if (trim(to_lower_ascii(orb_mode)) == 'variable_year') then
       orb_cyear = orb_iyear + (year - orb_iyear_align)
    end if

    call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp, .false.)

  end subroutine initialize_orbit_from_drv

  integer function start_ymd_to_year(ymd)

    integer, intent(in) :: ymd

    start_ymd_to_year = ymd / 10000

  end function start_ymd_to_year

  function strip_comment(text) result(out)

    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: bang_pos

    bang_pos = index(text, '!')
    if (bang_pos > 0) then
       out = text(:bang_pos-1)
    else
       out = text
    end if

  end function strip_comment

  function strip_quotes(text) result(out)

    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    character(len=len(text)) :: tmp
    integer :: n

    tmp = trim(adjustl(text))
    n = len_trim(tmp)
    if (n >= 2) then
       if ((tmp(1:1) == "'" .and. tmp(n:n) == "'") .or. (tmp(1:1) == '"' .and. tmp(n:n) == '"')) then
          out = tmp(2:n-1)
          return
       end if
    end if
    out = tmp

  end function strip_quotes

  function to_lower_ascii(text) result(out)

    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i
    integer :: ich

    out = text
    do i = 1, len(text)
       ich = iachar(text(i:i))
       if (ich >= iachar('A') .and. ich <= iachar('Z')) then
          out(i:i) = achar(ich + 32)
       else
          out(i:i) = text(i:i)
       end if
    end do

  end function to_lower_ascii

  subroutine infer_caseid_from_pwd(caseid)

    character(len=*), intent(inout) :: caseid

    character(len=1024) :: cwd
    character(len=1024) :: trimmed_cwd
    character(len=1024) :: base_name
    integer :: status
    integer :: env_len

    call get_environment_variable('PWD', cwd, length=env_len, status=status)
    if (status /= 0 .or. env_len <= 0) return

    trimmed_cwd = trim(cwd(:env_len))
    call strip_trailing_slashes(trimmed_cwd)
    if (len_trim(trimmed_cwd) == 0) return

    base_name = path_basename(trimmed_cwd)
    if (trim(base_name) == 'run') then
      trimmed_cwd = path_dirname(trimmed_cwd)
      call strip_trailing_slashes(trimmed_cwd)
      if (len_trim(trimmed_cwd) == 0) return
      base_name = path_basename(trimmed_cwd)
    end if

    if (len_trim(base_name) > 0) then
       caseid = trim(base_name)
    end if

  end subroutine infer_caseid_from_pwd

  subroutine strip_trailing_slashes(path)

    character(len=*), intent(inout) :: path
    integer :: n

    n = len_trim(path)
    do while (n > 1 .and. path(n:n) == '/')
       path(n:n) = ' '
       n = n - 1
    end do

  end subroutine strip_trailing_slashes

  function path_basename(path) result(base)

    character(len=*), intent(in) :: path
    character(len=len(path)) :: base
    integer :: slash_pos
    integer :: n

    n = len_trim(path)
    slash_pos = 0
    do while (n > 0)
       if (path(n:n) == '/') then
          slash_pos = n
          exit
       end if
       n = n - 1
    end do

    if (slash_pos == 0) then
       base = trim(path)
    else
       base = trim(path(slash_pos+1:len_trim(path)))
    end if

  end function path_basename

  function path_dirname(path) result(parent)

    character(len=*), intent(in) :: path
    character(len=len(path)) :: parent
    integer :: slash_pos
    integer :: n

    n = len_trim(path)
    slash_pos = 0
    do while (n > 0)
       if (path(n:n) == '/') then
          slash_pos = n
          exit
       end if
       n = n - 1
    end do

    if (slash_pos <= 1) then
       parent = '/'
    else
       parent = trim(path(:slash_pos-1))
    end if

  end function path_dirname

  subroutine validate_drv_timemgr(nlfile, calendar, start_ymd, start_tod, find_status, read_status)

    character(len=*), intent(in) :: nlfile
    character(len=*), intent(in) :: calendar
    integer         , intent(in) :: start_ymd, start_tod
    integer         , intent(in) :: find_status, read_status

    integer :: mm
    character(len=1024) :: errmsg
    character(len=64) :: c_find_status, c_read_status
    character(len=64) :: c_start_ymd, c_start_tod

    mm = mod(start_ymd / 100, 100)
    if (find_status /= 0 .or. read_status /= 0 .or. len_trim(calendar) == 0 .or. mm < 1 .or. mm > 12) then
       write(c_find_status,'(i0)') find_status
       write(c_read_status,'(i0)') read_status
       write(c_start_ymd,'(i0)') start_ymd
       write(c_start_tod,'(i0)') start_tod
       errmsg = 'elm_offline_driver invalid drv_in seq_timemgr_inparm: ' // &
            'nlfile=' // trim(nlfile) // ' ' // &
            'find_status=' // trim(c_find_status) // ' ' // &
            'read_status=' // trim(c_read_status) // ' ' // &
            'calendar=' // trim(calendar) // ' ' // &
            'start_ymd=' // trim(c_start_ymd) // ' ' // &
            'start_tod=' // trim(c_start_tod)
       call endrun(trim(errmsg))
    end if

  end subroutine validate_drv_timemgr

  subroutine add_one_day(start_ymd, calendar, stop_ymd)

    integer         , intent(in)  :: start_ymd
    character(len=*), intent(in)  :: calendar
    integer         , intent(out) :: stop_ymd

    integer :: year
    integer :: month
    integer :: day
    integer :: ndays
    character(len=len(calendar)) :: cal_lower

    year = start_ymd / 10000
    month = mod(start_ymd / 100, 100)
    day = mod(start_ymd, 100)
    cal_lower = to_lower_ascii(calendar)

    select case (month)
    case (1, 3, 5, 7, 8, 10, 12)
       ndays = 31
    case (4, 6, 9, 11)
       ndays = 30
    case (2)
       ndays = 28
       if (trim(cal_lower) == 'gregorian') then
          if (is_gregorian_leap_year(year)) ndays = 29
       end if
    case default
       call endrun('elm_offline_driver add_one_day received invalid month in start_ymd')
    end select

    day = day + 1
    if (day > ndays) then
       day = 1
       month = month + 1
       if (month > 12) then
          month = 1
          year = year + 1
       end if
    end if

    stop_ymd = year * 10000 + month * 100 + day

  end subroutine add_one_day

  subroutine derive_stop_from_option(start_ymd, start_tod, calendar, stop_option, stop_n, stop_ymd, stop_tod)

    integer         , intent(in)  :: start_ymd
    integer         , intent(in)  :: start_tod
    character(len=*), intent(in)  :: calendar
    character(len=*), intent(in)  :: stop_option
    integer         , intent(in)  :: stop_n
    integer         , intent(out) :: stop_ymd
    integer         , intent(out) :: stop_tod

    integer :: year, month, day
    integer :: next_ymd
    character(len=len(stop_option)) :: opt

    opt = trim(to_lower_ascii(stop_option))
    stop_tod = start_tod

    select case (trim(opt))
    case ('ndays', 'nday')
       call add_days(start_ymd, calendar, max(stop_n,1), stop_ymd)
    case ('nhours', 'nhour')
       stop_ymd = start_ymd
       stop_tod = start_tod + max(stop_n,1) * 3600
       do while (stop_tod >= 86400)
          stop_tod = stop_tod - 86400
          call add_days(stop_ymd, calendar, 1, next_ymd)
          stop_ymd = next_ymd
       end do
    case ('nmonths', 'nmonth')
       call split_ymd(start_ymd, year, month, day)
       month = month + max(stop_n,1)
       do while (month > 12)
          month = month - 12
          year = year + 1
       end do
       stop_ymd = year * 10000 + month * 100 + min(day, days_in_month(year, month, calendar))
    case ('nyears', 'nyear')
       call split_ymd(start_ymd, year, month, day)
       year = year + max(stop_n,1)
       stop_ymd = year * 10000 + month * 100 + min(day, days_in_month(year, month, calendar))
    case default
       call add_one_day(start_ymd, calendar, stop_ymd)
    end select

  end subroutine derive_stop_from_option

  subroutine add_days(start_ymd, calendar, nadd, out_ymd)

    integer         , intent(in)  :: start_ymd
    character(len=*), intent(in)  :: calendar
    integer         , intent(in)  :: nadd
    integer         , intent(out) :: out_ymd

    integer :: i
    integer :: year, month, day
    integer :: ndays

    call split_ymd(start_ymd, year, month, day)
    do i = 1, nadd
       ndays = days_in_month(year, month, calendar)
       day = day + 1
       if (day > ndays) then
          day = 1
          month = month + 1
          if (month > 12) then
             month = 1
             year = year + 1
          end if
       end if
    end do
    out_ymd = year * 10000 + month * 100 + day

  end subroutine add_days

  subroutine split_ymd(ymd, year, month, day)

    integer, intent(in) :: ymd
    integer, intent(out) :: year, month, day

    year = ymd / 10000
    month = mod(ymd / 100, 100)
    day = mod(ymd, 100)

  end subroutine split_ymd

  integer function days_in_month(year, month, calendar)

    integer         , intent(in) :: year, month
    character(len=*), intent(in) :: calendar
    character(len=len(calendar)) :: cal_lower

    cal_lower = to_lower_ascii(calendar)
    select case (month)
    case (1, 3, 5, 7, 8, 10, 12)
       days_in_month = 31
    case (4, 6, 9, 11)
       days_in_month = 30
    case (2)
       days_in_month = 28
       if (trim(cal_lower) == 'gregorian') then
          if (is_gregorian_leap_year(year)) days_in_month = 29
       end if
    case default
       call endrun('elm_offline_driver days_in_month received invalid month')
    end select

  end function days_in_month

  logical function is_gregorian_leap_year(year)

    integer, intent(in) :: year

    is_gregorian_leap_year = .false.
    if (mod(year,4) /= 0) return
    if (mod(year,100) /= 0) then
       is_gregorian_leap_year = .true.
    else if (mod(year,400) == 0) then
       is_gregorian_leap_year = .true.
    end if

  end function is_gregorian_leap_year

end program elm_offline_driver
