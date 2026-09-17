module MicrobeMethaneMod

  !-----------------------------------------------------------------------
  ! Lifecycle boundary for the microbial decomposition and revised methane
  ! backend. Phase 1 intentionally owns no prognostic state and stops during
  ! initialization. Later phases will add state only behind
  ! use_microbe_methane and will read scientific parameters from ELM's
  ! standard parameter file.
  !-----------------------------------------------------------------------

  use abortutils  , only : endrun
  use CH4varcon   , only : allowlakeprod
  use decompMod   , only : bounds_type
  use elm_varctl  , only : use_microbe_methane
  use pio         , only : file_desc_t
  use shr_log_mod , only : errMsg => shr_log_errMsg

  implicit none
  private
  save

  type, public :: microbe_methane_type
   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitCold
     procedure, private :: InitHistory
     procedure, private :: ReadParams
     procedure, public  :: Restart
  end type microbe_methane_type

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    class(microbe_methane_type)   :: this
    type(bounds_type), intent(in) :: bounds

    if (.not. use_microbe_methane) return

    if (allowlakeprod) then
       call endrun(msg=' ERROR: use_microbe_methane=.true. requires allowlakeprod=.false.'//&
            errMsg(__FILE__, __LINE__))
    end if

    call this%InitAllocate(bounds)
    call this%ReadParams()
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

    ! ReadParams is a deliberate Phase 1 stop. Retain this second guard so
    ! later work cannot expose a half-connected timestep path accidentally.
    call endrun(msg=' ERROR: use_microbe_methane is not ready for timestepping.'//&
         errMsg(__FILE__, __LINE__))

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    class(microbe_methane_type)   :: this
    type(bounds_type), intent(in) :: bounds

    ! No state is allocated in the dormant Phase 1 type.

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine ReadParams(this)
    class(microbe_methane_type) :: this

    ! All scientific parameters will be added to and read from ELM's standard
    ! parameter NetCDF in later phases. Do not add a module-specific file.
    call endrun(msg=' ERROR: use_microbe_methane is a Phase 1 development-only option; '//&
         'microbial parameters and timestep coupling are not implemented yet.'//&
         errMsg(__FILE__, __LINE__))

  end subroutine ReadParams

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    class(microbe_methane_type)   :: this
    type(bounds_type), intent(in) :: bounds

    ! History fields are added with prognostic state in a later phase.

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    class(microbe_methane_type)   :: this
    type(bounds_type), intent(in) :: bounds

    ! Cold-start state is added with prognostic state in a later phase.

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    class(microbe_methane_type)     :: this
    type(bounds_type), intent(in)   :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in)    :: flag

    if (use_microbe_methane) then
       call endrun(msg=' ERROR: microbial-methane restart state is not implemented in Phase 1.'//&
            errMsg(__FILE__, __LINE__))
    end if

  end subroutine Restart

end module MicrobeMethaneMod
