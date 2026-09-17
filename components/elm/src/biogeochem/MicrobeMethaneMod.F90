module MicrobeMethaneMod

  !-----------------------------------------------------------------------
  ! Lifecycle boundary for the microbial decomposition and revised methane
  ! backend. Phase 1 intentionally owns no prognostic state and stops during
  ! initialization. Later phases will add state only behind
  ! use_microbe_methane and will read scientific parameters from ELM's
  ! standard parameter file.
  !-----------------------------------------------------------------------

  use decompMod   , only : bounds_type
  use elm_varctl  , only : use_microbe_methane
  use pio         , only : file_desc_t

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

    ! Phase 2 keeps the established CH4 type active and owns no revised
    ! methane state. These lifecycle calls are connected in Phase 3.

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

    ! Revised methane parameters will be read from ELM's standard parameter
    ! NetCDF in Phase 3. Do not add a module-specific file.

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

    ! No revised methane restart state exists in Phase 2. The established CH4
    ! type owns restart state until the Phase 3 backend is connected.

  end subroutine Restart

end module MicrobeMethaneMod
