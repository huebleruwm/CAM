module radiation

! stub module

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl

implicit none
private
save

public :: &
   radiation_readnl,         &
   radiation_do
! Top of valid pressure range (Pa) for this radiation scheme
! in local thermo. equilibrium. No limits for this
! scheme, so
! will be set to zero below
public :: p_top_for_equil_rad

real(r8), public, protected :: nextsw_cday = -1._r8 ! future radiation calday for surface models
real(r8) :: p_top_for_equil_rad = 0._r8

!========================================================================================
contains
!========================================================================================

subroutine radiation_readnl(nlfile)

   ! this stub can be called, but does nothing

   character(len=*), intent(in) :: nlfile

end subroutine radiation_readnl

!========================================================================================

function radiation_do(op, timestep)

   ! Returns true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value
   !---------------------------------------------------------------------------

   radiation_do = .false.

end function radiation_do

!========================================================================================

end module radiation

