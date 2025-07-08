module rrtmgp_inputs_cam

!--------------------------------------------------------------------------------
! Transform data for inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.  Add an extra layer if CAM's top is below 1 Pa.
! The vertical indexing increases from top to bottom of atmosphere in both
! CAM and RRTMGP arrays.   
!--------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use cam_logfile, only: iulog

use physconst,        only: stebol, pi

use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc
use camsrfexch,       only: cam_in_t

use radconstants,     only: nradgas, gaslist, nswbands, nlwbands

use rad_constituents, only: rad_cnst_get_gas

use cloud_rad_props,  only: get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
                            get_ice_optics_sw,    ice_cloud_get_rad_props_lw,    &
                            get_snow_optics_sw,   snow_cloud_get_rad_props_lw,   &
                            get_grau_optics_sw,   grau_cloud_get_rad_props_lw
                                 
use mcica_subcol_gen, only: mcica_subcol_sw

use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw

use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
use ccpp_optical_props,      only: ty_optical_props_2str_ccpp, ty_optical_props_1scl_ccpp

use cam_history_support,   only: fillvalue
use cam_abortutils,   only: endrun
use error_messages,   only: alloc_err
use radiation_utils, only: get_sw_spectral_boundaries_ccpp
use spmd_utils, only: masterproc
use spmd_utils, only: iam

implicit none
private
save

public :: &
   rrtmgp_inputs_cam_init, &
   rrtmgp_get_gas_mmrs, &
   rrtmgp_set_gases_sw, &
   rrtmgp_set_aer_lw,   &
   rrtmgp_set_aer_sw


! This value is to match the arbitrary small value used in RRTMG to decide
! when a quantity is effectively zero.
real(r8), parameter :: tiny = 1.0e-80_r8
real(r8) :: sw_low_bounds(nswbands)
real(r8) :: sw_high_bounds(nswbands)
integer :: ktopcam
integer :: ktoprad
integer :: idx_sw_diag
integer :: idx_nir_diag
integer :: idx_uv_diag
integer :: idx_sw_cloudsim
integer :: idx_lw_diag
integer :: idx_lw_cloudsim

! Mapping from RRTMG shortwave bands to RRTMGP.  Currently needed to continue using
! the SW optics datasets from RRTMG (even thought there is a slight mismatch in the
! band boundaries of the 2 bands that overlap with the LW bands).
integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
   [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

!==================================================================================================
contains
!==================================================================================================

!==================================================================================================
subroutine rrtmgp_inputs_cam_init(ktcam, ktrad, idx_sw_diag_in, idx_nir_diag_in, idx_uv_diag_in, &
               idx_sw_cloudsim_in, idx_lw_diag_in, idx_lw_cloudsim_in)
      
   ! Note that this routine must be called after the calls to set_wavenumber_bands which set
   ! the sw/lw band boundaries in the radconstants module.

   integer, intent(in) :: ktcam
   integer, intent(in) :: ktrad
   integer, intent(in) :: idx_sw_diag_in
   integer, intent(in) :: idx_nir_diag_in
   integer, intent(in) :: idx_uv_diag_in
   integer, intent(in) :: idx_sw_cloudsim_in
   integer, intent(in) :: idx_lw_diag_in
   integer, intent(in) :: idx_lw_cloudsim_in
   character(len=512) :: errmsg
   integer :: errflg

   ktopcam = ktcam
   ktoprad = ktrad
   idx_sw_diag = idx_sw_diag_in
   idx_nir_diag = idx_nir_diag_in
   idx_uv_diag = idx_uv_diag_in
   idx_sw_cloudsim = idx_sw_cloudsim_in
   idx_lw_diag = idx_lw_diag_in
   idx_lw_cloudsim = idx_lw_cloudsim_in

   ! Initialize the module data containing the SW band boundaries.
   call get_sw_spectral_boundaries_ccpp(sw_low_bounds, sw_high_bounds, 'cm^-1', errmsg, errflg)
   if (errflg /= 0) then
      call endrun('rrtmgp_inputs_cam_init: error during get_sw_spectral_boundaries_ccpp - message: '//errmsg)
   end if

end subroutine rrtmgp_inputs_cam_init

!=========================================================================================

function get_molar_mass_ratio(gas_name) result(massratio)

   ! return the molar mass ratio of dry air to gas based on gas_name

   character(len=*),intent(in) :: gas_name
   real(r8)                    :: massratio

   ! local variables
   real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
   real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
   real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
   real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
   real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
   real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
   real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
   real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

   character(len=*), parameter :: sub='get_molar_mass_ratio'
   !----------------------------------------------------------------------------

   select case (trim(gas_name)) 
      case ('H2O') 
         massratio = amdw
      case ('CO2')
         massratio = amdc
      case ('O3')
         massratio = amdo
      case ('CH4')
         massratio = amdm
      case ('N2O')
         massratio = amdn
      case ('O2')
         massratio = amdo2
      case ('CFC11')
         massratio = amdc1
      case ('CFC12')
         massratio = amdc2
      case default
         call endrun(sub//": Invalid gas: "//trim(gas_name))
   end select

end function get_molar_mass_ratio

!=========================================================================================

subroutine rad_gas_get_vmr(icall, gas_name, state, pbuf, nlay, numactivecols, gas_concs, idxday)

   ! Set volume mixing ratio in gas_concs object.
   ! The gas_concs%set_vmr method copies data into internally allocated storage.

   integer,                     intent(in) :: icall      ! index of climate/diagnostic radiation call
   character(len=*),            intent(in) :: gas_name
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)
   integer,                     intent(in) :: nlay           ! number of layers in radiation calculation
   integer,                     intent(in) :: numactivecols  ! number of columns, ncol for LW, nday for SW

   type(ty_gas_concs_ccpp),     intent(inout) :: gas_concs  ! the result is VMR inside gas_concs

   integer, optional,          intent(in) :: idxday(:)   ! indices of daylight columns in a chunk

   ! Local variables
   integer :: i, idx(numactivecols)
   integer :: istat
   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)
   real(r8), allocatable :: mmr(:,:)
   real(r8) :: massratio

   ! For ozone profile above model
   real(r8) :: P_top, P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rad_gas_get_vmr'
   !----------------------------------------------------------------------------

   ! set the column indices; when idxday is provided (e.g. daylit columns) use them, otherwise just count.
   do i = 1, numactivecols
      if (present(idxday)) then
         idx(i) = idxday(i)
      else
         idx(i) = i
      end if
   end do

   ! gas_mmr points to a "chunk" in either the state or pbuf objects.  That storage is
   ! dimensioned (pcols,pver).
   call rad_cnst_get_gas(icall, gas_name, state, pbuf, gas_mmr)

   ! Copy into storage for RRTMGP
   allocate(mmr(numactivecols, nlay), stat=istat)
   call alloc_err(istat, sub, 'mmr', numactivecols*nlay)
   allocate(gas_vmr(numactivecols, nlay), stat=istat)
   call alloc_err(istat, sub, 'gas_vmr', numactivecols*nlay)

   do i = 1, numactivecols
      mmr(i,ktoprad:) = gas_mmr(idx(i),ktopcam:)
   end do

   ! If an extra layer is being used, copy mmr from the top layer of CAM to the extra layer.
   if (nlay == pverp) then
      mmr(:,1) = mmr(:,2)
   end if

   ! special case: H2O is specific humidity, not mixing ratio. Use r = q/(1-q):
   if (gas_name == 'H2O') then 
      mmr = mmr / (1._r8 - mmr)
   end if  

   ! convert MMR to VMR, multipy by ratio of dry air molar mas to gas molar mass.
   massratio = get_molar_mass_ratio(gas_name)
   gas_vmr = mmr * massratio

   ! special case: Setting O3 in the extra layer:
   ! 
   ! For the purpose of attenuating solar fluxes above the CAM model top, we assume that ozone 
   ! mixing decreases linearly in each column from the value in the top layer of CAM to zero at 
   ! the pressure level set by P_top. P_top has been set to 50 Pa (0.5 hPa) based on model tuning 
   ! to produce temperatures at the top of CAM that are most consistent with WACCM at similar pressure levels. 

   if ((gas_name == 'O3') .and. (nlay == pverp)) then
      P_top = 50.0_r8
      do i = 1, numactivecols
            P_int = state%pint(idx(i),1) ! pressure (Pa) at upper interface of CAM
            P_mid = state%pmid(idx(i),1) ! pressure (Pa) at midpoint of top layer of CAM
            alpha = log(P_int/P_top)
            beta =  log(P_mid/P_int)/log(P_mid/P_top)
      
            a =  ( (1._r8 + alpha) * exp(-alpha) - 1._r8 ) / alpha
            b =  1._r8 - exp(-alpha)
   
            if (alpha .gt. 0) then             ! only apply where top level is below 80 km
               chi_mid = gas_vmr(i,1)          ! molar mixing ratio of O3 at midpoint of top layer
               chi_0 = chi_mid /  (1._r8 + beta)
               chi_eff = chi_0 * (a + b)
               gas_vmr(i,1) = chi_eff
            end if
      end do
   end if

   errmsg = gas_concs%gas_concs%set_vmr(gas_name, gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR, gas_concs%set_vmr: '//trim(errmsg))
   end if

   deallocate(gas_vmr)
   deallocate(mmr)

end subroutine rad_gas_get_vmr

!==================================================================================================

subroutine rrtmgp_get_gas_mmrs(icall, state, pbuf, nlay, gas_mmrs)

   ! Retrieve mass mixing ratios for radiatively active gases from rad_constituents

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   real(r8),                    intent(out)   :: gas_mmrs(:,:,:)

   ! local variables
   integer :: i, ncol
   real(r8), pointer :: gas_mmr(:,:)
   character(len=*), parameter :: sub = 'rrtmgp_get_gas_mmrs'
   !--------------------------------------------------------------------------------

   ncol = state%ncol
   do i = 1, nradgas
      call rad_cnst_get_gas(icall, gaslist(i), state, pbuf, gas_mmr)
      gas_mmrs(:,:,i) = gas_mmr(:ncol,:)
   end do
end subroutine rrtmgp_get_gas_mmrs

!==================================================================================================

subroutine rrtmgp_set_gases_sw( &
   icall, state, pbuf, nlay, nday, &
   idxday, gas_concs)

   ! Return gas_concs with gas volume mixing ratio on DAYLIT columns.
   ! Set all gases in radconstants gaslist.

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   integer,                     intent(in)    :: nday
   integer,                     intent(in)    :: idxday(:)
   type(ty_gas_concs_ccpp),     intent(inout) :: gas_concs

   ! local variables
   integer :: i
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_sw'
   !----------------------------------------------------------------------------

   ! use the optional argument idxday to specify which columns are sunlit
    do i = 1,nradgas
      call rad_gas_get_vmr(icall, gaslist(i), state, pbuf, nlay, nday, gas_concs, idxday=idxday)
   end do

end subroutine rrtmgp_set_gases_sw

!==================================================================================================

subroutine rrtmgp_set_aer_lw(icall, state, pbuf, aer_lw)

   ! Load LW aerosol optical properties into the RRTMGP object.

   ! Arguments
   integer,                     intent(in) :: icall
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   type(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw

   ! Local variables
   integer :: ncol

   ! Aerosol LW absorption optical depth
   real(r8) :: aer_lw_abs (pcols,pver,nlwbands)

   character(len=*), parameter :: sub = 'rrtmgp_set_aer_lw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! Get aerosol longwave optical properties.
   call aer_rad_props_lw(icall, state, pbuf, aer_lw_abs)

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   aer_lw%optical_props%tau = 0.0_r8

   aer_lw%optical_props%tau(:ncol, ktoprad:, :) = aer_lw_abs(:ncol, ktopcam:, :)

   errmsg = aer_lw%optical_props%validate()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: aer_lw%optical_props%validate: '//trim(errmsg))
   end if
end subroutine rrtmgp_set_aer_lw

!==================================================================================================

subroutine rrtmgp_set_aer_sw( &
   icall, state, pbuf, nday, idxday, nnite, idxnite, aer_sw)

   ! Load SW aerosol optical properties into the RRTMGP object.

   ! Arguments
   integer,                     intent(in) :: icall
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   integer,  intent(in) :: nday
   integer,  intent(in) :: idxday(:)
   integer,  intent(in) :: nnite          ! number of night columns
   integer,  intent(in) :: idxnite(pcols) ! indices of night columns in the chunk

   type(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw

   ! local variables
   integer  :: i

   ! The optical arrays dimensioned in the vertical as 0:pver.
   ! The index 0 is for the extra layer used in the radiation
   ! calculation.  The index ktopcam assumes the CAM vertical indices are
   ! in the range 1:pver, so using this index correctly ignores vertical
   ! index 0.  If an "extra" layer is used in the calculations, it is
   ! provided and set in the RRTMGP aerosol object aer_sw.
   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! asymmetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! forward scattered fraction * w * tau
                                                  ! aer_tau_w_f is not used by RRTMGP.
   character(len=*), parameter :: sub = 'rrtmgp_set_aer_sw'
   !--------------------------------------------------------------------------------

   ! Get aerosol shortwave optical properties.
   ! Make outfld calls for aerosol optical property diagnostics.
   call aer_rad_props_sw( &
      icall, state, pbuf, nnite, idxnite, &
      aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)

   ! The aer_sw object is only initialized if nday > 0.
   if (nday > 0) then

      ! aerosol optical properties need to be re-ordered from the RRTMG spectral bands
      ! (as assumed in the optics datasets) to the RRTMGP band order.
      aer_tau(:,:,:)     = aer_tau(    :,:,rrtmg_to_rrtmgp_swbands)
      aer_tau_w(:,:,:)   = aer_tau_w(  :,:,rrtmg_to_rrtmgp_swbands)
      aer_tau_w_g(:,:,:) = aer_tau_w_g(:,:,rrtmg_to_rrtmgp_swbands)
                  
      ! If there is an extra layer in the radiation then this initialization
      ! will provide default values.
      aer_sw%optical_props%tau = 0.0_r8
      aer_sw%optical_props%ssa = 1.0_r8
      aer_sw%optical_props%g   = 0.0_r8

      ! CAM fields are products tau, tau*ssa, tau*ssa*asy
      ! Fields expected by RRTMGP are computed by
      ! aer_sw%optical_props%tau = aer_tau
      ! aer_sw%optical_props%ssa = aer_tau_w / aer_tau
      ! aer_sw%optical_props%g   = aer_tau_w_g / aer_taw_w
      ! aer_sw arrays have dimensions of (nday,nlay,nswbands)

      do i = 1, nday
         ! set aerosol optical depth, clip to zero
         aer_sw%optical_props%tau(i,ktoprad:,:) = max(aer_tau(idxday(i),ktopcam:,:), 0._r8)
         ! set value of single scattering albedo
         aer_sw%optical_props%ssa(i,ktoprad:,:) = merge(aer_tau_w(idxday(i),ktopcam:,:)/aer_tau(idxday(i),ktopcam:,:), &
                                          1._r8, aer_tau(idxday(i),ktopcam:,:) > 0._r8)
         ! set value of asymmetry
         aer_sw%optical_props%g(i,ktoprad:,:) = merge(aer_tau_w_g(idxday(i),ktopcam:,:)/aer_tau_w(idxday(i),ktopcam:,:), &
                                        0._r8, aer_tau_w(idxday(i),ktopcam:,:) > tiny)
      end do

      ! impose limits on the components
      aer_sw%optical_props%ssa = min(max(aer_sw%optical_props%ssa, 0._r8), 1._r8)
      aer_sw%optical_props%g   = min(max(aer_sw%optical_props%g, -1._r8), 1._r8)

   end if

end subroutine rrtmgp_set_aer_sw

!==================================================================================================

end module rrtmgp_inputs_cam
