
module radheat
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
!
! Original version: B.A. Boville
! Change weighting function for RRTMG: A J Conley
!-----------------------------------------------------------------------

! Use a cubic polynomial over the domain from minimum pressure to maximum pressure
! Cubic polynomial is chosen so that derivative is zero at minimum and maximum pressures
!   and is monotonically increasing from zero at minimum pressure to one at maximum pressure

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver
  use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
  use physconst,       only: gravit,cpair,mwco2
  use air_composition, only: cpairv
  use perf_mod
  use cam_logfile,     only: iulog

  implicit none
  private
  save

! Public interfaces
  public  &
       radheat_readnl,        &!
       radheat_register,      &!
       radheat_init,          &!
       radheat_timestep_init, &!
       radheat_tend            ! return net radiative heating

  public :: radheat_disable_waccm ! disable waccm heating in the upper atm


! Private variables for merging heating rates
  real(r8):: qrs_wt(pver)             ! merge weight for cam solar heating
  real(r8):: qrl_wt(pver)             ! merge weight for cam long wave heating

  logical :: waccm_heating
  logical :: waccm_heating_on = .true.

  ! sw merge region
  ! highest altitude (lowest  pressure) of merge region (Pa)
  real(r8) :: min_pressure_sw= 5._r8
  ! lowest  altitude (lowest  pressure) of merge region (Pa)
  real(r8) :: max_pressure_sw=50._r8

  ! lw merge region
  ! highest altitude (lowest  pressure) of merge region (Pa)
  real(r8) :: min_pressure_lw= 5._r8
  ! lowest  altitude (highest pressure) of merge region (Pa)
  real(r8) :: max_pressure_lw=50._r8

  integer :: ntop_qrs_cam             ! top level for pure cam solar heating

!===============================================================================
contains
!===============================================================================
subroutine radheat_readnl(nlfile)

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! No options for this version of radheat; this is just a stub.

end subroutine radheat_readnl

!================================================================================================

  subroutine radheat_register

  ! No options for this version of radheat; this is just a stub.

  end subroutine radheat_register

!================================================================================================
!================================================================================================
!================================================================================================

  subroutine radheat_init(pref_mid)

    
    use nlte_fomichev,    only: nlte_fomichev_init
    use cam_history,      only: add_default, addfld
    use phys_control,     only: phys_getopts
    use physics_buffer,   only : physics_buffer_desc

    ! args

    real(r8),          intent(in) :: pref_mid(pver) ! mid point reference pressure
    ! local vars
    real(r8) :: co2_mw, o1_mw, o2_mw, o3_mw, no_mw, n2_mw ! molecular weights

    real(r8) :: delta_merge_sw      ! range of merge region
    real(r8) :: midpoint_sw         ! midpoint of merge region
    real(r8) :: delta_merge_lw      ! range of merge region
    real(r8) :: midpoint_lw         ! midpoint of merge region
    real(r8) :: psh(pver)           ! pressure scale height
    integer  :: k
    logical :: camrt

    character(len=16) :: rad_pkg
    logical :: history_scwaccm_forcing
    logical :: history_waccm
    logical :: nlte_limit_co2

!-----------------------------------------------------------------------


    call phys_getopts(radiation_scheme_out=rad_pkg, &
                      history_waccm_out=history_waccm, &
                      history_scwaccm_forcing_out=history_scwaccm_forcing)
    camrt = rad_pkg == 'CAMRT' .or. rad_pkg == 'camrt'

    ! set max/min pressures for merging regions.

    if (camrt) then
       min_pressure_sw = 1e5_r8*exp(-10._r8)
       max_pressure_sw = 1e5_r8*exp(-9._r8)
       min_pressure_lw = 1e5_r8*exp(-10._r8)
       max_pressure_lw = 1e5_r8*exp(-8.57_r8)
    else
       min_pressure_sw =  5._r8
       max_pressure_sw = 50._r8
       min_pressure_lw = 10._r8
       max_pressure_lw = 50._r8
    endif

    delta_merge_sw = max_pressure_sw - min_pressure_sw
    delta_merge_lw = max_pressure_lw - min_pressure_lw

    midpoint_sw = (max_pressure_sw + min_pressure_sw)/2._r8
    midpoint_lw = (max_pressure_lw + min_pressure_lw)/2._r8

    do k=1,pver

       ! pressure scale heights for camrt merging (waccm4)
       psh(k)=log(1e5_r8/pref_mid(k))

       if ( pref_mid(k) .le. min_pressure_sw  ) then
          qrs_wt(k) = 0._r8
       else if( pref_mid(k) .ge. max_pressure_sw) then
          qrs_wt(k) = 1._r8
       else
          if (camrt) then
             ! camrt
             qrs_wt(k) = 1._r8 - tanh( (psh(k) - 9._r8)/.75_r8 )
          else
             ! rrtmg
             qrs_wt(k) = 0.5_r8 + 1.5_r8*((pref_mid(k)-midpoint_sw)/delta_merge_sw) &
                       - 2._r8*((pref_mid(k)-midpoint_sw)/delta_merge_sw)**3._r8
          endif
       endif

       if ( pref_mid(k) .le. min_pressure_lw  ) then
          qrl_wt(k)= 0._r8
       else if( pref_mid(k) .ge. max_pressure_lw) then
          qrl_wt(k)= 1._r8
       else
          if (camrt) then
             ! camrt
             qrl_wt(k) = 1._r8 - tanh( (psh(k) - 8.57_r8) / 0.71_r8 )
          else
             ! rrtmg
             qrl_wt(k) = 0.5_r8 + 1.5_r8*((pref_mid(k)-midpoint_lw)/delta_merge_lw) &
                       - 2._r8*((pref_mid(k)-midpoint_lw)/delta_merge_lw)**3._r8
          endif
       endif

    end do
    
    co2_mw = mwco2
    o1_mw  = 16._r8
    o2_mw  = 32._r8
    o3_mw  = 48._r8
    no_mw  = 30._r8
    n2_mw  = 28._r8
    nlte_limit_co2 = .true.
! Initialize Fomichev parameterization
    call nlte_fomichev_init (co2_mw, n2_mw, o1_mw, o2_mw, o3_mw, no_mw, nlte_limit_co2)
    

    ! determine upppermost level that is purely solar heating (no MLT chem heationg)
    ntop_qrs_cam = 0
    do k=pver,1,-1
       if (qrs_wt(k)==1._r8) ntop_qrs_cam = k
    enddo


    ! WACCM heating if top-most layer is above merge region
    waccm_heating = (pref_mid(1) .le. min_pressure_sw)

    if (masterproc) then
       write(iulog,*) 'WACCM Heating is computed (true/false): ',waccm_heating
    end if

  end subroutine radheat_init

!================================================================================================

  subroutine radheat_timestep_init (state, pbuf2d)

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in):: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)



  end subroutine radheat_timestep_init

!================================================================================================

  subroutine radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
       fsnt, flns, flnt, asdir, net_flx)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!
! This routine provides the waccm hook for computing nonLTE cooling and
! eUV heating.
!-----------------------------------------------------------------------

    use cam_history,           only: outfld
    use nlte_fomichev,         only: nlte_fomichev_calc

    use physics_buffer,        only : physics_buffer_desc
    use constituents,          only: cnst_get_ind
    use calculate_net_heating, only: calculate_net_heating_run
    use cam_abortutils,        only: endrun

!+++arh
    use rad_constituents,      only: rad_cnst_get_gas

! Arguments
    type(physics_state), intent(in)  :: state             ! Physics state variables

    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
    real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating
    real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating
    real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
    real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
    real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
    real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
    real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
    real(r8),            intent(out) :: net_flx(pcols)

! Local variables
    integer  :: i, k , ico2
    integer  :: ncol                                ! number of atmospheric columns
    integer  :: lchnk                               ! chunk identifier
    real(r8) :: qrl_mrg(pcols,pver)                 ! merged LW heating
    real(r8) :: qrl_mlt(pcols,pver)                 ! M/LT longwave heating rates
    real(r8) :: qrs_mrg(pcols,pver)                 ! merged SW heating
    real(r8) :: qrs_mlt(pcols,pver)                 ! M/LT solar heating rates
    real(r8) :: qout(pcols,pver)                    ! temp for outfld call
    real(r8) :: dcoef(6)                            ! for tidal component of heating
    real(r8) :: tau_newt                            ! time scale for 'IR' relaxation

    real(r8) :: qrlfomichev(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: o3cool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: co2cool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: c2scool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: xco2mmr(pcols,pver) ! CO2
    real(r8) :: xo2mmr(pcols,pver)  ! O2
    real(r8) :: xo3mmr(pcols,pver)  ! O3
    real(r8) :: xommr(pcols,pver)   ! O
    real(r8) :: xhmmr(pcols,pver)   ! H
    real(r8) :: xnommr(pcols,pver)  ! NO
    real(r8) :: xn2mmr(pcols,pver)  ! N2

!+++arh
    integer  :: icall
    real(r8), pointer  :: gas_mmr(:,:)
    character(len=512) :: errmsg
    integer            :: errflg

!-----------------------------------------------------------------------

    ncol  = state%ncol
    lchnk = state%lchnk
    
    call cnst_get_ind( 'CO2', ico2 )
    xco2mmr = state%q(:,:,ico2)

    call physics_ptend_init(ptend, state%psetcols, 'radheat', ls=.true.)

    qrs_mlt(:,:) = 0._r8
 
! Merge cam solar heating for lower atmosphere with M/LT heating
!    call merge_qrs (ncol, qrs, qrs_mlt, qrs_mrg, cpairv(:,:,lchnk))
!    qout(:ncol,:) = qrs_mrg(:ncol,:)/cpair
!++jtb
! just use RRTMG's
    qrs_mrg(:,:) = qrs(:,:)
    
!++jtb

#if 0
!  think I need to scale w/ cpair   
    tau_newt = 1._r8 * 86400._r8 ! 1 day Newtonian timescale
    do k = 1, pver
       do i = 1, ncol
          qrl_mlt(i,k) = -cpair * ( state%t(i,k) - 200._r8 ) / tau_newt
       end do
    end do
#else

!+++arh
   icall = 0

! can't find gas name 'O'
   !call rad_cnst_get_gas(icall,'O    ', state, pbuf, gas_mmr)
   !xommr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   !nullify(gas_mmr)
   xommr(:pcols,:pver) = 0._r8

   call rad_cnst_get_gas(icall,'O2   ', state, pbuf, gas_mmr)
   xo2mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

   call rad_cnst_get_gas(icall,'O3   ', state, pbuf, gas_mmr)
   xo3mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

! can't find gas name 'N2'
   !call rad_cnst_get_gas(icall,'N2   ', state, pbuf, gas_mmr)
   !xn2mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   !nullify(gas_mmr)
   !++ jtb set to 0.7547 ('standard' mass mixing ratio of nitrogen gas in atmos)
   xn2mmr(:pcols,:pver) =  0.7547_r8 !0._r8

   call rad_cnst_get_gas(icall,'CO2  ', state, pbuf, gas_mmr)
   xco2mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

    !xommr   = 0._r8
    !xo3mmr  = 0._r8
    !xo2mmr  = 0.2320_r8
    !xn2mmr  = 0.7547_r8
    !xhmmr   = 0._r8
    !xnommr  = 0._r8
    call nlte_fomichev_calc (lchnk,ncol,state%pmid,state%pint,state%t, &
         xo2mmr,xommr,xo3mmr,xn2mmr,xco2mmr,qrlfomichev,co2cool,o3cool,c2scool)
    qrl_mlt = qrlfomichev
#endif

!   Merge cam long wave heating for lower atmosphere with M/LT (nlte) heating
    call merge_qrl (ncol, qrl, qrl_mlt, qrl_mrg)
    qout(:ncol,:) = qrl_mrg(:ncol,:)/cpair

    ! REMOVECAM no longer need once CAM is retired and pcols doesn't exist
    net_flx = 0._r8
    ptend%s = 0._r8
    ! END_REMOVECAM

    call calculate_net_heating_run(ncol, pver, ptend%s(:ncol,:), qrl_mrg(:ncol,:), qrs_mrg(:ncol,:), &
            gravit, state%pdel(:ncol,:), net_flx(:ncol), errmsg, errflg)

  end subroutine radheat_tend

!================================================================================================
  subroutine radheat_disable_waccm()
    waccm_heating_on = .false.
  end subroutine radheat_disable_waccm
!================================================================================================

  subroutine merge_qrs (ncol, hcam, hmlt, hmrg, cpair)
!
!  Merges short wave heating rates
!
    implicit none

!-----------------Input arguments-----------------------------------
    integer ncol

    real(r8), intent(in)  :: hmlt(pcols,pver)                ! Upper atmosphere heating rates
    real(r8), intent(in)  :: hcam(pcols,pver)                ! CAM heating rate
    real(r8), intent(out) :: hmrg(pcols,pver)                ! merged heating rates
    real(r8), intent(in)  :: cpair(pcols,pver)               ! Specific heat of dry air

!-----------------Local workspace------------------------------------

    integer k

    do k = 1, pver
       hmrg(:ncol,k) = qrs_wt(k)*hcam(:ncol,k) + (1._r8 - qrs_wt(k))*cpair(:ncol,k)*hmlt(:ncol,k)
    end do

  end subroutine merge_qrs

!==================================================================================================

  subroutine merge_qrl (ncol, hcam, hmlt, hmrg)
!
!  Merges long wave heating rates
!
!-----------------Input arguments-----------------------------------
    integer ncol

    real(r8), intent(in)  :: hmlt(pcols,pver)               ! Upper atmosphere heating rates
    real(r8), intent(in)  :: hcam(pcols,pver)               ! CAM heating rate
    real(r8), intent(out) :: hmrg(pcols,pver)               ! merged heating rates

!-----------------Local workspace------------------------------------

    integer k

!--------------------------------------------------------------------

    do k = 1, pver
       hmrg(:ncol,k) = qrl_wt(k) * hcam(:ncol,k) + (1._r8-qrl_wt(k)) * hmlt(:ncol,k)
    end do

  end subroutine merge_qrl

end module radheat
