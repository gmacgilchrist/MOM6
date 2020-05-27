!> Regrid columns for the continuous isopycnal (rho) coordinate
module coord_utils

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, DEGREE_MAX

implicit none ; private

public build_scalar_column, sort_scalar_k_1d, copy_finite_thicknesses, old_inflate_layers_1d

!> Build a scalar coordinate column
!!
!! 1. Positions of target scalar values (for interfaces) are found by interpolation.
subroutine build_scalar_column(CS, nz, depth, h, phi, z_interface, ksort, &
                            h_neglect, h_neglect_edge)
  type(rho_CS),        intent(in)    :: CS !< coord_* control structure
  integer,             intent(in)    :: nz !< Number of levels on source grid (i.e. length of  h, phi)
  real,                intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real, dimension(nz), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: phi  !< Scalar for source column
  real, dimension(CS%nk+1), &
                       intent(inout) :: z_interface !< Absolute positions of interfaces
  real, dimension(nz), intent(out)   :: ksort !< k-indices for monotonically increasing column
  real,      optional, intent(in)    :: h_neglect !< A negligibly small width for the purpose
                                             !! of cell reconstructions [H ~> m or kg m-2]
  real,      optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                             !! of edge value calculations [H ~> m or kg m-2]

  ! Local variables
  integer :: k, count_nonzero_layers
  integer, dimension(nz) :: mapping
  real, dimension(nz) :: h_nv     ! Thicknesses of non-vanishing layers [H ~> m or kg m-2]
  real, dimension(nz+1) :: xTmp   ! Temporary positions [H ~> m or kg m-2]
  real, dimension(CS%nk) :: h_new ! New thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk+1) :: x1

  ! Construct source column with vanished layers removed (stored in h_nv)
  call copy_finite_thicknesses(nz, h, CS%min_thickness, count_nonzero_layers, h_nv, mapping)

  if (count_nonzero_layers > 1) then
    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + h_nv(k)
    enddo

    do k = 1,count_nonzero_layers
      phi(k) = phi(mapping(k))
    enddo
    
    ! Determine k-indices for a monotonically increasing column
    ! Do not change the scalar distribution in the column
    call sort_scalar_k_1d(nz, phi, ksort)

    ! Based on source column density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, phi, count_nonzero_layers, &
                                    h_nv, xTmp, CS%target_density, CS%nk, h_new, &
                                    x1, h_neglect, h_neglect_edge)

    ! Inflate vanished layers
    call old_inflate_layers_1d(CS%min_thickness, CS%nk, h_new)

    ! Comment: The following adjustment of h_new, and re-calculation of h_new via x1 needs to be removed
    x1(1) = 0.0 ; do k = 1,CS%nk ; x1(k+1) = x1(k) + h_new(k) ; enddo
    do k = 1,CS%nk
      h_new(k) = x1(k+1) - x1(k)
    enddo

  else ! count_nonzero_layers <= 1
    if (nz == CS%nk) then
      h_new(:) = h(:) ! This keeps old behavior
    else
      h_new(:) = 0.
      h_new(1) = h(1)
    endif
  endif

  ! Return interface positions
  if (CS%integrate_downward_for_e) then
    ! Remapping is defined integrating from zero
    z_interface(1) = 0.
    do k = 1,CS%nk
      z_interface(k+1) = z_interface(k) - h_new(k)
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    z_interface(CS%nk+1) = -depth
    do k = CS%nk,1,-1
      z_interface(k) = z_interface(k+1) + h_new(k)
    enddo
  endif

end subroutine build_scalar_column

!------------------------------------------------------------------------------
!> Return the index of a sorted array of scalar values
subroutine sort_scalar_k_1d(nz, phi, ksort)
  integer,   intent(in)                :: nz    !< number of levels on source grid
  real, dimension(nz), &
                           intent(in)  :: phi  !< Array of scalar quantity to be sorted
  integer, dimension(nz), &
                           intent(out) :: ksort !< An array of indicies for a 
                                                  !! monotonically increasing scalar
!------------------------------------------------------------------------------
! Check each water column to see if a given scalar is monotonically increasing.
! If not, return an array of the sorted indices (bubble sort algorithm).
! No need to return the sorted scalar array itself.
!------------------------------------------------------------------------------

  ! Local variables
  integer   :: k
  real      :: P0, P1       ! temperatures
  logical   :: monotonic

  ! Repeat swapping of indices until complete
  do
    monotonic = .true.
    do k = 1,nz-1
      ! Gather information of scalar value in current and next cells
      P0 = phi(k)  ; P1 = phi(k+1)
      ! If the scalar value of the current cell is larger than the scalar
      ! below it, we swap the cell indices
      if ( P0 > P1 ) then
        ksort(k) = k+1 ; ksort(k+1) = k
        monotonic = .false.
      endif
    enddo  ! k

    if ( monotonic ) exit
  enddo

end subroutine sort_scalar_k_1d

!> Copy column thicknesses with vanished layers removed
subroutine copy_finite_thicknesses(nk, h_in, thresh, nout, h_out, mapping)
  integer,                intent(in)  :: nk      !< Number of layer for h_in, T_in, S_in
  real, dimension(nk),    intent(in)  :: h_in    !< Thickness of input column [H ~> m or kg m-2] or [Z ~> m]
  real,                   intent(in)  :: thresh  !< Thickness threshold defining vanished
                                                 !! layers [H ~> m or kg m-2] or [Z ~> m]
  integer,                intent(out) :: nout    !< Number of non-vanished layers
  real, dimension(nk),    intent(out) :: h_out   !< Thickness of output column [H ~> m or kg m-2] or [Z ~> m]
  integer, dimension(nk), intent(out) :: mapping !< Index of k-out corresponding to k-in
  ! Local variables
  integer :: k, k_thickest
  real :: thickness_in_vanished ! Summed thicknesses in discarded layers [H ~> m or kg m-2] or [Z ~> m]
  real :: thickest_h_out        ! Thickness of the thickest layer [H ~> m or kg m-2] or [Z ~> m]

  ! Build up new grid
  nout = 0
  thickness_in_vanished = 0.0
  thickest_h_out = h_in(1)
  k_thickest = 1
  do k = 1, nk
    mapping(k) = nout ! Note k>=nout always
    h_out(k) = 0.  ! Make sure h_out is set everywhere
    if (h_in(k) > thresh) then
      ! For non-vanished layers
      nout = nout + 1
      mapping(nout) = k
      h_out(nout) = h_in(k)
      if (h_out(nout) > thickest_h_out) then
        thickest_h_out = h_out(nout)
        k_thickest = nout
      endif
    else
      ! Add up mass in vanished layers
      thickness_in_vanished = thickness_in_vanished + h_in(k)
    endif
  enddo

  ! No finite layers
  if (nout <= 1) return

  ! Adjust for any lost volume in vanished layers
  h_out(k_thickest) = h_out(k_thickest) + thickness_in_vanished

end subroutine copy_finite_thicknesses

!------------------------------------------------------------------------------
!> Inflate vanished layers to finite (nonzero) width
subroutine old_inflate_layers_1d( min_thickness, nk, h )

  ! Argument
  real,               intent(in)    :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]
  integer,            intent(in)    :: nk  !< Number of layers in the grid
  real, dimension(:), intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variable
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta
  real      :: correction
  real      :: maxThickness

  ! Count number of nonzero layers
  count_nonzero_layers = 0
  do k = 1,nk
    if ( h(k) > min_thickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    endif
  enddo

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers == nk ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers == 0 ) then
    do k = 1,nk
      h(k) = min_thickness
    enddo
    return
  endif

  ! Inflate zero layers
  correction = 0.0
  do k = 1,nk
    if ( h(k) <= min_thickness ) then
      delta = min_thickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    endif
  enddo

  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,nk
    if ( h(k) > maxThickness ) then
      maxThickness = h(k)
      k_found = k
    endif
  enddo

  h(k_found) = h(k_found) - correction

end subroutine old_inflate_layers_1d

end module coord_utils
