! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

! Trajectory module template file
! If you would like to make custom trajectory I/O module, use this file as a template.
module trajectory
  type handle
     integer :: iohandle
  end type handle
 
contains

  subroutine init_trajectory()
    implicit none
  end subroutine init_trajectory

  subroutine finish_trajectory()
    implicit none
  end subroutine finish_trajectory

  ! Open trajectory and returns handle as htraj. 
  ! Should open fail, the program abends.
  subroutine open_trajectory(htraj, fname)
    use utility, only: newunit
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    ! use form="FORMATTED" for human-redable file format
    open(unit=newunit(htraj%iohandle), file=fname, action="READ", form="UNFORMATTED")

    ! Read header info, skip header, etc.

  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)

    ! release extra resources if necessary

  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real, intent(out) :: crd(3, natom)
    real, intent(out) :: cell(3, 3)
    integer, intent(out) :: status

    read(htraj%iohandle, err = 999, eof = 999) ! read cells, coodinates, etc.
    ! if necessary, angles_to_cell_vector() in module utility will help

    status = 0
    return

999 status = 1
    return
  end subroutine read_trajectory

end module trajectory
