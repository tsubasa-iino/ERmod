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

! module that governs trajectory I/O

module trajectory
   use, intrinsic :: iso_c_binding
   type handle
      type(c_ptr) :: vmdhandle
   end type handle

   interface
      subroutine vmdfio_init_traj() bind(C)
      end subroutine

      subroutine vmdfio_fini_traj() bind(C)
      end subroutine

      function vmdfio_open_traj(handle, fname) bind(C) result(retstatus)
         use, intrinsic :: iso_c_binding
         type(c_ptr), intent(out) :: handle
         character(kind=c_char), dimension(*), intent(in) :: fname
         integer(kind=c_int) :: retstatus
      end function

      subroutine vmdfio_close_traj(handle) bind(C)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value, intent(in) :: handle
      end subroutine

      function vmdfio_read_traj_step(handle, xout, box, natoms_aux) bind(C) result(retstatus)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value, intent(in) :: handle
         real(kind=c_float), intent(out) :: xout(*)
         real(kind=c_float), intent(out) :: box(3, 3)
         integer(kind=c_int), value, intent(in) :: natoms_aux
         integer(kind=c_int) :: retstatus
      end function
   end interface

contains
   subroutine init_trajectory()
      implicit none
      call vmdfio_init_traj()
   end subroutine init_trajectory

   subroutine finish_trajectory()
      implicit none
      call vmdfio_fini_traj()
   end subroutine finish_trajectory

   ! Open trajectory and returns handle as htraj.
   ! Should open fail, the program abends.
   subroutine open_trajectory(htraj, fname)
      use, intrinsic :: iso_c_binding
      implicit none
      type(handle), intent(inout) :: htraj
      character(len=*), intent(in) :: fname

      integer :: status

      status = vmdfio_open_traj(htraj%vmdhandle, trim(fname) // c_null_char)
      if(status /= 0) then
         stop "vmdfio_open_traj: unable to open trajectory. HISTORY must be a symlink"
      endif
   end subroutine open_trajectory

   ! Close trajectory specified by handle
   subroutine close_trajectory(htraj)
      implicit none
      type(handle), intent(inout) :: htraj

      call vmdfio_close_traj(htraj%vmdhandle)
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
#ifdef DP
      real(kind=4) :: crd_tmp(3, natom)
      real(kind=4) :: cell_tmp(3, 3)
#endif
      integer, intent(out) :: status

#ifdef DP
      status = vmdfio_read_traj_step(htraj%vmdhandle, crd_tmp, cell_tmp, natom)
      crd = real(crd_tmp, kind=8)
      cell = real(cell_tmp, kind=8)
#else
      status = vmdfio_read_traj_step(htraj%vmdhandle, crd, cell, natom)
#endif
   end subroutine read_trajectory

end module trajectory
