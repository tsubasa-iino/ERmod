! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! As a special exception, you may use this file as part of a free software
! without restriction.  Specifically, if other files instantiate
! templates or use macros or inline functions from this file, or you compile
! this file and link it with other files to produce an executable, this
! file does not by itself cause the resulting executable to be covered by
! the GNU General Public License.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!

module utility
   implicit none

   interface
      function hash_double(arr, n) bind(C) result(hashed)
         use, intrinsic :: iso_c_binding
         type(c_ptr), intent(in) :: arr
         integer(kind=c_size_t), value, intent(in) :: n
         integer(kind=c_int64_t) :: hashed
      end function
      function hash_float(arr, n) bind(C) result(hashed)
         use, intrinsic :: iso_c_binding
         type(c_ptr), intent(in) :: arr
         integer(kind=c_size_t), value, intent(in) :: n
         integer(kind=c_int64_t) :: hashed
      end function
   end interface

contains

#ifdef HAVE_TRANSFER
   integer(8) function hash(arr, size)
      implicit none
      integer, intent(in) :: size
      real, intent(in) :: arr(size)
      integer(8) :: ret
      integer :: i
      ret = 0

      do i = 1, size
         ret = ishftc(ret, 7)
         ret = ieor(ret, transfer(arr(i), ret))
      end do

      hash = ret
   end function hash

#else
! (== not HAVE_TRANSFER)
   integer(8) function hash(arr, size)
      use, intrinsic :: iso_c_binding
      implicit none

      integer, intent(in) :: size
      real, intent(in), target :: arr(size)
      integer(kind=c_size_t) :: size_c

      ! Since Fortran do not have template instantiation-like type operation,
      ! and because we may use "real as real(8)" option,
      ! we have no way to selectively call hash_double or hash_float with different types
      ! instead we choose to use generic type of c_ptr (void*) and use C_LOC to achieve this.
      size_c = size
      select case(kind(arr))
       case(4)
         hash = hash_float(c_loc(arr), size_c)
       case(8)
         hash = hash_double(c_loc(arr), size_c)
       case default
         stop "Unexpected real variable size"
      end select

   end function hash
#endif

   ! The following function is a snippet from Fortran wiki and in public domain.
   !
   ! This is a simple function to search for an available unit.
   ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
   ! The UNIT value is returned by the function, and also by the optional
   ! argument. This allows the function to be used directly in an OPEN
   ! statement, and optionally save the result in a local variable.
   ! If no units are available, -1 is returned.
   integer function newunit(unit)
      implicit none
      integer, intent(out), optional :: unit
      ! local
      integer, parameter :: LUN_MIN=100, LUN_MAX=110
      logical :: opened
      integer :: lun
      ! begin
      newunit=-1
      do lun=LUN_MIN,LUN_MAX
         inquire(unit=lun,opened=opened)
         if (.not. opened) then
            newunit=lun
            exit
         end if
      end do
      if (present(unit)) unit=newunit
   end function newunit

   ! convert cell-length & (alpha, beta, gamma) to cell vectors
   subroutine angles_to_cell_vector(cell_len, angles, out_cell_vectors)
      use engmain, only: PI
      implicit none
      real(kind=8), intent(in) :: cell_len(3)
      real(kind=8), intent(in) :: angles(3)
      real(kind=8), intent(out) :: out_cell_vectors(3, 3)
      real(kind=8) :: x, y, u, v, w
      real(kind=8) :: sina, cosa, sinb, cosb, sing, cosg
      real(kind=8) :: blen, clen

      ! all calculations are done in double precision irrespective to the real precisions

      call calc_exact_sincos(real(angles(1), 8), sina, cosa) ! alpha for b-c axes
      call calc_exact_sincos(real(angles(2), 8), sinb, cosb) ! beta for a-c axes
      call calc_exact_sincos(real(angles(3), 8), sing, cosg) ! gamma for a-b axes

      ! ~a = (1, 0, 0)
      ! ~b = (x, y, 0)
      ! ~c = (u, v, w)
      ! ~a.~b = x = cos gamma
      ! |~a*~b| = y = sin gamma
      ! ~a.~c = u = cos beta
      ! ~b.~c = xu + yv = cos alpha

      x = cosg
      y = sing
      u = cosb
      v = (cosa - x * u) / y
      w = sqrt(1.0d0 - u * u - v * v)

      out_cell_vectors(1, 1) = cell_len(1)
      out_cell_vectors(2, 1) = 0.0
      out_cell_vectors(3, 1) = 0.0

      blen = real(cell_len(2), 8)
      out_cell_vectors(1, 2) = real(blen * x)
      out_cell_vectors(2, 2) = real(blen * y)
      out_cell_vectors(3, 2) = 0.0

      clen = real(cell_len(3), 8)
      out_cell_vectors(1, 3) = real(clen * u)
      out_cell_vectors(2, 3) = real(clen * v)
      out_cell_vectors(3, 3) = real(clen * w)

   contains
      subroutine calc_exact_sincos(angle_in_degree, sinres, cosres)
         real(8), intent(in) :: angle_in_degree
         real(8), intent(out) :: sinres, cosres
         real(8) :: angle_in_radian

         if (angle_in_degree == 90.0) then
            sinres = 1.0
            cosres = 0.0
         else
            angle_in_radian = angle_in_degree * PI / 180.0d0
            sinres = sin(angle_in_degree)
            cosres = cos(angle_in_degree)
         end if
      end subroutine
   end subroutine angles_to_cell_vector

   character(len=16) function itoa(x)
      integer, intent(in) :: x
      character(len=16) :: buf
      write(buf,"(I16)") x
      itoa = buf
   end function itoa
end module utility

