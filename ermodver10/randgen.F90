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


! I wished to move xoshiro into the Fortran code, but it was impossible because
! Fortran standard do not have unsigned integer
! Example usage:
! use randgen
! type(randstate) :: s
!
! call s%init(seed)
! print *, s%next_real()
! print *, s%next_int31()
! FIXME TODO: you can't copy and assign randstate variables.
! To fix, we will need to implement the assignment operator.
! In other words:
! copy = s  !doesn't work!
module randgen
   use, intrinsic :: iso_c_binding
   implicit none
   type, public :: randstate
      type(c_ptr), private :: st = c_null_ptr
   contains
      ! I know in principle I should make dtor to be a "final" keyword, 
      ! But finalization timing is currently compiler dependent (https://stackoverflow.com/questions/59985499/are-fortrans-final-subroutines-reliable-enough-for-practical-use)
      ! This result in the function return value (i.e., module procedure function used as ctor) may or may not be finalized and may end up double-free.
      ! Thus I decided to leave the dtor handled manually for a while.
      procedure, public :: init, dtor
      procedure, public :: next_double, next_real, next_int31
      procedure, public :: jump, long_jump
   end type

   private
contains
   subroutine init(this, seed)
      implicit none
      class(randstate) :: this
      integer(8), intent(in) :: seed
      integer(kind=c_int64_t) :: seed_in
      interface
         function xoshiro256ss_init_state_with_seed(seed) bind(C)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int64_t), intent(in), value :: seed
            type(c_ptr) :: xoshiro256ss_init_state_with_seed
         end function
      end interface

      seed_in = int(seed, c_int64_t)
      this%st = xoshiro256ss_init_state_with_seed(seed_in)
      end subroutine

   subroutine dtor(this)
      implicit none
      class(randstate) :: this
      interface
         subroutine xoshiro256ss_fini(st) bind(C)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: st
         end subroutine
      end interface

      call xoshiro256ss_fini(this%st)
      this%st = c_null_ptr
   end subroutine

   ! returns [0..1)
   function next_double(this)
      implicit none
      class(randstate) :: this
      real(8) :: next_double
      real(kind=c_double) :: xvalue
      interface
         function xoshiro256ss_next_double(st) bind(C)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: st
            real(kind=c_double) :: xoshiro256ss_next_double
         end function
      end interface

      xvalue = xoshiro256ss_next_double(this%st)
      next_double = xvalue
   end function next_double

   function next_real(this)
      implicit none
      class(randstate) :: this
      real :: next_real
      real(8) :: temp

      temp = next_double(this)
      next_real = real(temp)
   end function next_real

   ! returns 0 <= x < 2^31
   function next_int31(this)
      implicit none
      class(randstate) :: this
      integer :: next_int31
      integer(kind=c_int32_t) :: xvalue
      interface
         function xoshiro256ss_next_int31(st) bind(C)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: st
            integer(kind=c_int32_t) :: xoshiro256ss_next_int31
         end function
      end interface

      xvalue = xoshiro256ss_next_int31(this%st)
      next_int31 = xvalue
   end function next_int31

   subroutine jump(this)
      implicit none
      class(randstate) :: this
      interface
         subroutine xoshiro256ss_jump(st) bind(C)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: st
         end subroutine
      end interface

      call xoshiro256ss_jump(this%st)
   end subroutine

   subroutine long_jump(this)
      implicit none
      class(randstate) :: this
      interface
         subroutine xoshiro256ss_long_jump(st) bind(C)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: st
         end subroutine
      end interface

      call xoshiro256ss_long_jump(this%st)
   end subroutine
end module randgen
