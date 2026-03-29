! -*- F90 -*-

! From Fortran 2003 and onwards, ISO_C_BINDING module is supported.
! The module allows safe interoperability to C codes.

module ermod_mpi_client_impl
  use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int
  interface
     type(c_ptr) function ermod_connect() bind(C)
       import
     end function ermod_connect

     subroutine ermod_disconnect(p) bind(C)
       import
       type(c_ptr), value, intent(in) :: p
     end subroutine ermod_disconnect
     
     integer(kind=c_int) function ermod_send_trajectory(p, natoms, cell, coords) bind(C)
       import
       type(c_ptr), value, intent(in) :: p
       integer(kind=c_int), value, intent(in) :: natoms
       real(kind=c_double), intent(in) :: cell(3,3)
       real(kind=c_double), intent(in) :: coords(3, natoms)
     end function ermod_send_trajectory
  end interface
end module ermod_mpi_client_impl

module ermod_mpi_client
  use, intrinsic :: iso_c_binding, only: c_ptr
  implicit none

  type(c_ptr) :: h
contains

  subroutine connect(status)
    use, intrinsic :: iso_c_binding, only: c_associated
    use ermod_mpi_client_impl
    implicit none
    logical, intent(out) :: status
    h = ermod_connect()
    status = c_associated(h)
    return
  end subroutine connect
  
  subroutine disconnect()
    use ermod_mpi_client_impl
    implicit none
    
    call ermod_disconnect(h)
  end subroutine disconnect
    
  subroutine send_trajectory(natoms, cell, coords, status)
    use, intrinsic :: iso_c_binding, only: c_loc
    use ermod_mpi_client_impl
    implicit none
    integer, intent(in) :: natoms
    real(8), intent(in) :: cell(3,3), coords(3, natoms)
    integer, intent(out) :: status
    integer(kind=c_int) :: tmp

    tmp = ermod_send_trajectory(h, natoms, cell, coords)
    status = tmp
  end subroutine send_trajectory

end module ermod_mpi_client

