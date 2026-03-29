! -*- F90 -*-

module ermod_mpi_client
  implicit none
  integer(8) :: vp

contains
  subroutine connect(status)
    implicit none
    logical, intent(out) :: status
    external ermod_connect

    vp = 0
    call ermod_connect(vp)
    if(vp /= 0) then
       status = .false.
       return
    endif
    status = .true.
    return
  end subroutine connect

  subroutine disconnect()
    implicit none
    external ermod_disconnect
    
    call ermod_disconnect(vp)
  end subroutine disconnect
    
  subroutine send_trajectory(natoms, cell, coords, status)
    implicit none
    integer, intent(in) :: natoms
    real(8), intent(in) :: cell(3,3), coords(3, natoms)
    integer, intent(out) :: status

    external ermod_send_trajctory
    
    call ermod_send_trajctory(vp, natoms, cell, coords, status)
  end subroutine send_trajectory
end module ermod_mpi_client

