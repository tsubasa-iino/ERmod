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

module realcal
   implicit none
   integer :: nsolu_atom, nsolv_atom
   integer, allocatable :: block_solu(:, :), block_solv(:, :)
   integer, allocatable :: belong_solu(:), belong_solv(:)
   integer, allocatable :: atomno_solu(:), atomno_solv(:)
   integer, allocatable :: counts_solu(:, :, :), counts_solv(:, :, :)
   integer, allocatable :: psum_solu(:), psum_solv(:)
   real, allocatable :: sitepos_solu(:, :), sitepos_solv(:, :)
   integer, allocatable :: ljtype_solu(:, :), ljtype_solv(:, :)

   integer :: max_solu_block, max_solv_block

   real, allocatable :: ljeps_lowlj(:, :), ljsgm2_lowlj(:, :), dist_lowlj(:, :)
   integer, allocatable :: belong_lowlj(:, :)
   real, allocatable :: ljeps_switch(:, :), ljsgm2_switch(:, :), dist_switch(:, :)
   integer, allocatable :: belong_switch(:, :)
   real, allocatable :: charge_el(:, :), dist_el(:, :)
   integer, allocatable :: belong_el(:, :)
   real, allocatable :: e_t(:, :)

   ! subcell_neighbour unit length and its number 
   real :: block_unit_axes(3), block_unit_axes_inv(3)
   integer :: block_size(3)

   ! "straight" coordinate system
   real, allocatable :: sitepos_normal(:, :)
   real :: cell_normal(3, 3), invcell_normal(3), cell_len_normal(3)
   logical :: is_cuboid
   real, parameter :: check_rotate = 1e-8, cuboid_thres = 1e-8, cutoff_thres = 1e-4

contains
   ! The big super subroutine for allocating - calculate - deallocating combo.
   ! After calling this routine, [uvengy] stores the interaction energy between solute and solvent.
   subroutine realcal_proc(target_solu, tagpt, slvmax, uvengy)
      use engmain, only: numsite
!$    use omp_lib, only: omp_get_max_threads
      implicit none
      integer, intent(in) :: target_solu, tagpt(:), slvmax
      real, intent(out) :: uvengy(0:slvmax)
      real, allocatable :: eng(:, :)
      integer :: lsize, i, j
      integer :: npar

      ! print *, "DEBUG: relcal_proc called"
      ! FIXME: fix calling convention & upstream call tree
      ! to calculate several solutes at once
      nsolu_atom = numsite(target_solu)
      nsolv_atom = count_solv(target_solu, tagpt, slvmax)

      allocate(block_solu(3, nsolu_atom), block_solv(3, nsolv_atom))
      allocate(belong_solu(nsolu_atom), belong_solv(nsolv_atom))
      allocate(atomno_solu(nsolu_atom), atomno_solv(nsolv_atom))
      allocate(counts_solu(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
      allocate(counts_solv(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))
      allocate(psum_solu(0:block_size(1) * block_size(2) *  block_size(3)))
      allocate(psum_solv(0:block_size(1) * block_size(2) *  block_size(3)))

      call set_solv_atoms(target_solu, tagpt, slvmax)
      call set_solu_atoms(target_solu)

      ! assertion
      ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged"

      call blockify(nsolu_atom, atomno_solu, block_solu)
      call blockify(nsolv_atom, atomno_solv, block_solv)

      ! We sort both blocks (solute, solvent). This is because the two close atoms will interact with similar atoms.
      ! We speedup the calculation using the access locality.
      call sort_block(block_solu, nsolu_atom, belong_solu, atomno_solu, counts_solu, psum_solu)
      call sort_block(block_solv, nsolv_atom, belong_solv, atomno_solv, counts_solv, psum_solv)

      allocate(sitepos_solv(3, nsolv_atom + 1))
      sitepos_solv(1:3, 1:nsolv_atom) = sitepos_normal(1:3, atomno_solv(1:nsolv_atom))
      sitepos_solv(1:3, nsolv_atom+1) = 0 ! sentinel

      max_solu_block = maxval(counts_solu)
      max_solv_block = 0
      do i = 0, block_size(2) - 1
         do j = 0, block_size(3) - 1
            max_solv_block = max(max_solv_block, sum(counts_solv(:, i, j)))
         end do
      end do
      lsize = max_solu_block * max_solv_block

      npar = 1
!$    npar = omp_get_max_threads()

      allocate(ljeps_lowlj(lsize, npar), ljsgm2_lowlj(lsize, npar), dist_lowlj(lsize, npar), belong_lowlj(lsize, npar))
      allocate(ljeps_switch(lsize, npar), ljsgm2_switch(lsize, npar), dist_switch(lsize, npar), belong_switch(lsize, npar))
      allocate(charge_el(lsize, npar), dist_el(lsize, npar), belong_el(lsize, npar))
      allocate(e_t(lsize, npar))
      ! assertion
      ! if (.not. all(belong_solu(:) == target_solu)) stop "realcal_blk: target_solu bugged after sorting"

      ! calculate the energy on parallel and sum up after that.
      allocate(eng(1:slvmax, npar))
      eng(:, :) = 0.0
      call get_pair_energy(eng)

      uvengy(1:slvmax) = sum(eng(1:slvmax, 1:npar), 2)

      deallocate(ljeps_lowlj, ljsgm2_lowlj, dist_lowlj, belong_lowlj)
      deallocate(ljeps_switch, ljsgm2_switch, dist_switch, belong_switch)
      deallocate(charge_el, dist_el, belong_el)
      deallocate(e_t)

      deallocate(eng)
      deallocate(block_solu, belong_solu, atomno_solu, counts_solu, psum_solu)
      deallocate(block_solv, belong_solv, atomno_solv, counts_solv, psum_solv)
      deallocate(sitepos_solv)
   end subroutine realcal_proc

   subroutine realcal_prepare
      use engmain, only: numatm, sitepos, boxshp, SYS_PERIODIC
      implicit none

      allocate(sitepos_normal(3, numatm))
      sitepos_normal(:, :) = sitepos(:, :)

      if(boxshp == SYS_PERIODIC) then
         ! "Straighten" and normalize cell
         call normalize_cell

         ! normalize coordinate system, put atoms in 0 <= r < L
         call normalize_coordinate

         ! In realcal, first system are split into blocks, then calculation is perfomed block-wise.
         call set_block_info
      endif
   end subroutine realcal_prepare

   subroutine realcal_cleanup
      implicit none
      deallocate(sitepos_normal)
   end subroutine realcal_cleanup


   ! Calculate i-j interaction energy in the bare 1/r form
   function realcal_bare(i, j) result(pairep)
      use engmain, only:  boxshp, numsite, &
         elecut, lwljcut, upljcut, cltype, screen, charge, specatm, &
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, SYS_PERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
      implicit none
      integer, intent(in) :: i, j
      real :: pairep
      integer :: is, js, ismax, jsmax, ati, atj
      real :: reelcut, rst, dis2, invr2, invr3, invr6
      real :: eplj, epcl, xst(3), half_cell(3)
      real :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
      real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
      real :: repA, repB, repC, attA, attB, attC
      integer :: ljtype_i, ljtype_j
      !
      if(i == j) stop "cannot happen: two particle arguments should not be the same"
      if(cltype /= EL_COULOMB) stop "cannot happen: realcal_bare is called only when cltype is 'bare coulomb'."

      if(boxshp == SYS_NONPERIODIC) reelcut=huge(reelcut)
      if(boxshp == SYS_PERIODIC) then
         reelcut = elecut
         half_cell(:) = 0.5 * cell_len_normal(:)
      else
         ! suppress warnings
         half_cell(:) = 0.0
      endif

      pairep = 0.0
      ismax = numsite(i)
      jsmax = numsite(j)

      if(ljswitch == LJSWT_FRC_CHM) then       ! force switch (CHARMM type)
         lwljcut3 = lwljcut ** 3
         upljcut3 = upljcut ** 3
         lwljcut6 = lwljcut3 * lwljcut3
         upljcut6 = upljcut3 * upljcut3
      endif
      if(ljswitch == LJSWT_FRC_GMX) then       ! force switch (GROMACS type)
         call calc_gmx_switching_force_params(12, lwljcut, upljcut, repA, repB, repC)
         call calc_gmx_switching_force_params(6,  lwljcut, upljcut, attA, attB, attC)
      endif

      do is = 1, ismax
         do js = 1, jsmax
            ati = specatm(is,i)
            atj = specatm(js,j)
            ljtype_i = ljtype(ati)
            ljtype_j = ljtype(atj)
            xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
            if(boxshp == SYS_PERIODIC) then    ! when the system is periodic
               if(is_cuboid) then
                  xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
               else
                  ! Note some ops can be skipped because cell_normal is upper triangular
                  xst(:) = xst(:) - cell_normal(:, 3) * anint(xst(3) * invcell_normal(3))
                  xst(:) = xst(:) - cell_normal(:, 2) * anint(xst(2) * invcell_normal(2))
                  xst(:) = xst(:) - cell_normal(:, 1) * anint(xst(1) * invcell_normal(1))
               end if
            endif
            dis2 = sum(xst(1:3) ** 2)
            rst = sqrt(dis2)
            if(rst > upljcut) then
               eplj = 0.0
            else
               ljeps = ljene_mat(ljtype_i, ljtype_j)
               ljsgm2 = ljlensq_mat(ljtype_i, ljtype_j)

               invr2 = ljsgm2 / dis2
               invr6 = invr2 * invr2 * invr2
               select case(ljswitch)
                case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
                  eplj = 4.0 * ljeps * invr6 * (invr6 - 1.0)
                  if(rst > lwljcut) then
                     select case(ljswitch)
                      case(LJSWT_POT_CHM)                  ! CHARMM type
                        lwljcut2 = lwljcut ** 2
                        upljcut2 = upljcut ** 2
                        swth = (2.0 * dis2 + upljcut2 - 3.0 * lwljcut2)      &
                           * ((dis2 - upljcut2) ** 2)                      &
                           / ((upljcut2 - lwljcut2) ** 3)
                      case(LJSWT_POT_GMX)                  ! GROMACS type
                        swfac = (rst - lwljcut) / (upljcut - lwljcut)
                        swth = 1.0 - 10.0 * (swfac ** 3)                     &
                           + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5)
                      case default
                        stop "Unknown ljswitch"
                     end select
                     eplj = swth * eplj
                  endif
                case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
                  ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                  if(rst <= lwljcut) then
                     vdwa = invr6 * invr6 - ljsgm6 *ljsgm6 / (lwljcut6 * upljcut6)
                     vdwb = invr6 - ljsgm6 / (lwljcut3 * upljcut3)
                  else
                     invr3 = sqrt(invr6)
                     ljsgm3 = sqrt(ljsgm6)
                     vdwa = upljcut6 / (upljcut6 - lwljcut6)                 &
                        * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
                     vdwb = upljcut3 / (upljcut3 - lwljcut3)                 &
                        * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
                  endif
                  eplj = 4.0 * ljeps * (vdwa - vdwb)
                case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
                  ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
                  if(rst <= lwljcut) then
                     vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * repC
                     vdwb = invr6 - ljsgm6 * attC
                  else
                     swfac = rst - lwljcut
                     vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 *                &
                        (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
                     vdwb = invr6 - ljsgm6 *                                 &
                        (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
                  endif
                  eplj = 4.0 * ljeps * (vdwa - vdwb)
                case default
                  stop "Unknown ljswitch"
               end select
            endif
            if(rst >= reelcut) then
               epcl = 0.0
            else
               epcl = charge(ati) * charge(atj) / rst
            endif
            pairep = pairep + eplj + epcl
         end do
      end do
   end function realcal_bare

   ! self-energy part, no LJ calculation performed
   function realcal_self(i) result(pairep)
      use engmain, only: numsite, screen, cltype, charge, specatm, &
         EL_COULOMB, PI8
      implicit none
      integer, intent(in) :: i
      real(kind=8) :: pairep
      integer :: is, js, ismax, ati, atj
      real(kind=8) :: chgi, epcl, rst, dis2
      real :: xst(3), half_cell(3)

      pairep = 0.0
      if(cltype == EL_COULOMB) return

      half_cell(:) = 0.5 * cell_len_normal(:)

      ismax=numsite(i)

      !$omp parallel do private(is,js,ati,atj,epcl,rst,dis2,xst,chgi) reduction(+:pairep)
      do is = 1, ismax
         ati = specatm(is, i)
         ! When solute molecule is huge (> ~10,000 atoms), the error of cumulative addition of pairep would be large for the single-precision build.
         ! To reduce the cumulative error, chgi is declared as double precision for the first variable in the cumulative summation
         ! so that the summation is performed in double precision.
         chgi = charge(ati)

         ! Atom residual
         ! self (the same ati arguments for two charge variables below)
         epcl = - chgi * chgi * screen / sqrt(PI8)
         pairep = pairep + epcl

         do js = is + 1, ismax
            atj = specatm(js, i)

            xst(:) = sitepos_normal(:,ati) - sitepos_normal(:,atj)
            if(is_cuboid) then
               xst(:) = half_cell(:) - abs(half_cell(:) - abs(xst(:)))
            else
               ! Note some ops can be skipped because cell_normal is upper triangular
               xst(:) = xst(:) - cell_normal(:, 3) * anint(xst(3) * invcell_normal(3))
               xst(:) = xst(:) - cell_normal(:, 2) * anint(xst(2) * invcell_normal(2))
               xst(:) = xst(:) - cell_normal(:, 1) * anint(xst(1) * invcell_normal(1))
            end if

            ! real(..., 8) is to reduce the error of cumulative addition for the single-precision build.
            dis2 = sum( real(xst(1:3), 8) ** 2)
            rst = sqrt(dis2)

            ! distinct (different ati and atj arguments for two charge variables)
            epcl = - chgi * charge(atj) * erf(screen * rst) / rst
            pairep = pairep + epcl
         enddo
      enddo

   end function realcal_self


   integer function count_solv(solu, tagpt, slvmax)
      use engmain, only: numsite
      implicit none
      integer, intent(in) :: solu, tagpt(:), slvmax
      integer :: i, j, cnt
      cnt = 0
      do i = 1, slvmax
         j = tagpt(i)
         if(j == solu) cycle
         cnt = cnt + numsite(j)
      end do
      count_solv = cnt
   end function count_solv

   subroutine set_solu_atoms(solu)
      use engmain, only: numsite, mol_begin_index
      implicit none
      integer, intent(in) :: solu
      integer :: i
      do i = 1, numsite(solu)
         atomno_solu(i) = mol_begin_index(solu) + (i - 1)
         belong_solu(i) = solu
      end do
   end subroutine set_solu_atoms

   subroutine set_solv_atoms(solu, tagpt, slvmax)
      use engmain, only: numsite, mol_begin_index
      implicit none
      integer, intent(in) :: solu, tagpt(:), slvmax
      integer :: i, j, k, cnt
      cnt = 1
      do i = 1, slvmax
         j = tagpt(i)
         if(j == solu) cycle
         do k = 1, numsite(j)
            atomno_solv(cnt) = mol_begin_index(j) + (k - 1)
            belong_solv(cnt) = i
            cnt = cnt + 1
         end do
      end do
   end subroutine set_solv_atoms

   subroutine set_block_info
      use engmain, only: block_threshold
      implicit none

      ! get the length of axes
      ! assumes cell's 1st axis being x-axis, 2nd axis on x-y plane
      ! set block size
      block_size(:) = ceiling(cell_len_normal(:) / block_threshold)
      block_unit_axes(:) = cell_len_normal(:) / block_size(:)
      block_unit_axes_inv(:) = 1. / block_unit_axes(:)
   end subroutine set_block_info

   subroutine blockify(natom, atomlist, blk)
      implicit none
      integer, intent(in) :: natom, atomlist(:)
      integer, intent(out) :: blk(:, :)
      integer :: i, j, a

      do i = 1, natom
         a = atomlist(i)
         blk(:, i) = floor(sitepos_normal(:, a) * block_unit_axes_inv(:))
         blk(:, i) = max(0, blk(:, i)) ! In the worst case, due to round-off error it can be < 0
         blk(:, i) = min(block_size(:) - 1, blk(:, i)) ! similar
      end do
   end subroutine blockify

   subroutine sort_block(blk, nmol, belong, atomno, counts, psum)
      implicit none
      integer, intent(inout) :: blk(:, :)
      integer, intent(in) :: nmol
      integer, intent(inout) :: belong(:)
      integer, intent(inout) :: atomno(:)
      integer, intent(inout) :: counts(0:block_size(1) - 1, 0:block_size(2) - 1, 0:block_size(3) - 1)
      integer, intent(out) :: psum(0:block_size(1) * block_size(2) * block_size(3))
      integer, allocatable :: buffer(:, :) ! FIXME: ugly!
      integer, allocatable :: pnum(:, :, :)
      integer :: a, b, c
      integer :: i, j, k, partialsum, pos

      counts(:, :, :) = 0

      do i = 1, nmol
         a = blk(1, i)
         b = blk(2, i)
         c = blk(3, i)
         if(a < 0 .or. b < 0 .or. c < 0) STOP "INVL"
         counts(a, b, c) = counts(a, b, c) + 1
      end do

      allocate(pnum(0:block_size(1)-1, 0:block_size(2)-1, 0:block_size(3)-1))

      partialsum = 0
      do k = 0, block_size(3) - 1
         do j = 0, block_size(2) - 1
            do i = 0, block_size(1) - 1
               pnum(i, j, k) = partialsum
               partialsum = partialsum + counts(i, j, k)
            end do
         end do
      end do

      allocate(buffer(5, nmol))
      do i = 1, nmol
         a = blk(1, i)
         b = blk(2, i)
         c = blk(3, i)
         pos = pnum(a, b, c) + 1  ! buffer (and output) index starts from 1
         pnum(a, b, c) = pos
         buffer(1:3, pos) = blk(1:3, i)
         buffer(4, pos) = belong(i)
         buffer(5, pos) = atomno(i)
      end do
      deallocate(pnum)

      blk(1:3, :) = buffer(1:3, :)
      belong(:) = buffer(4, :)
      atomno(:) = buffer(5, :)

      deallocate(buffer)

      partialsum = 0
      pos = 0
      do k = 0, block_size(3) - 1
         do j = 0, block_size(2) - 1
            do i = 0, block_size(1) - 1
               psum(pos) = partialsum + 1 ! output index starts from 1
               pos = pos + 1
               partialsum = partialsum + counts(i, j, k)
            end do
         end do
      end do
      psum(pos) = partialsum + 1
      ! print *, psum
   end subroutine sort_block

   ! Get all-to-all pair interaction, using order of blocks.
   subroutine get_pair_energy(energy_vec)
      ! calculate for each subcell
      ! cut-off by subcell distance
      use engmain, only: upljcut, elecut
      implicit none
      real, intent(inout) :: energy_vec(:, :)
      integer :: u1, u2, u3
      integer :: vbs(3)
      integer :: upos, vpos_base, vpos_line_end, vpos_begin, vpos_end
      integer :: pbc_solushift(3)
      integer :: psx, psy, psz, bpx, bpy, bpz
      integer :: zmin, zmax, ymin, ymax, xmin, xmax
      real :: p(3), cutxyz, cutxy2, cutx2
      !$omp parallel &
      !$omp   private(u1, u2, u3, upos, vpos_begin, vpos_end, pbc_solushift) &
      !$omp   private(p, cutxyz, cutxy2, cutx2) &
      !$omp   private(psx, psy, psz, bpx, bpy, bpz, zmin, zmax, ymin, ymax, xmin, xmax) &
      !$omp   shared(energy_vec)
      !$omp single
      do psz = -1, 1
         pbc_solushift(3) = psz
         do psy = -1, 1
            pbc_solushift(2) = psy
            do psx = -1, 1
               pbc_solushift(1) = psx
               do u3 = 0, block_size(3) - 1
                  do u2 = 0, block_size(2) - 1
                     do u1 = 0, block_size(1) - 1
                        ! At this moment solute block (u1, u2, u3) is considered.
                        upos = u1 + block_size(1) * (u2 + block_size(2) * u3)
                        if(psum_solu(upos + 1) == psum_solu(upos)) cycle

                        !$omp task
                        ! find interacting solvent box, with solute shifted by (pvx, pvy, pvz) <dot> cell vector
                        ! note this interacting solvent box may target the same box twice;
                        ! but the cutoff operation shall prevent the double counting.
                        ! leftbottom-like position
                        p(:) = matmul(cell_normal, real(pbc_solushift)) + (/ u1, u2, u3 /) * block_unit_axes(:)
                        ! At this moment p(:) is periodic image (of 27).

                        cutxyz = max(upljcut, elecut)
                        zmin = floor((p(3) - cutxyz) * block_unit_axes_inv(3))
                        zmax = floor((p(3) + block_unit_axes(3) + cutxyz) * block_unit_axes_inv(3))
                        zmin = max(zmin, 0)
                        zmax = min(zmax, block_size(3) - 1)
                        do bpz = zmin, zmax
                           cutxy2 = cutxyz ** 2 - &
                              range_distance(bpz * block_unit_axes(3), &
                                             (bpz + 1) * block_unit_axes(3), &
                                             p(3), p(3) + block_unit_axes(3)) ** 2
                           if(cutxy2 <= 0) cycle
                           ymin = floor((p(2) - sqrt(cutxy2)) * block_unit_axes_inv(2))
                           ymax = floor((p(2) + block_unit_axes(2) + sqrt(cutxy2)) * block_unit_axes_inv(2))
                           ymin = max(ymin, 0)
                           ymax = min(ymax, block_size(2) - 1)
                           do bpy = ymin, ymax
                              cutx2 = cutxy2 - &
                                 range_distance(bpy * block_unit_axes(2), &
                                                (bpy + 1) * block_unit_axes(2), &
                                                p(2), p(2) + block_unit_axes(2)) ** 2
                              if(cutx2 <= 0) cycle
                              xmin = floor((p(1) - sqrt(cutx2)) * block_unit_axes_inv(1))
                              xmax = floor((p(1) + block_unit_axes(1) + sqrt(cutx2)) * block_unit_axes_inv(1))
                              xmin = max(xmin, 0)
                              xmax = min(xmax, block_size(1) - 1)
                              vpos_begin = xmin + block_size(1) * (bpy + block_size(2) * bpz)
                              vpos_end = xmax + block_size(1) * (bpy + block_size(2) * bpz) + 1 ! [vpos_begin, vpos_end)
                              if(xmin > xmax) cycle
                              call get_pair_energy_block(upos, vpos_begin, vpos_end, energy_vec, pbc_solushift)
                           end do
                        end do
                        !$omp end task
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end single
      !$omp end parallel

      contains
         function range_distance(mina, maxa, minb, maxb)
            real :: range_distance
            real, intent(in) :: mina, maxa, minb, maxb
            real :: minaa, maxaa, minbb, maxbb
            real :: temp
            
            if(mina <= minb) then
               minaa = mina; maxaa = maxa
               minbb = minb; maxbb = maxb
            else
               minaa = minb; maxaa = maxb
               minbb = mina; maxbb = maxa
            end if
            ! Now minaa is the minimum. There can be 3!=6 patterns (ignoring equality):
            ! 1. minaa < minbb < maxaa < maxbb
            ! 2. minaa < minbb < maxbb < maxaa
            ! 3. minaa < maxaa < minbb < maxbb
            ! 4. minaa < maxaa < maxbb < minbb    but impossible due to minbb < maxbb
            ! 5. minaa < maxbb < minbb < maxaa    but impossible due to minbb < maxbb
            ! 6. minaa < maxbb < maxaa < minbb    but` impossible due to minbb < maxbb
            if(minbb <= maxaa .and. maxaa <= maxbb) then
               ! case 1, overlapping
               range_distance = 0
            elseif(minbb <= maxbb .and. maxbb <= maxaa) then
               ! case 2, included
               range_distance = 0
            else
               range_distance = minbb - maxaa
            endif
         end function range_distance
   end subroutine get_pair_energy

   ! Computational kernel to calculate pair-energy between particles between solute-block and contiguous solvent-blocks
   subroutine get_pair_energy_block(upos, vpos_b, vpos_e, energy_vec, pbcushift)
      use engmain, only: cltype, boxshp, &
         upljcut, lwljcut, elecut, screen, charge,&
         ljswitch, ljtype, ljtype_max, ljene_mat, ljlensq_mat, &
         SYS_NONPERIODIC, EL_COULOMB, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX
!$    use omp_lib, only: omp_get_thread_num
      implicit none
      integer, intent(in) :: upos, vpos_b, vpos_e, pbcushift(3)
      real, intent(inout) :: energy_vec(:, :)
      integer :: ui, vi, ua, va, i, curp
      integer :: n_lowlj, n_switch, n_el
      integer :: belong_u, belong_v, ljtype_u, ljtype_v
      real :: crdu(3), crdv(3), d(3), dist, r, dist_next, invr2, invr3, invr6
      real :: lwljcut2, upljcut2, lwljcut3, upljcut3, lwljcut6, upljcut6
      real :: ljeps, ljsgm2, ljsgm3, ljsgm6, vdwa, vdwb, swfac
      real :: repA, repB, repC, attA, attB, attC
      real :: elecut2, half_cell(3)
      real :: pbcushift_real(3)

      if(cltype == EL_COULOMB) stop "realcal%get_pair_energy_block: cltype assertion failure"
      if(boxshp == SYS_NONPERIODIC) stop "realcal%get_pair_energy_block: boxshp assertion failure"

      curp = 1
!$    curp = omp_get_thread_num() + 1

      half_cell(:) = 0.5 * cell_len_normal(:)

      n_lowlj = 0
      n_switch = 0
      n_el = 0

      lwljcut2 = lwljcut ** 2
      upljcut2 = upljcut ** 2
      if(ljswitch == LJSWT_FRC_CHM) then       ! force switch (CHARMM type)
         lwljcut3 = lwljcut ** 3
         upljcut3 = upljcut ** 3
         lwljcut6 = lwljcut3 * lwljcut3
         upljcut6 = upljcut3 * upljcut3
      endif
      if(ljswitch == LJSWT_FRC_GMX) then       ! force switch (GROMACS type)
         call calc_gmx_switching_force_params(12, lwljcut, upljcut, repA, repB, repC)
         call calc_gmx_switching_force_params(6,  lwljcut, upljcut, attA, attB, attC)
      endif

      elecut2 = elecut ** 2
      pbcushift_real(:) = matmul(cell_normal, real(pbcushift))

      ! Here we first calculate and list distance between particles in the first loop,
      ! then calculate actual energy in the second loop.
      ! This was faster when I first wrote this code, but for any modern compiler + architecture
      ! I am not sure whether it's needed.
      ! TODO optimize:
      ! if you sort / reorder atomno, ljtype etc.
      ! this loop can be vectorized
      ! TODO optimize:
      ! software pipelining
      do ui = psum_solu(upos), psum_solu(upos + 1) - 1
         ua = atomno_solu(ui)
         belong_u = belong_solu(ui) ! FIXME: not used in later calculation
         ljtype_u = ljtype(ua)
         crdu(:) = sitepos_normal(:, ua)

         ! hide latency by calculating distance of next coordinate set
         d(:) = crdu(:) - sitepos_solv(:, psum_solv(vpos_b)) + pbcushift_real(:)

         dist_next = sum(d(:) ** 2)

         do vi = psum_solv(vpos_b), psum_solv(vpos_e) - 1
            va = atomno_solv(vi)
            belong_v = belong_solv(vi)
            ljtype_v = ljtype(va)

            crdv(:) = sitepos_solv(:, vi + 1)

            d(:) = crdu(:) - crdv(:) + pbcushift_real(:)

            ! assumes that only a single image matters for both electrostatic and LJ.
            ! if the box is very small and strongly anisotropic,
            ! there is a risk that second nearest image still being inside the cutoff length.
            ! But it's not considered in this case ...

            dist = dist_next
            dist_next = sum(d(:) ** 2) ! CHECK: any sane compiler will expand and unroll

            ! lines up all variables, to enable vectorization in 2nd phase
            ljeps = ljene_mat(ljtype_v, ljtype_u)
            if(ljeps > 0) then
               if(dist <= lwljcut2) then
                  n_lowlj = n_lowlj + 1
                  ljeps_lowlj (n_lowlj, curp) = ljeps
                  ljsgm2_lowlj(n_lowlj, curp) = ljlensq_mat(ljtype_v, ljtype_u)
                  dist_lowlj(n_lowlj, curp) = dist
                  belong_lowlj(n_lowlj, curp) = belong_v
               elseif(dist <= upljcut2) then
                  n_switch = n_switch + 1
                  ljeps_switch (n_switch, curp) = ljeps
                  ljsgm2_switch(n_switch, curp) = ljlensq_mat(ljtype_v, ljtype_u)
                  dist_switch(n_switch, curp) = dist
                  belong_switch(n_switch, curp) = belong_v
               end if
            end if

            if(dist <= elecut2) then
               n_el = n_el + 1
               charge_el(n_el, curp) = charge(ua) * charge(va)
               dist_el(n_el, curp) = dist
               belong_el(n_el, curp) = belong_v
            end if
         end do
      end do

      ! 2nd phase: calculate actual values with vectorized loop

      ! TODO optimize:
      ! explicit vectorization may increase performance

      ! LJ inside low cutoff
      do i = 1, n_lowlj
         ljeps = ljeps_lowlj(i, curp)
         ljsgm2 = ljsgm2_lowlj(i, curp)
         dist = dist_lowlj(i, curp)
         invr2 = ljsgm2 / dist
         invr6 = invr2 * invr2 * invr2
         select case(ljswitch)
          case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
            e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)
          case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
            ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
            vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 / (lwljcut6 * upljcut6)
            vdwb = invr6 - ljsgm6 / (lwljcut3 * upljcut3)
            e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
          case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
            ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
            vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * repC
            vdwb = invr6 - ljsgm6 * attC
            e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
          case default
            stop "Unknown ljswitch"
         end select
      end do
      do i = 1, n_lowlj
         energy_vec(belong_lowlj(i, curp), curp) = energy_vec(belong_lowlj(i, curp), curp) + e_t(i, curp)
      end do

      ! LJ switching region
      do i = 1, n_switch
         ljeps = ljeps_switch(i, curp)
         ljsgm2 = ljsgm2_switch(i, curp)
         dist = dist_switch(i, curp)
         invr2 = ljsgm2 / dist
         invr6 = invr2 * invr2 * invr2
         select case(ljswitch)
          case(LJSWT_POT_CHM)                   ! potential switch (CHRAMM type)
            e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)             &
               * (2.0 * dist + upljcut2 - 3.0 * lwljcut2)        &
               * ((dist - upljcut2) ** 2) / ((upljcut2 - lwljcut2) ** 3)
          case(LJSWT_POT_GMX)                   ! potential switch (GROMACS type)
            swfac = (sqrt(dist) - lwljcut) / (upljcut - lwljcut)
            e_t(i, curp) = 4.0 * ljeps * invr6 * (invr6 - 1.0)             &
               * (1.0 - 10.0 * (swfac ** 3)                      &
               + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5) )
          case(LJSWT_FRC_CHM)                   ! force switch (CHARMM type)
            invr3 = sqrt(invr6)
            ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
            ljsgm3 = sqrt(ljsgm6)
            vdwa = upljcut6 / (upljcut6 - lwljcut6)          &
               * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
            vdwb = upljcut3 / (upljcut3 - lwljcut3)          &
               * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
            e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
          case(LJSWT_FRC_GMX)                   ! force switch (GROMACS type)
            ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
            swfac = sqrt(dist) - lwljcut
            vdwa = invr6 * invr6 - ljsgm6 * ljsgm6 * &
               (repA * (swfac ** 3) + repB * (swfac ** 4) + repC)
            vdwb = invr6 - ljsgm6 * &
               (attA * (swfac ** 3) + attB * (swfac ** 4) + attC)
            e_t(i, curp) = 4.0 * ljeps * (vdwa - vdwb)
          case default
            stop "Unknown ljswitch"
         end select
      end do
      do i = 1, n_switch
         energy_vec(belong_switch(i, curp), curp) = energy_vec(belong_switch(i, curp), curp) + e_t(i, curp)
      end do

      ! ewald electrostatic
      do i = 1, n_el
         r = sqrt(dist_el(i, curp))
         e_t(i, curp) = charge_el(i, curp) * (1.0 - erf(screen * r)) / r
      end do
      do i = 1, n_el
         energy_vec(belong_el(i, curp), curp) = energy_vec(belong_el(i, curp), curp) + e_t(i, curp)
      end do
   end subroutine get_pair_energy_block

   ! Rotate box and coordinate so that the cell(:,:) is upper triangular:
   ! cell(2, 1) = cell(3, 1) = cell(3, 2) = 0
   ! (1st axis to be aligned to x-axis, 2nd axis to be within xy-plane)
   ! updates cell_normal and sitepos_normal
   subroutine make_cell_uppertriangular
      use engmain, only: cell, sitepos
      use mpiproc, only: warning
      implicit none

      real :: q(3,3), r(3, 3), temp(3, 3), axis_len
      integer :: i, j

      if (cell(2, 1) == 0.0 .and. cell(3, 1) == 0.0 .and. cell(3, 2) == 0.0) then
         ! already upper triangular
         ! As long as we use modern MD software the program shall go this path.
         cell_normal(:, :) = cell(:, :)
         sitepos_normal(:, :) = sitepos(:, :)
         return
      end if

      call warning("cell")

      ! modified Gram-Schmidt
      q(:, :) = 0.
      r(:, :) = 0.
      temp(:, :) = cell(:, :)
      do i = 1, 3
         axis_len = sqrt(dot_product(temp(:, 1), temp(:, 1)))
         q(:, i) = temp(:, i) / axis_len
         r(i, i) = axis_len
         do j = i + 1, 3
            ! remove inner products
            r(i, j) = dot_product(q(:, i), temp(:, j))
            temp(:, j) = temp(:, j) - r(i, j) * q(:, i)
         end do
      end do

      ! Now cell = QR
      ! cell_normal = R = Q^T Q R
      cell_normal(:, :) = r(:, :)

      ! sitepos_normal = Q^T X
      do i = 1, size(sitepos, 2)
         sitepos_normal(:, i) = matmul(transpose(q), sitepos(:, i))
      end do
   end subroutine make_cell_uppertriangular

   ! Cell vectors are normalized by either inverting the cell vector or by adding integer times other vectors
   subroutine normalize_cell_vector
      use engmain, only: elecut, upljcut
      use mpiproc, only: warning
      integer :: i
      real :: cutoff

      do i = 1, 3
         if (cell_normal(1, 1) < 0.0) then
            cell_normal(:, i) = -cell_normal(:, i)
         end if
      end do

      ! check cell size restrictions
      if (abs(cell_normal(1, 2)) > 0.5 * cell_normal(1, 1) + cutoff_thres) then
         call warning("cel2")
         cell_normal(:, 2) = cell_normal(:, 2) - &
            cell_normal(:, 1) * anint(cell_normal(1, 2) * invcell_normal(1))
      end if
      if (abs(cell_normal(1, 3)) > 0.5 * cell_normal(1, 1) + cutoff_thres) then
         call warning("cel2")
         cell_normal(:, 3) = cell_normal(:, 3) - &
            cell_normal(:, 1) * anint(cell_normal(1, 3) * invcell_normal(1))
      end if
      if (abs(cell_normal(2, 3)) > 0.5 * cell_normal(2, 2) + cutoff_thres) then
         call warning("cel2")
         cell_normal(:, 3) = cell_normal(:, 3) - &
            cell_normal(:, 2) * anint(cell_normal(2, 3) * invcell_normal(2))
      end if

      cutoff = min(elecut, upljcut)
      ! check cutoff restrictions
      do i = 1, 3
         if (cutoff > cell_normal(i, i) * 0.5) then
            stop "One of axis in periodic cell is too small. This is either the box is too small compared to the cell, " // &
               "or the periodic cell is too skewed."
         endif
      end do
   end subroutine normalize_cell_vector

   subroutine normalize_cell
      use engmain, only: cell, sitepos
      implicit none
      integer :: i

      call make_cell_uppertriangular

      call normalize_cell_vector

      do i = 1, 3
         cell_len_normal(i) = abs(cell_normal(i, i))
      end do
      invcell_normal(:) = 1 / cell_len_normal(:)

      if(abs(cell(1, 2)) > cuboid_thres .or. &
         abs(cell(1, 3)) > cuboid_thres .or. &
         abs(cell(2, 3)) > cuboid_thres ) then
         is_cuboid = .false.
      else
         is_cuboid = .true.
      end if
   end subroutine normalize_cell

   subroutine normalize_coordinate
      use engmain, only: sitepos
      implicit none
      integer :: n, i

      n = size(sitepos, 2)

      ! move all particles inside the cuboid spanned by (0 .. cell(1, 1)), (0 .. cell(2,2)), (0 .. cell(3,3)).
      ! Shift with the order Z -> Y -> X
      do i = 1, n
         ! TODO optimize: some cell_normals are guaranteed to be zero
         ! shift Z(&X,Y)
         sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 3) * floor(invcell_normal(3) * sitepos_normal(3, i))
         ! shift Y(&X)
         sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 2) * floor(invcell_normal(2) * sitepos_normal(2, i))
         ! shift X
         sitepos_normal(1:3, i) = sitepos_normal(1:3, i) - &
            cell_normal(:, 1) * floor(invcell_normal(1) * sitepos_normal(1, i))
         
         ! Note that due to the round-off error, there may be a case s.t. sitepos_normal(:, i) < 0 or > L.
      end do
   end subroutine normalize_coordinate

   ! get the coefficients for gromacs force switching
   subroutine calc_gmx_switching_force_params(pow, lwljcut, upljcut, coeffA, coeffB, coeffC)
      implicit none
      integer, intent(in) :: pow
      real, intent(in) :: lwljcut, upljcut
      real, intent(out) :: coeffA, coeffB, coeffC
      real :: dfljcut

      dfljcut = upljcut - lwljcut
      coeffA = - real(pow) * (real(pow + 4) * upljcut                   &
         - real(pow + 1) * lwljcut)                  &
         / ((upljcut ** (pow + 2)) * (dfljcut ** 2)) / 3.0
      coeffB =   real(pow) * (real(pow + 3) * upljcut                   &
         - real(pow + 1) * lwljcut)                  &
         / ((upljcut ** (pow + 2)) * (dfljcut ** 3)) / 4.0
      coeffC = 1.0 / (upljcut ** pow) - coeffA * (dfljcut ** 3)         &
         - coeffB * (dfljcut ** 4)
   end subroutine calc_gmx_switching_force_params

end module realcal
