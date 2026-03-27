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

module sysvars
   implicit none

   character(len=5) :: clcond = 'merge'
   character(len=3) :: uvread = 'yes',    slfslt = 'yes',   ljlrc = 'not'
   character(len=3) :: infchk = 'not',    meshread = 'not', cumuint = 'not'
   character(len=3) :: write_mesherror = 'cnd'

   character(len=3) :: extsln = 'lin'    ! extrapolation of omega
   character(len=3) :: slncor = 'not'    ! energy correlation in solution
   character(len=3) :: refmerge = 'yes'
   character(len=3) :: readwgtfl = 'yes' ! reading weight_soln and weight_refs

   character(len=3) :: invmtrx = 'reg'   ! gce, reg, or evd
   character(len=4) :: zerosft = 'eczr'
   character(len=3) :: wrtzrsft = 'not'  ! writing the zero shifts
   character(len=4) :: wgtfnform = 'harm'
   character(len=3) :: wgtf2smpl = 'yes'

   character(len=3) :: normalize = 'yes', showdst= 'not'

   integer :: numprm = 0                           ! initialized to 0
   integer :: numprm_def_inf_yes = 11   ! default numprm at infchk = 'yes'
   integer :: numprm_def_inf_not = 5    ! default numprm at infchk = 'not'

   integer :: numsln = 0, numref = 0, numdiv = 0   ! initialized to 0
   integer :: maxsln, maxref, numrun, prmmax
   integer :: numslv, ermax

   real :: inptemp = 300.0              ! temperature in Kelvin, initialized
   real :: temp, kT, slfeng
   real :: avevolume = 0.0              ! average volume of system, initialized

   integer :: pickgr = 3
   integer :: msemin = 1, msemax = 5
   real :: mesherr = 0.1                ! allowed mesh error in kcal/mol

   integer :: extthres_soln = 1, extthres_refs = 1
   integer :: minthres_soln = 0, minthres_refs = 0

   ! iterative normalization of correlation matrix
   real :: norm_error = 1.0e-8          ! allowed error
   integer :: itrmax = 100              ! maximum number of iteration

   ! maximum numbers of lines of engsln and engref
   integer :: ermax_limit = 15000

   ! number of digits for suffixes of engsln, engref, and corref
   integer :: digits_of_suffix = 2

   real, parameter :: zero = 0.0
   real :: tiny = 1.0e-8
   integer :: large = 500000

   character(len=1024) :: solndirec = 'soln'
   character(len=1024) :: refsdirec = 'refs'
   character(len=1024) :: wgtslnfl  = 'weight_soln'
   character(len=1024) :: wgtreffl  = 'weight_refs'
   character(len=1024) :: slndnspf  = 'engsln'
   character(len=1024) :: slncorpf  = 'corsln'
   character(len=1024) :: refdnspf  = 'engref'
   character(len=1024) :: refcorpf  = 'corref'
   character(len=1024) :: aveuvfile = 'aveuv.tt'
   character(len=1024) :: engmeshfile = 'EngMesh'
   character(len=1024) :: cumuintfl = 'cumsfe'
   character(len=10), parameter :: numbers='0123456789'

   real, dimension(:),     allocatable :: nummol
   integer, dimension(:),  allocatable :: rduvmax, rduvcore
   real, dimension(:),     allocatable :: rdcrd, rddst, rddns
   real, dimension(:,:),   allocatable :: rdslc, rdcor
   integer, dimension(:),  allocatable :: rdspec
   real, dimension(:,:,:), allocatable :: chmpt
   real, dimension(:),     allocatable :: aveuv
   real, dimension(:,:),   allocatable :: uvene, blockuv
   integer, dimension(:),  allocatable :: svgrp, svinf
   real, dimension(:),     allocatable :: wgtsln, wgtref

   logical :: force_calculation = .false., strict_ewald_parameters = .false.
   logical :: check_parameters_er = .true.
   logical :: suffix_of_engsln_is_tt = .false., suffix_of_engref_is_tt = .false.

   namelist /fevars/ clcond, numprm, numsln, numref, numdiv, &
      uvread, slfslt, infchk, meshread, &
      refmerge, extsln, slncor, readwgtfl, &
      invmtrx, zerosft, wrtzrsft, wgtfnform, wgtf2smpl, &
      inptemp, ljlrc, avevolume, &
      pickgr, write_mesherror, msemin, msemax, mesherr, &
      normalize, showdst, &
      extthres_soln, extthres_refs, minthres_soln, minthres_refs, &
      solndirec, refsdirec, wgtslnfl, wgtreffl, &
      slndnspf, slncorpf, refdnspf, refcorpf, &
      aveuvfile, engmeshfile, cumuint, cumuintfl, &
      ermax_limit, norm_error, itrmax, tiny, large, &
      force_calculation, strict_ewald_parameters, check_parameters_er

contains

   subroutine init_sysvars
      implicit none
      character(len=*), parameter :: parmfname = 'parameters_fe'
      integer, parameter :: iounit = 191
      integer :: ioerr, i
      integer :: count_soln = 0, count_refs = 0

      open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
      if(ioerr /= 0) goto 99
      read(iounit, nml = fevars)
      close(iounit)
99    continue

      if(clcond == 'merge') then
         call check_tt_or_numeric_files(solndirec, slndnspf, suffix_of_engsln_is_tt, count_soln)
         call check_tt_or_numeric_files(refsdirec, refdnspf, suffix_of_engref_is_tt, count_refs)

         open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
         if(ioerr /= 0) goto 91
         read(iounit, nml = fevars)
         close(iounit)
91       continue

         if((numsln <= 0) .or. (numsln > count_soln)) numsln = count_soln
         if((numref <= 0) .or. (numref > count_refs)) numref = count_refs

         if((numdiv <= 0) .or. (numdiv >= numsln)) numdiv = numsln
         if(mod(numsln, numdiv) /= 0) then
            do i = numdiv + 1, numsln      ! find the larger and closest divisor
               if(mod(numsln, i) == 0) exit
            enddo
            numdiv = i
         endif
         if(refmerge == 'not') then        ! see subroutine datread for refmerge
            if(numdiv > numref) stop " With refmerge = 'not', numdiv needs to be not larger than numref"
            if(mod(numref, numdiv) /= 0) then
               write(6, "(A,i2,A,i2,A)") " Note: only ", numdiv * (numref / numdiv), &
               &" files out of ", numref, " engref and corref files prepared"
               numref = numdiv * (numref / numdiv)
            endif
         endif

         ! check consistency in parameters_er
         if( check_parameters_er ) call check_params()
      endif

      if(numprm <= 0) then                 ! default setting
         if(infchk == 'yes') then
            numprm = numprm_def_inf_yes    ! default numprm at infchk = 'yes'
         else
            numprm = numprm_def_inf_not    ! default numprm at infchk = 'not'
         endif
      endif

      if(pickgr < msemin) stop " Incorrect setting: pickgr < msemin not allowed"
      if(pickgr > msemax) stop " Incorrect setting: pickgr > msemax not allowed"
      if(pickgr > numprm) stop " Incorrect setting: pickgr > numprm not allowed"

      contains

      ! To set default values for suffix_of_eng*_is_tt, we check and count number of files
      subroutine check_tt_or_numeric_files(dirname, trunk, use_tt, nfiles)
         use iso_fortran_env, only : ERROR_UNIT
         implicit none
         character(len=*), intent(in) :: dirname, trunk
         logical, intent(out) :: use_tt
         integer, intent(out) :: nfiles
         integer, parameter :: sufmax = 99
         character(len=3) :: file_suf
         character(len=1024) :: opnfile
         integer :: count_suf
         logical :: file_exist

         ! First check tt file exists
         opnfile = trim(dirname) // '/' // trim(trunk) // "." // get_suffix(1, .true.)
         inquire(file = opnfile, exist = file_exist)
         if (file_exist) then
            ! .tt files exists. In this case .01-.99 files shall NOT exist
            do count_suf = 1, sufmax
               opnfile = trim(dirname) // '/' // trim(trunk) // "." // get_suffix(count_suf)
               inquire(file = opnfile, exist = file_exist)
               if (file_exist) then
                  write (ERROR_UNIT, *) trim(trunk) // ".tt is not supposed to coexist with " &
                     // trim(trunk) // ".01, " &
                     // trim(trunk) // ".02, ..."
                  stop "Inconsistent files"
               end if
            end do
            nfiles = 1
            use_tt = .true.
            return
         end if

         nfiles = 0
         use_tt = .false.
         ! .tt file does not exist. In this case return number of existing files with numeric suffixes
         do count_suf = 1, sufmax
            opnfile = trim(dirname) // '/' // trim(trunk) // '.' // get_suffix(count_suf)
            inquire(file = opnfile, exist = file_exist)
            if( file_exist ) then
               nfiles = count_suf
            else
               if (count_suf == 1) then
                  write (ERROR_UNIT, *) "Neither " // trim(trunk) // ".01 nor " // trim(trunk) &
                     // ".tt exists. Perhaps " &
                     // trim(dirname) // " part is not calculated yet?"
                  stop "File not found"
               endif
               exit ! exit from loop with the number of existing files
            endif
         enddo
      end subroutine

   end subroutine init_sysvars

   function get_suffix(n, suffix_is_tt) result(res)
      implicit none
      integer, intent(in) :: n
      logical, intent(in), optional :: suffix_is_tt
      logical :: suffix_is_tt_work
      character(len=digits_of_suffix) :: res

      if (present(suffix_is_tt)) then
         suffix_is_tt_work = suffix_is_tt
      else
         suffix_is_tt_work = .false.
      endif

      if (.not. suffix_is_tt_work) then
         res = zero_padded_str(n, digits_of_suffix)
      else
         res = repeat('t', digits_of_suffix)
      endif

      return
   end function

   ! TODO (shun): upon next major release move this to utility or utility_file module
   ! Returns zero-padded string of [n]
   function zero_padded_str(n, digits) result(res)
      implicit none
      integer, intent(in) :: n, digits
      character(len=digits) :: res
      character(len=25) :: formatbuf ! max (I0.x) len is assumed to be ceil(log10(2^64+1)) + 5

      if(digits <= 0) stop "Invalid args to zero_padded_str()"

      write (formatbuf, "(A4,I0,A1)") "(I0.", digits, ")"
      write (res, trim(formatbuf)) n
      return
   end function

   ! Check parameter consistency (issue #17)
   subroutine check_params
      use engmain, engmain_force_calculation=>force_calculation
      implicit none
      integer :: boxshp_s, estype_s, insposition_s, insstructure_s
      real :: lwreg_s, upreg_s, lwstr_s, upstr_s, inptemp_s
      integer :: ljformat_s, ljswitch_s, cmbrule_s
      real :: lwljcut_s, upljcut_s, elecut_s, screen_s, ewtoler_s
      integer :: splodr_s, cltype_s, ms1max_s, ms2max_s, ms3max_s
      logical :: inconsistent
      
      ! load both parameters_er and compare varibles
      call dummy_init              ! dummy initialization for soln parameters
      call init_params(solndirec)
      boxshp_s = boxshp ; estype_s = estype
      insposition_s = insposition ; insstructure_s = insstructure
      lwreg_s = lwreg ; upreg_s = upreg
      lwstr_s = lwstr ; upstr_s = upstr
      inptemp_s = inptemp
      ljformat_s = ljformat ; ljswitch_s = ljswitch ; cmbrule_s = cmbrule
      lwljcut_s = lwljcut ; upljcut_s = upljcut
      elecut_s = elecut ; screen_s = screen ; ewtoler_s = ewtoler
      splodr_s = splodr ; cltype_s = cltype
      ms1max_s = ms1max ; ms2max_s = ms2max ; ms3max_s = ms3max

      call dummy_init              ! dummy initialization for refs parameters
      call init_params(refsdirec)

      inconsistent = .false.

      if((boxshp /= boxshp_s) .or. (estype /= estype_s)) then
         print *, "box shape / system type inconsistent"
         inconsistent = .true.
      end if
      if((cltype /= cltype_s) .or. &
         (lwljcut /= lwljcut_s) .or. &
         (upljcut /= upljcut_s) .or. &
         (ljswitch /= ljswitch_s) .or. &
         (ljformat /= ljformat_s) .or. &
         (cmbrule /= cmbrule_s) .or. &
         (splodr /= splodr_s) .or. & ! splodr is technically ok, but different splodrs look *very* wrong...
         (elecut /= elecut_s)) then
         print *, "Nonbond calculation conditions (electrostatic, LJ) inconsistent"
         inconsistent = .true.
      end if
      if(((insposition /= insposition_s) .and. &
          ((insposition /= INSPOS_RANDOM) .and. (insposition /= INSPOS_NOCHANGE) .and. (insposition /= INSPOS_GAUSS))) .or. &
         ((insstructure /= insstructure_s) .and. (insstructure /= INSSTR_NOREJECT)) .or. &
         (lwreg /= lwreg_s) .or. &
         (upreg /= upreg_s) .or. &
         (lwstr /= lwstr_s) .or. &
         (upstr /= upstr_s)) then
         print *, "Insertion conditions inconsistent between solution and reference."
         print *, " Test particle insertion parameters are also used in the solution system to reject" // &
            " snapshots in the solution system. We thus need to use same parameters for:" // &
            "insposition, insstructure, lwreg, upreg, lwstr, and upstr"
         inconsistent = .true.
      end if

      if(inconsistent) then
         if(force_calculation) then
            print *, "Proceeding the calculation because force_calculation is set"
         else
            print *, "==>Aborting because solution and reference systems do not match"
            print *, "This typically mean you need to rerun the simulation with matching conditions"
            print *, 'If you dare to proceed, set "force_calculation=.true." in parameters_fe'
            stop "Stopping the calculation due to inconsistency"
         end if
      end if

      if((screen /= screen_s) .or. &
         (ms1max /= ms1max_s) .or. &
         (ms2max /= ms2max_s) .or. &
         (ms3max /= ms3max_s) .or. &
         (ewtoler /= ewtoler_s)) then
         print *, "Some of Ewald parameters are inconsistent"
         ! error message appears only when screen is set in parameters_er
         if((screen_s > zero) .or. (screen > zero)) then
            print *, "beta =", screen_s, "(soln) /", screen, "(refs)"
         endif
         print *, "grid =", ms1max_s, ",", ms2max_s, ",", ms3max_s, "(soln) /", &
            ms1max, ",", ms2max, ",", ms3max, "(refs)"
         ! ewtoler is not used in calcutions of interaction energies in erdst
         ! error message appears only when ewtoler is set in parameters_er
         if((ewtoler_s > zero) .or. (ewtoler > zero)) then
            print *, "Ewald tolerance =", ewtoler_s, "(soln)", ewtoler, "(refs)"
         endif
         print *, "This is mostly harmless and slvfe continues;"
         print *, " but keep in mind that the F.E. estimation may become unstable due to this difference"
         print *, " (if you cannot prevent this inconsistency, keep ewald tolerance low in the SIMULATION software)"
         if(strict_ewald_parameters) then
            stop "Aborting due to strict condition"
         end if
      end if

      contains
         subroutine dummy_init
             integer, parameter :: dummy_int = -9   ! dummy value, negative
             real, parameter :: dummy_real = -99    ! dummy value, negative
             boxshp = dummy_int ; estype = dummy_int
             insposition = dummy_int ; insstructure = dummy_int
             lwreg = dummy_real ; upreg = dummy_real
             lwstr = dummy_real ; upstr = dummy_real
             inptemp = dummy_real
             ljformat = dummy_int ; ljswitch = dummy_int ; cmbrule = dummy_int
             lwljcut = dummy_real ; upljcut = dummy_real
             elecut = dummy_real ; screen = dummy_real ; ewtoler = dummy_real
             splodr = dummy_int ; cltype = dummy_int
             ms1max = dummy_int ; ms2max = ms1max ; ms3max = ms1max
         end subroutine dummy_init
   end subroutine check_params
end module sysvars
