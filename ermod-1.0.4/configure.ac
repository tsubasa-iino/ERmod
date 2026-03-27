#                                               -*- Autoconf -*-

# autoconf <= 2.64 has a bug in AC_OPENMP
AC_PREREQ([2.69])
#AC_INIT([ermod],[0.3.%rev%],[Matubayasi Laboratory])
AC_INIT([ermod],[1.0.4],[Matubayasi Laboratory])
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR([vmdfio.c])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug-symbol],
	AS_HELP_STRING([--disable-debug-symbol],
		[Disable debug symbol]))

AS_IF([test x$enable_debug_symbol = xno],
      [ac_cv_prog_cc_g=no
       ac_cv_prog_fc_g=no])

# Workaround Intel MPI inconsistency.
# Intel MPI has ifx but some old Intel MPI does not provide mpiifx.
# This leads to the inconsistent FC and MPIFC.
# To workaround, we search ifx only if mpiifx is available.
AC_CHECK_PROG(HAS_MPIIFX, mpiifx, yes, no)
AS_IF([test x"$HAS_MPIIFX" = xyes -o x$enable_mpi = xno],
      [ifx_or_empty=ifx; icx_or_empty=icx]
      [])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC([$icx_or_empty icc fcc gcc cc])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
# The order is: ifx (only if mpiifx exists), ifort, Fujitsu, XLf, GNU fortran.
# Removed f90/f95 as a fallback compiler because there are too many issues related to this fallback.
# If you want exotic fortran environment use FC= directly and/or contact devs.
AC_PROG_FC([$ifx_or_empty ifort frt xlf90 xlf gfortran])
AC_PROG_LN_S

# Check for buggy variant of ifort.
AS_IF([test x$FC = xifort],
      [version=`$FC --version`
       case "x$version" in # (
       *12.1.*):
         AC_MSG_ERROR([This program does not run correctly with Intel Compiler 12.1.x. Plese update the compiler to the newer version.])
       ;;
       esac
       ])

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

: ${enable_openmp=no}
AC_ARG_ENABLE([openmp],
              AS_HELP_STRING([--enable-openmp],
                             [use OpenMP (default = no)]))
AC_OPENMP

AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR(["Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around"])]
      [])

AC_ARG_ENABLE([vmdplugin],
	AS_HELP_STRING([--disable-vmdplugin],
		[disable VMD plugin]))

AC_DEFUN([MY_CHECK_BLAS], [
  AC_MSG_CHECKING([whether $1 is available])
  save_LIBS=$LIBS
  LIBS="$LIBS $2"
  AC_LINK_IFELSE([AC_LANG_CALL([], [sgemm])],
                 [blas_ok=yes], [blas_ok=no])
  LIBS=$save_LIBS
  AC_MSG_RESULT([$blas_ok])
])

AC_ARG_WITH([mkl],
  AS_HELP_STRING([--with-mkl],
    [Use Intel MKL for BLAS and FFT]))

AS_IF([test "x$with_mkl" != "x" -a "x$with_mkl" != "xyes"],
      [AC_MSG_FAILURE([Specifying version to --with-mkl=(version) is now deprecated because of Intel MKl update])])

AS_IF([test "x$with_mkl" != "x"],
  [case x`uname -m` in # (
  xx86_64) :
  mkl_plat=em64t
  mkl_platnew=intel64
  ;; # (
  xi386) :
  mkl_plat=32
  mkl_platnew=ia32
  AC_MSG_WARN([i386 platform is no longer actively maintained; if you have trouble please contact ERmod usergroup])
  ;; # (
  *) :
  ;;
  esac

  AS_IF([test "x$MKLROOT" = "x"],
    [AC_MSG_FAILURE([Could not get MKL version. Please set MKLROOT environment variable, or run "source /path/to/mkl/xxxx.x/env/vars.sh"])])
  mklroot=$MKLROOT
  AC_MSG_NOTICE([Using MKL installed at $MKLROOT])
  for dummy in 1; do # dummy loop to enable "break"
    with_blas="-qmkl=parallel"
    MY_CHECK_BLAS([-qmkl=parallel works (only works for Intel compilers)], [$with_blas])
    AS_IF([test "x$blas_ok" = "xyes"],
          [FCFLAGS="-qmkl=parallel $FCFLAGS"; LDFLAGS="-qmkl=parallel $LDFLAGS"; break])

    with_blas="-mkl=parallel" # Intel 2021 manpage says -mkl is deprecated, but 2021.2 does not accept -qmkl and only accepts -mkl (ifort's known bug)
    MY_CHECK_BLAS([-mkl=parallel works (only works for old Intel compilers)], [$with_blas])
    AS_IF([test "x$blas_ok" = "xyes"],
          [FCFLAGS="-mkl=parallel $FCFLAGS"; LDFLAGS="-mkl=parallel $LDFLAGS"; break])

    # check 2024-format
    with_blas="-L$mklroot/lib -lmkl_rt -lpthread -lm -ldl"
    MY_CHECK_BLAS([MKL 2024+ or LIB_PATH is set],
                  [$with_blas])
    AS_IF([test "x$blas_ok" = "xyes"],
          [FCFLAGS="-I$mklroot/include $FCFLAGS"; LDFLAGS="$with_blas $LDFLAGS"; with_lapack=""; break])

    # check 2017-format
    with_blas="-L$mklroot/lib/$mkl_platnew -lmkl_rt -lpthread -lm -ldl"
    MY_CHECK_BLAS([MKL 2017-2023],
                  [$with_blas])
    AS_IF([test "x$blas_ok" = "xyes"],
          [FCFLAGS="-I$mklroot/include $FCFLAGS"; LDFLAGS="$with_blas $LDFLAGS"; with_lapack=""; break])
    AC_MSG_FAILURE([Failed to find appropriate MKL. If you are using old mkl, please update])
  done
])

AC_ARG_WITH([fftw],
  AS_HELP_STRING([--with-fftw[[=(fftw root, default to /usr)]]],
    [Set fftw root directory]))

case x$with_fftw in # (
  xyes) :
  with_fftw=/usr
  ;; # (
  xno) :
  with_fftw=
  ;; # (
  *) :
  # do nothing
  ;;
esac

# FIXME: AC_SUBST should be used
AS_IF([test "x$with_fftw" != "x"],
  [FFTWROOT=$with_fftw
   CFLAGS="$CFLAGS -I$FFTWROOT/include"
   FCFLAGS="$FCFLAGS -I$FFTWROOT/include"
   LDFLAGS="$LDFLAGS -L$FFTWROOT/lib -L$FFTWROOT/lib64"])

ACX_BLAS
ACX_LAPACK([lapack_exist=yes],
           [lapack_exist=no])

################################

AC_ARG_ENABLE([double],
    [AS_HELP_STRING([--disable-double],
        [Disable conversion of single precision real to double precision (real(8))])],
    [disable_double=yes],
    [disable_double=no])

AC_DEFUN([MY_FC_REAL_DOUBLE],[
AC_CACHE_CHECK([how to read real as real(8) as default],
my_fc_real_double_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
# for default, gfortran, ifort, xlf, Fujitsu f90, Fujitsu f90 (old), Portland, Cray
for flag in "" "-fdefault-real-8 -fdefault-double-8" "-fdefault-real-8" "-autodouble" "-qrealsize=8" "-CcdRR8" "-Ad" "-Mr8" "-s default64"; do
FCFLAGS="$save_FCFLAGS $flag"
# This is very tricky indeed...
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [dnl
      implicit none
      real :: a
      integer, parameter :: len = kind(a)
      integer, parameter :: x(len) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
      print *, "Success"
])], [my_check_realasdouble=yes], [my_check_realasdouble=no])
  if test $my_check_realasdouble = yes; then
    break
  fi
done
if test $my_check_realasdouble = no; then
  AC_MSG_ERROR([unable to find compiler flag to read real as real(8)])
fi
FCFLAGS="$save_FCFLAGS $flag"
if test "x$flag" = x; then
  echo -n "(none)"
else
  echo -n $flag
fi
AC_LANG_POP(Fortran)
])])

AS_IF([test "x$disable_double" = "xno"],
	    [MY_FC_REAL_DOUBLE
	     FCFLAGS="$FCFLAGS -DDP"])

################################

AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether Fortran compiler supports transfer()])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none
      integer(8) :: a = 1
      real(8) :: b = 1.0d0
      
      a = transfer(b, a)
      print *, a
])], [my_check_have_transfer=yes], [my_check_have_transfer=no])

AC_MSG_RESULT([$my_check_have_transfer])
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_TRANSFER], [test "x$my_check_have_transfer" = "xyes"])

################################

AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether Fortran compiler supports flush()])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none

      open(10,file='/tmp/configure.ermod.flush.test')
      write(10,*) 342
      flush(10)
      close(10)
])], [my_check_have_flush=yes], [my_check_have_flush=no])

AC_MSG_RESULT([$my_check_have_flush])
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_FLUSH], [test "x$my_check_have_flush" = "xyes"])

################################

# set to f90 as default
AC_FC_SRCEXT(f90)

# select FFT library
AC_ARG_WITH([fft],
	AS_HELP_STRING([--with-fft=(fftw|mkl)],
		[specify which FFT library to be used]),
	[],
	[with_fft=no]
)

# Search for fftw3 - libfftw3.a
FFTWLIBS=
AC_CHECK_LIB([fftw3], [dfftw_execute], [have_fftw3=yes;FFTWLIBS=-lfftw3], [have_fftw3=no], [-lm])
have_fftw3_mp=no
AS_IF([test "$enable_openmp" != no],
  [AC_CHECK_LIB([fftw3_threads], [dfftw_init_threads], [have_fftw3_mp=yes;fftw3_mp=-lfftw3_threads], [], [-lfftw3 -lm])
   AC_CHECK_LIB([fftw3_omp], [dfftw_init_threads], [have_fftw3_mp=yes;fftw3_mp=-lfftw3_omp], [], [-lfftw3 -lm $OPENMP_FCFLAGS])]
   FFTWLIBS="$fftw3_mp $FFTWLIBS")
AC_SUBST([FFTWLIBS])

AS_IF([test "x$with_fft" == xfftw], 
            [AS_IF([test "x$have_fftw3" == xno],
                         [AC_MSG_FAILURE(
                                 [--with-fft=fftw is specified, but could not find fftw3])])])

AS_IF([test "x$with_fft" == xno],
	    [AS_IF([test "x$have_fftw3" == xyes],
	    		 [with_fft=fftw])
             AS_IF([test "x$with_mkl" != "x"],
                         [with_fft=mkl])])

AS_IF([test "x$with_fft" != xfftw && test "x$with_fft" != xmkl],
	    [AC_MSG_FAILURE(
		[invalid --with-fft=<fft-type> specification, please check with --help])])

# check whether mkl_dfti.f90 is accesible
# FIXME: This may fail to detect with gfortran 4.4.0-4.4.3
AS_IF([test "x$with_fft" = xmkl], [dnl
  AC_FC_SRCEXT(F90)
  AC_MSG_CHECKING([mkl_dfti.f90 is accessible])
  AC_COMPILE_IFELSE([#include "mkl_dfti.f90"
program main 
end program], [echo yes], [echo no
  AC_MSG_FAILURE([could not find mkl_dfti.f90 in include path])])
  AC_FC_SRCEXT(f90)])

AC_ARG_ENABLE([mpi],
	AS_HELP_STRING([--disable-mpi],
		[disable MPI even if exists]),
	AS_IF([test x$enable_mpi = xno],
		     [use_mpi=no]),
	[])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# FIXME: check cpp flags

AC_ARG_ENABLE([defaults],
	AS_HELP_STRING([--disable-defaults],
		[(deprcated)]))

AS_IF([test x$enable_defaults != xno],
CFLAGS="$CFLAGS"
FCFLAGS="$FCFLAGS")

AC_ARG_ENABLE([perf],
	AS_HELP_STRING([--enable-perf],
		[(internal use only)]),
	AS_IF([test x$enable_perf = xyes],
		     [FCFLAGS="$FCFLAGS -fno-omit-frame-pointer -DPERF"
                      CFLAGS="$CFLAGS -fno-omit-frame-pointer"
                      LDFLAGS="$LDFLAGS -lprofiler"]))

AM_CONDITIONAL(USE_MPI, [test "x$use_mpi" = "xyes"])

AS_IF([test "x$use_mpi" = "xyes"],
	    [mpi_result="enabled"
	     FC_CANONICAL="$FC"
	     FC="$MPIFC"
             FCFLAGS="-DMPI $FCFLAGS"],
	    [mpi_result="disabled"
	     FC_CANONICAL="$FC"])
AC_SUBST([FC_CANONICAL])

# Check FC_CANONICAL and FC are essentially the same compiler, otherwise confusing error message appears
save_FC=$FC
FC=$FC_CANONICAL
my_mpifc_fc_equiv=yes
AC_LANG_PUSH(Fortran)
FC=$FC_CANONICAL
AC_MSG_CHECKING([whether MPI- and non-MPI compilers have equivalent backends (non-MPI)])
AC_COMPILE_IFELSE([      program main
      end program], [], [my_mpifc_fc_equiv=no])
AC_MSG_RESULT([$my_mpifc_fc_equiv])
AC_LANG_POP(Fortran)
AC_LANG_PUSH(Fortran)
FC=$save_FC
AC_MSG_CHECKING([whether MPI- and non-MPI compilers have equivalent backends (MPI)])
AC_COMPILE_IFELSE([      program main
      end program], [], [my_mpifc_fc_equiv=no])
AC_MSG_RESULT([$my_mpifc_fc_equiv])
AC_LANG_POP(Fortran)
AS_IF([test $my_mpifc_fc_equiv = no], [
AC_MSG_ERROR([MPIFC($FC) and FC($FC_CANONICAL) seems to use different compiler backends. To fix this problem:
1. Please check whether these programs use the same compiler and the same version. For example, try "$FC --version" and "$FC_CANONICAL --version".
2. If these are not equal, find your way to specify the backend compiler for MPI (or change your compiler for non-MPI). Some examples are:
  OpenMPI: "export OMPI_FC=/path/to/compiler" (bash/zsh) or "setenv OMPI_FC /path/to/compiler" (csh/tcsh) will set the backend compiler of mpif90.
  MPICH: use e.g. ./configure MPIFC="mpif90 -fc=/path/to/compiler"
3. Run configure again with "./configure MPIFC=... FC=..." to explicitly specify the compiler.
])])

# optimization flags
AC_DEFUN([MY_FC_OPT_ACCEPT_NATIVE],[
AC_CACHE_CHECK([fortran compiler accepts -march=native],
my_fc_opt_accept_native_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
flag="-march=native"
FCFLAGS="$save_FCFLAGS $flag"
AC_COMPILE_IFELSE([      program main
      end program], [my_fc_opt_accept_native=yes], [my_fc_opt_accept_native=no])
echo -n $my_fc_opt_accept_native
FCFLAGS=$save_FCFLAGS
AC_LANG_POP(Fortran)
])])

AC_ARG_ENABLE([opt],
	AS_HELP_STRING([--disable-opt],
		[disable automatic edit to optimization flags]))

AS_IF([test x$enable_opt = xyes || test x$enable_opt = x],
  [FCFLAGS="$FCFLAGS"
   AS_IF([test x$ac_compiler_gnu = xyes],
        [optflag="-O3 -ffast-math -fno-finite-math-only"
         MY_FC_OPT_ACCEPT_NATIVE
         AS_IF([test x$my_fc_opt_accept_native = xyes],
                optflag="$optflag -march=native")],
        [case x$FC in #(
           xmpiifort|xifort|xmpiifx|xifx)
             optflag="-O3 -no-prec-div -xHost"
             ;; #(
           *)
             optflag="-O"
             ;;
         esac])
   FCFLAGS="$FCFLAGS $optflag"])

# Check OpenMP
AS_IF([test "$enable_openmp" != "no"],
  [FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
   # Check OpenMP version
   AC_MSG_CHECKING([whether OpenMP version >= 3.0])
   AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([], [dnl
!\$omp task
  print *, "test"
!\$omp end task
])], 
     [AC_MSG_RESULT([yes])], 
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([OpenMP must be >= 3.0])])])

# AM_CONDITIONAL([MPI], [test x$use_mpi = xyes])
AM_CONDITIONAL([USE_FFTW], [test x$with_fft = xfftw])
AM_CONDITIONAL([USE_MKL], [test x$with_fft = xmkl])
AM_CONDITIONAL([USE_VMDPLUGIN], [test x$enable_vmdplugin = xyes || test x$enable_vmdplugin = x])

AS_IF([test "x$lapack_exist" != xyes],
	    [AC_MSG_ERROR(
		[dnl
LAPACK library was not found.
Try specifying --with-lapack=<lib> options or set FCFLAGS/FLIBS environment variables to fix the problem.])])

lapack_result="$LAPACK_LIBS $BLAS_LIBS"

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AS_IF([test "$enable_openmp" != "no" && test "$with_fft" = "fftw" && test "$have_fftw3_mp" != "yes"],
  [AC_MSG_WARN([dnl
OpenMP is enabled, but your FFTW does not support multithreading.
This results in lower FFTW performance.
It is recommended to recompile FFTW library with --enable-omp or --enable-threads; see the FFTW manual for details.])])
  
case "x/$mklroot" in
  */2020.?.?/*|*/2021.?.?/*|*/2022.[[01]].?/*)
    AS_IF([test "x$enable_openmp" = "xyes" -a "x$with_fft" = "xmkl"],
          [AC_MSG_ERROR([MKL -2022.1.0's FFT is suspected to have problem with OpenMP. Please choose one from (1) disabling OpenMP (recommended), (2) using newer MKL, or (3) compiling FFTW then adding --with-fft=fftw to configure. Note current ERmod works the best without OpenMP.])])
  ;;
esac

AC_MSG_RESULT([dnl
LAPACK: $lapack_result
FFT: $with_fft
MPI: $mpi_result
OpenMP: $enable_openmp
FC(for erdst)=$FC
FC_CANONICAL(for slvfe)=$FC_CANONICAL
FCFLAGS=$FCFLAGS
LDFLAGS=$LDFLAGS
LIBS=$LIBS

Configuration complete; type "make" to start compiling])

