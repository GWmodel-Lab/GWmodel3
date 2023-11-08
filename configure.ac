AC_INIT([GWmodel3], 0.2-2)

AC_LANG(C++)
AC_REQUIRE_CPP

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test "${GSL_CONFIG}" != ""; then
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_ARG_ENABLE([cuda], [AS_HELP_STRING([--enable-cuda],[install cuda (default no)])],
[ case "${enableval}" in
    yes) ENABLE_CUDA=1 ;;
    no)  ENABLE_CUDA=0 ;;
    *)   AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
esac])


if test ${ENABLE_CUDA} -eq 1; then
    AC_MSG_CHECKING("Checking environment variable CUDA_HOME")
    if test -z "${CUDA_HOME}"; then
        AC_MSG_RESULT("CUDA_HOME not set, stop finding nvcc")
    else
        AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
    fi


    if test -z "${CUDA_HOME}"; then
        if test ${ENABLE_CUDA} -eq 1; then
            echo "Please make sure your cuda environment is under the /usr/local directory"
            exit 1
        fi
    fi


    if test ${ENABLE_CUDA} -eq 1; then
        AC_MSG_CHECKING("whether this is the 64 bit linux version of CUDA")
        AC_CHECK_FILE([${CUDA_HOME}/lib64/libcudart.so], [BIT64="yes"])
        if test -z "${BIT64}"; then
            AC_MSG_RESULT("no -- using ${CUDA_HOME}/lib for CUDA libs")
            cu_libdir="${CUDA_HOME}/lib"
            AC_SUBST(CUDA_LIBDIR,["$cu_libdir"])
        else
            AC_MSG_RESULT("yes -- using ${CUDA_HOME}/lib64 for CUDA libs")
            cu_libdir="${CUDA_HOME}/lib64"
            AC_SUBST(CUDA_LIBDIR,["$cu_libdir"])
        fi
    else
    AC_SUBST(CUDA_LIBDIR,[""])
    fi
fi


if test ${ENABLE_CUDA} -eq 1; then
    AC_MSG_RESULT(You have CUDA environment and we will compile the cuda function for you)
    # cu_sources="BandwidthWeightKernel.cu CRSDistanceKernel.cu"
    # cu_sharedlibs="BandwidthWeightKernel.o CRSDistanceKernel.o"
    # cpp_sources_cuda="GWRBasicGpuTask.cpp IGWmodelCUDA.cpp"
    # cpp_sharedlibs_cuda="CGWmodelCUDA.o IGWmodelCUDA.o"
    # OBJECTS_CUDA="${cu_sharedlibs} ${cpp_sharedlibs_cuda}"

    CUDA_DEFS="-DENABLE_CUDA"
    CUDA_LIBS="-L${CUDA_LIBDIR} -Wl,-rpath=${CUDA_LIBDIR} -lcudart -lcublas "
    # AC_SUBST(cu_sources)
    # AC_SUBST(cu_sharedlibs)
    # AC_SUBST(cpp_sources_cuda)
    # AC_SUBST(cpp_sharedlibs_cuda)
    # AC_SUBST(OBJECTS_CUDA)
    AC_SUBST(CUDA_DEFS)
    AC_SUBST(TARGET_CUDA_LIB)
    AC_SUBST(CUDA_LIBS)
fi


AC_MSG_RESULT(building the nvcc command line)
NVCC="${CUDA_HOME}/bin/nvcc"
AC_MSG_RESULT(building the cuda include path)
CUDA_INCL="${CUDA_HOME}/include"


AC_SUBST(GWmodel_CUDA_LIB_PATH,["~/.cache/GWmodel"])


AC_MSG_NOTICE([Building Makevars])
AC_SUBST(CU_INCL)
AC_SUBST(NVCC)
AC_SUBST(PKG_LIBS)
AC_SUBST(CUDA_INCL)
AC_SUBST(ENABLE_CUDA)

AC_CONFIG_FILES([src/Makevars])
AC_OPENMP
AC_OUTPUT
