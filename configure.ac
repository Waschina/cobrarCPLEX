AC_INIT([cobrarCPLEX], [0.1.1], [s.waschina@nutrinf.uni-kiel.de])

FOUND_ALL="no"

# Check if global variables are already set, if so we are done
if test [-n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS"] ; then
   FOUND_ALL="yes"
else
   if test [-n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS"] ; then
      AC_MSG_ERROR([All of PKG_CPPFLAGS and PKG_LIBS must be given together])
   fi
fi

dnl get arguments
AC_ARG_WITH([cplex-dir],
    [AS_HELP_STRING([--with-cplex-dir],[specify location of CPLEX installation directory])],
    CPLEX_DIR=$withval,
    CPLEX_DIR="")


# Define default CPLEX path
DEFAULT_CPLEX_PATH="/opt/ibm/ILOG/CPLEX_Studio"

# Function to check for the latest CPLEX installation in the default path
check_latest_cplex_path() {
    latest_version_dir=$(ls -d ${DEFAULT_CPLEX_PATH}* 2>/dev/null | sort -V | tail -n 1)
    if test -d "$latest_version_dir/cplex/include" ; then
        CPLEX_DIR="$latest_version_dir/cplex"
        return 0
    fi
    return 1
}

# If no CPLEX_DIR is provided and not already set, check default path
if test ["$FOUND_ALL" = "no" -a -z "$CPLEX_DIR"] ; then
    check_latest_cplex_path
fi


if test ["$FOUND_ALL" = "no" -a -n "$CPLEX_DIR"] ; then
   PKG_CPPFLAGS="-I${CPLEX_DIR}/include -I${CPLEX_DIR}/include/ilcplex -I${CPLEX_DIR}/../concert/include"
   SYSTEM=`uname -m | tr '_' '-'`
   LIBFORMAT=`uname -s | tr A-Z a-z`
   if test ["$LIBFORMAT" = "darwin"]; then
        LIBFORMAT=osx
   fi
   PKG_LIBS="-L${CPLEX_DIR}/lib/${SYSTEM}_${LIBFORMAT}/static_pic -L${CPLEX_DIR}/../concert/lib/${SYSTEM}_${LIBFORMAT}/static_pic -lilocplex -lcplex -lconcert -lm -lpthread -ldl"

   FOUND_ALL="yes"
fi

if test "$FOUND_ALL" = "no" ; then
   AC_MSG_ERROR([No CPLEX installation found])
fi

AC_PROG_CXX
AC_PROG_CC

dnl substitute in src/Makevars
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

exit 0

