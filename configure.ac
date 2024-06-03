AC_INIT([cobrarCPLEX], [0.1.0], [s.waschina@nutrinf.uni-kiel.de])

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
AC_ARG_WITH([cplex-include],
    [AS_HELP_STRING([--with-cplex-include],[specify location of CPLEX C++ include directory])],
    CPLEX_INCLUDE=$withval,
    CPLEX_INCLUDE="")

AC_ARG_WITH([cplex-lib],
    [AS_HELP_STRING([--with-cplex-lib],[specify location of CPLEX callable library and other libraries])],
    CPLEX_LIB=$withval,
    CPLEX_LIB="")

AC_ARG_WITH([cplex-dir],
    [AS_HELP_STRING([--with-cplex-dir],[specify location of CPLEX installation directory])],
    CPLEX_DIR=$withval,
    CPLEX_DIR="/opt/ibm/ILOG/CPLEX_Studio2211/cplex")

if test ["$FOUND_ALL" = "no" -a -n "$CPLEX_DIR"] ; then
   PKG_CPPFLAGS="-I${CPLEX_DIR}/include -I${CPLEX_DIR}/include/ilcplex -I${CPLEX_DIR}/../concert/include"
   SYSTEM=`uname -m`
   LIBFORMAT=`uname -s | tr '[:upper:]' '[:lower:]'`
   PKG_LIBS="-L${CPLEX_DIR}/lib/x86-64_linux/static_pic -L${CPLEX_DIR}/../concert/lib/x86-64_linux/static_pic -lilocplex -lcplex -lconcert -lm -lpthread -ldl"

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

