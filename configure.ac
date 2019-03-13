AC_PREREQ([2.69])
AC_INIT(HHG, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME, respecting an environment variable if one is set 
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
## Use R to set CXX and CXXFLAGS

CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
AC_LANG(C++)


##CXX=$(${R_HOME}/bin/R CMD config CXX)
##CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
##AC_LANG(C++)
AC_REQUIRE_CPP


#AC_INIT([V_CPPFLAGS],' -pthread')
#AC_INIT([V_LIBS],' -pthread')

case "uname" in
*SunOS*)
V_CPPFLAGS=" -D_REENTRANT -ICstd -ICrun"
V_LIBS=" -lpthread"

AC_SUBST(V_CPPFLAGS)
AC_SUBST(V_LIBS)

;;
*)
V_CPPFLAGS=" -pthread"
V_LIBS=" -pthread"

AC_SUBST(V_CPPFLAGS)
AC_SUBST(V_LIBS)

;;
esac


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
