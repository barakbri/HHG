AC_INIT([V_CPPFLAGS]," -pthread")
AC_INIT([V_LIBS]," -pthread")

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


AC_CONFIG_FILES([src/Makevars.in])
AC_OUTPUT

