#  Default macros for compiling C and Fortran source code.

CC = gcc
CFLAGS = -fpic -O
FC = g77
FFLAGS = -fpic -fno-second-underscore -O
#FFLAGS = -fno-second-underscore -O
LD = ld
LDFLAG1 = -lm  -shared  -lf2c # (for g77)
#LDFLAG2 =  -lF90 -lCEPCF90 -L/opt/intel/compiler70/ia32/lib

F_ROUTINES = aarray.f agauss.f akrank.f arisk.f bhk.f bin.f \
   bj.f buckly.f cens.f chol.f coeff.f coxreg.f data1.f data2.f \
   datain.f datreg.f em.f emalgo.f factor.f gamma.f gehan.f grdprb.f \
   kmadj.f kmdif.f kmestm.f kmprnt.f lrank.f matinv.f mulvar.f pchisq.f \
   peto.f plestm.f pwlcxn.f r3.f r4.f r5.f r6.f ran1.f rearrn.f regres.f \
   rmills.f schmit.f sort1.f sort2.f spearho.f sprman.f stat.f sumry.f \
   syminv.f tie.f twokm.f twost.f unpack.f wlcxn.f xdata.f xvar.f 

F_UNIVAR_LIB_ROUTINES = twost_l.f
C_UNIVAR_LIB_ROUTINES = lunivar.c
F_UNIVAR_LIB_TEST_ROUTINES = lib_univar.f

F_BIVAR_LIB_ROUTINES = coxreg_l.f bhk_l.f sprman_l.f em_l.f bj_l.f mulvar_l.f
C_BIVAR_LIB_ROUTINES = lbivar.c
F_BIVAR_LIB_TEST_ROUTINES = lib_bivar.f

#  List of platform specific Fortran routines.  The source tar file will
#  contain a version of these files for each set of platforms.  The names
#  given here are the general version of the file, e.g. pkg_open.f.  The
#  names in the tar file will be pkg_open.f_sun4, pkg_open.f_mips, etc.

PLATFORM_F = asurv_main.f bivar.f univar.f

#  List of executables.

EXECUTABLE = asurv

UNIVAR_LIBRARY = univar.so
BIVAR_LIBRARY  = bivar.so

################################################################################
#-------------------------------------------------------------------------------

#  Define files required for building the package.
#  ==============================================
#
#  This section defines the set of files produced from the source files
#  when the package is built and installed.

#  Use only .o, .c and .f suffix rules.

.SUFFIXES:
.SUFFIXES: .o .c .f

#  List of object files produced by compiling the source code and rules
#  for performing the compilations.

OBJECT_FILES = $(F_ROUTINES:.f=.o) $(PLATFORM_F:.f=.o) 


# Univar lib
UNIVAR_LIB_OBJECT_FILES = $(F_ROUTINES:.f=.o) \
	 $(F_UNIVAR_LIB_ROUTINES:.f=.o) $(C_UNIVAR_LIB_ROUTINES:.c=.o)

UNIVAR_LIB_TEST_OBJECT_FILES = $(F_ROUTINES:.f=.o) \
	 $(F_UNIVAR_LIB_ROUTINES:.f=.o) $(C_UNIVAR_LIB_ROUTINES:.c=.o) \
	 $(F_UNIVAR_LIB_TEST_ROUTINES:.f=.o)


# Bivar lib
BIVAR_LIB_OBJECT_FILES = $(F_ROUTINES:.f=.o) \
	 $(F_BIVAR_LIB_ROUTINES:.f=.o) $(C_BIVAR_LIB_ROUTINES:.c=.o)

BIVAR_LIB_TEST_OBJECT_FILES = $(F_ROUTINES:.f=.o) \
	 $(F_BIVAR_LIB_ROUTINES:.f=.o) \
	 $(F_BIVAR_LIB_TEST_ROUTINES:.f=.o)



.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<

#  Rules for building executables.

all: $(EXECUTABLE) $(UNIVAR_LIBRARY) $(BIVAR_LIBRARY) bivar_lib_test univar_lib_test


$(EXECUTABLE): $(OBJECT_FILES)
	$(FC) -o $@ $? # $(LDFLAG1) $(LDFLAG2)



$(UNIVAR_LIBRARY): $(UNIVAR_LIB_OBJECT_FILES)
	$(LD) -o $@ $? $(LDFLAG1) $(LDFLAG2)

univar_lib_test: $(UNIVAR_LIB_TEST_OBJECT_FILES)
	$(FC) -o $@ $? 


$(BIVAR_LIBRARY): $(BIVAR_LIB_OBJECT_FILES)
	$(LD) -o $@ $? $(LDFLAG1) $(LDFLAG2)

bivar_lib_test: $(BIVAR_LIB_TEST_OBJECT_FILES)
	$(FC) -o $@ $? 



clean: 
	rm -f *.o 

cleanall: clean
	rm -f $(EXECUTABLE) $(UNIVAR_LIBRARY) univar_lib_test bivar_lib_test \#*\# core*