bin_PROGRAMS = cbuccaneer csequins

cbuccaneer_SOURCES = \
buccaneer-build.cpp     buccaneer-join.cpp      buccaneer-sequence.cpp \
buccaneer-correct.cpp   buccaneer-lib.cpp       cbuccaneer.cpp \
buccaneer-filter.cpp    buccaneer-link.cpp      simplex-lib.cpp \
buccaneer-find.cpp      buccaneer-prot.cpp      simulate-lib.cpp \
buccaneer-grow.cpp      buccaneer-prune.cpp 	buccaneer-ncsbuild.cpp \
buccaneer-build.h       buccaneer-join.h        buccaneer-sequence.h \
buccaneer-correct.h     buccaneer-lib.h         simplex-lib.h \
buccaneer-filter.h      buccaneer-link.h        simulate-lib.h \
buccaneer-find.h        buccaneer-prot.h \
buccaneer-grow.h        buccaneer-prune.h	buccaneer-ncsbuild.h

csequins_SOURCES = \
csequins.cpp buccaneer-sequence.cpp simulate-lib.cpp buccaneer-prot.cpp \
buccaneer-lib.cpp


INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) \
	   -I$(includedir)

LIBS = $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) @LIBS@

LDADD = ../../clipper/contrib/libclipper-contrib.la ../../clipper/cif/libclipper-cif.la ../../clipper/mmdb/libclipper-mmdb.la ../../clipper/ccp4/libclipper-ccp4.la ../../clipper/phs/libclipper-phs.la ../../clipper/core/libclipper-core.la ../../clipper/minimol/libclipper-minimol.la

