bin_PROGRAMS = cbuccaneer

cbuccaneer_SOURCES = \
buccaneer-build.cpp     buccaneer-lib.cpp       cbuccaneer.cpp \
buccaneer-find.cpp      buccaneer-prot.cpp      simplex-lib.cpp \
buccaneer-grow.cpp      buccaneer-prune.cpp     simulate-lib.cpp \
buccaneer-join.cpp      buccaneer-sequence.cpp \
buccaneer-build.h       buccaneer-lib.h         simplex-lib.h \
buccaneer-find.h        buccaneer-prot.h        simulate-lib.h \
buccaneer-grow.h        buccaneer-prune.h \
buccaneer-join.h        buccaneer-sequence.h


INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) \
	   -I$(includedir)

LIBS = $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) @LIBS@

LDADD = ../../clipper/contrib/libclipper-contrib.la ../../clipper/cif/libclipper-cif.la ../../clipper/mmdb/libclipper-mmdb.la ../../clipper/ccp4/libclipper-ccp4.la ../../clipper/phs/libclipper-phs.la ../../clipper/core/libclipper-core.la ../../clipper/minimol/libclipper-minimol.la

