bin_PROGRAMS = cbuccaneer csequins

cbuccaneer_SOURCES = \
buccaneer-build.cpp    buccaneer-join.h        buccaneer-prune.cpp \
buccaneer-build.h      buccaneer-lib.cpp       buccaneer-prune.h \
buccaneer-correct.cpp  buccaneer-lib.h         buccaneer-sequence.cpp \
buccaneer-correct.h    buccaneer-link.cpp      buccaneer-sequence.h \
buccaneer-filter.cpp   buccaneer-link.h        buccaneer-util.cpp \
buccaneer-filter.h     buccaneer-ncsbuild.cpp  buccaneer-util.h \
buccaneer-find.cpp     buccaneer-ncsbuild.h    cbuccaneer.cpp \
buccaneer-find.h       buccaneer-prep.cpp      simplex-lib.cpp \
buccaneer-grow.cpp     buccaneer-prep.h        simplex-lib.h \
buccaneer-grow.h       buccaneer-prot.cpp      simulate-lib.cpp \
buccaneer-join.cpp     buccaneer-prot.h        simulate-lib.h

csequins_SOURCES = \
csequins.cpp buccaneer-sequence.cpp simulate-lib.cpp buccaneer-prot.cpp \
buccaneer-lib.cpp


INCLUDES = $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS) \
           $(CCP4_CXXFLAGS) -I$(includedir)

LIBS = $(CLIPPER_LIBS) $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) @LIBS@


