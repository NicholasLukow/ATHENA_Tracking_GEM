AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libGEM_g4alsupp.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(G4_MAIN)/include \
    -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libGEM_g4alsupp_la_LIBADD = \
  -lSubsysReco \
  -lg4detectors \
  -lg4testbench 

pkginclude_HEADERS = \
  GEM_support_Subsystem.h

libGEM_g4alsupp_la_SOURCES = \
  GEM_support_Detector.cc \
  GEM_support_Subsystem.cc \
  GEM_support_SteppingAction.cc

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4alsupp


testexternals_g4alsupp_SOURCES = testexternals.cc
testexternals_g4alsupp_LDADD = libGEM_g4alsupp.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
