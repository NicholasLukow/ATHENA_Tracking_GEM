AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4_modifieddrich.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4_modifieddrich_la_LIBADD = \
  -lfun4all \
  -lg4detectors 

pkginclude_HEADERS = \
  dRICH_Subsystem.h 

libg4_modifieddrich_la_SOURCES = \
  dRICH_Detector.cc \
  dRICH_SteppingAction.cc \
  dRICH_Subsystem.cc 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4_modifieddrich.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
