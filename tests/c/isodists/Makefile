#############################################################################
# Makefile for building: isocalc
#############################################################################

####### Compiler, tools and options IGNORE THESE
CC       = gcc
CXX      = g++
LEX      = flex
YACC     = yacc
CFLAGS   = -ansi -pipe -Wall -W -g -c
CXXFLAGS = -ansi -pipe -Wall -W -g -c
LEXFLAGS = 
YACCFLAGS= -d
LINK     = g++
LFLAGS   = 
AR       = ar rv
RANLIB   = ranlib
TAR      = tar cfz
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p
  
  
####### Source directories (there are two: this one and the one in the R package
VPATH = ./:../../../q2e/src
  
  
  
####### Output directory

OBJECTS_DIR = ./

####### Files
HEADERS = files.h \
SOURCES = isodists.c \
	  	isowrap.c \
	  	

TESTSOURCES = isodists.c \
	  	isowrap.c \
	  	testiso.c	  	
		
OBJECTS = isodists.o \
		isowrap.o \
		
TESTOBJECTS = isodists.o \
		isowrap.o \
		testiso.o

DIST = 	Makefile \

#IMB: IGNORE THESE!
FORMS = 
DESTDIR  = 
TARGET   = isocalc
TESTTARGET = testiso

first: all
####### Implicit rules

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

####### Build rules - IGNORE THESE

all: Makefile $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)
	-$(DEL_FILE) $(OBJECTS)

dist: 
	@$(TAR) image.tgz $(SOURCES) $(HEADERS) $(LIBSOURCES) $(LIBHEADERS) $(DIST);
#	@mkdir -p .tmp/image && $(COPY_FILE) --parents Makefile $(SOURCES) $(HEADERS) $(LIBSOURCES) $(LIBHEADERS) $(FORMS) $(DIST) .tmp/image/ && ( cd `dirname .tmp/image` && $(TAR) image.tgz image ) && $(MOVE) `dirname .tmp/image`/image.tgz . && $(DEL_FILE) -r .tmp/image

yaccclean:
lexclean:
clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) $(LIBOBJECTS)
	-$(DEL_FILE) *~ core *.core core*

####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) $(TARGET)


FORCE:

####### Compile

isodists.o: isodists.c constants.h files.h

isowrap.o: isowrap.c constants.h


####### Install

install:  

uninstall:  

