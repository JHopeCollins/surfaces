#-------------------------------------
# USER INPUTS:

# directories
INCLDEDIR = inc/#			header files
SOURCEDIR = src/#			class / function source files
SCRIPTDIR = script/#		main() function source files
PROGRMDIR = progrm/#		executables

# Class / function definition source files
CSOURCE = pseudosurf.cpp \
			 surf2.cpp \
			 ssurf.cpp

# main() function files
CSCRIPT = test_rbf.cpp \
			 test_rbf_grad.cpp \
        	 test_2Drbf.cpp \
        	 test_3Drbf.cpp \
        	 test_sharpbranch.cpp \
        	 test_filletbranch.cpp \
        	 test_filletpass.cpp \
		  	 pointcloud_2Dprofile.cpp \
		  	 pointcloud_3Dprofile.cpp \
        	 reconstruct_2Dprofile.cpp \
        	 reconstruct_3Dprofile.cpp

CCMP = g++

COPT = -g -Wall

LIBS = -llapack


#-------------------------------------
#  variable definitions

DIRS = $(INCLDEDIR) $(SOURCEDIR) $(SCRIPTDIR) $(PROGRMDIR)

INCLDE = $(addprefix -I,$(INCLDEDIR))

# full paths for source, script and executable files
SOURCE = $(addprefix $(SOURCEDIR),$(CSOURCE))
SCRIPT = $(addprefix $(SCRIPTDIR),$(CSCRIPT))
PROGRM = $(addprefix $(PROGRMDIR),$(CSCRIPT:.cpp=.out))

# names of scripts (no suffix)
PNAMES = $(CSCRIPT:.cpp=)

# object files
SOURCEOBJ = $(SOURCE:.cpp=.o)
SCRIPTOBJ = $(SCRIPT:.cpp=.o)

OBJS = $(SOURCEOBJ) $(SCRIPTOBJ)


#-------------------------------------
# compilation recipes

# default
%.o : %.cpp
	$(CCMP) $(COPT) $(INCLDE) -o $@ -c $<

# wraps all source object files into one
all.o : $(SOURCEOBJ)
	ld -r $^ -o $@

# each executable depends on its own object file, and all source objects
$(PROGRM) : $(PROGRMDIR)%.out : $(SCRIPTDIR)%.o $(SOURCEOBJ)
	$(CCMP) $(COPT) $(INCLDE) -o $@ $^ $(LIBS)


#-------------------------------------
# misc recipes

.PHONY : all clean echo obj mkdir names $(PNAMES)

# make <pname> will compile only the executable 'pname.out'
$(PNAMES) : % : $(PROGRMDIR)%.out

# make all executables
all: $(PROGRM)

# make source objects
obj: $(SOURCEOBJ)

#-------------------------------------
# command recipes

# create required directories
mkdir:
	mkdir $(DIRS)

# print names of all executables to standard output
names:
	@for name in $(PNAMES); do echo $$name; done

# delete all non-source files
clean:
	rm -f $(OBJS) $(PROGRM)


#-------------------------------------
