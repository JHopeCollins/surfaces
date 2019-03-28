INCDIR = inc/
SRCDIR = src/
RUNDIR = run/
SCRIPTDIR = scripts/

CSRC = vec.cpp \
		 pseudosurf.cpp \
		 rbf_interp.cpp \
		 rbf_surf.cpp \
		 surf2.cpp \
		 surf.cpp \
		 ssurf.cpp

CMAIN = test_rbf.cpp \
        test_rbf_grad.cpp \
        test_2Drbf.cpp \
        test_3Drbf.cpp \
        test_sharpbranch.cpp \
        test_filletbranch.cpp \
		  pointcloud_2Dprofile.cpp \
		  pointcloud_3Dprofile.cpp \
        reconstruct_2Dprofile.cpp \
        reconstruct_3Dprofile.cpp

CCMP = g++

COPT = -g -Wall

LIBS = -llapack

#-------------------------------------

INC = $(addprefix -I,$(INCDIR))

SRC =   $(addprefix $(SRCDIR),$(CSRC))

MAIN = $(addprefix $(SCRIPTDIR),$(CMAIN))

PROGS      = $(addprefix $(RUNDIR),$(CMAIN:.cpp=.out))
PROGSNAMES = $(CMAIN:.cpp=)

SOBJ=$(SRC:.cpp=.o)

MOBJ = $(MAIN:.cpp=.o)

#-------------------------------------

%.o : %.cpp
	$(CCMP) $(COPT) $(INC) -o $@ -c $<

#-------------------------------------

$(PROGS) : $(RUNDIR)%.out : $(SCRIPTDIR)%.o $(SOBJ)
	$(CCMP) $(COPT) $(INC) -o $@ $^ $(LIBS)

#-------------------------------------

.PHONY : all clean echo obj $(PROGSNAMES)

$(PROGSNAMES) : % : $(RUNDIR)%.out

all: $(PROGS)

obj: $(SOBJ)

echo:
	echo

clean:
	rm -f $(SOBJ) $(MOBJ) $(PROGS)
