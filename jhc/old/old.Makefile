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

CMAINx = pointcloud_2Dprofile.cpp

CCMP = g++

COPT = -g -Wall

LIBS = -llapack

#-------------------------------------

INC = $(addprefix -I,$(INCDIR))

SRC =   $(addprefix $(SRCDIR),$(CSRC))

MAIN = $(addprefix $(SCRIPTDIR),$(CMAIN))
#MAINx= $(addprefix $(SCRIPTDIR),$(CMAINx))

PROGS      = $(addprefix $(RUNDIR),$(CMAIN:.cpp=.out))
#PROGSNAMES = $(subst .cpp,,$(CMAIN))
PROGSNAMES = $(CMAIN:.cpp=)
#RUNx= $(addprefix $(RUNDIR),$(CMAINx:.cpp=.out))

SOBJ=$(SRC:.cpp=.o)

MOBJ = $(MAIN:.cpp=.o)
#MOBJx= $(MAINx:.cpp=.o)

#OBJS =$(SOBJ) $(MOBJ)
#OBJSx=$(SOBJ) $(MOBJx)

#-------------------------------------

#$(RUNx) : $(OBJSx) $(MOBJx)
#	$(CCMP) $(COPT) $(INC) -o $@ $(OBJSx) $(LIBS)

#$(RUNx) : $(SOBJ) $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o))
#	$(CCMP) $(COPT) $(INC) -o $@ $(SOBJ) $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o)) $(LIBS)

#$(RUNx) : $(SOBJ) $(subst $(RUNDIR),$(SCRIPTDIR),$(RUNx:.out=.o))
#	$(CCMP) $(COPT) $(INC) -o $@ $^ $(LIBS)

#.cpp.o:
#	$(CCMP) $(COPT) $(INC) -o $@ -c $<

%.o : %.cpp
	$(CCMP) $(COPT) $(INC) -o $@ -c $<

#-------------------------------------

#$(RUN) : $(OBJS)
#	$(CCMP) $(COPT) $(INC) -o $@ $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o)) $(SOBJ) $(LIBS)

#$(RUN) : $(MOBJ) $(SOBJ)
#	$(CCMP) $(COPT) $(INC) -o $@ $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o)) $(SOBJ) $(LIBS)

#$(RUN) : $(MOBJ) $(SOBJ)
#	$(CCMP) $(COPT) $(INC) -o $@ $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o)) $(SOBJ) $(LIBS)

#$(RUNDIR)%.out : $(SCRIPTDIR)%.o $(SOBJ)
#	$(CCMP) $(COPT) $(INC) -o $@ $(subst $(RUNDIR),$(SCRIPTDIR),$(@:.out=.o)) $(SOBJ) $(LIBS)

#$(RUNDIR)%.out : $(SCRIPTDIR)%.o $(SOBJ)
#	$(CCMP) $(COPT) $(INC) -o $@ $^ $(LIBS)

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
