.SUFFIXES:
.SUFFIXES: .cpp .o

INCDIR = ./inc/
SRCDIR = ./src/
RUNDIR = ./run/
SCRIPTDIR = ./scripts/

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
        reconstruct_2Dprofile.cpp \
        pointcloud_2Dprofile.cpp \
        reconstruct_3Dprofile.cpp \
        pointcloud_3Dprofile.cpp

CMAINx = pointcloud_3Dprofile.cpp

CCMP = g++

COPT = -g -Wall

LIBS = -llapack

#-------------------------------------

INC = $(addprefix -I,$(INCDIR))

SRC =   $(addprefix $(SRCDIR),$(CSRC))

MAIN = $(addprefix $(SCRIPTDIR),$(CMAIN))
MAINx= $(addprefix $(SCRIPTDIR),$(CMAINx))

RUN = $(addprefix $(RUNDIR),$(CMAIN:.cpp=.out))
RUNx= $(addprefix $(RUNDIR),$(CMAINx:.cpp=.out))

SOBJ=$(SRC:.cpp=.o)

MOBJx= $(MAINx:.cpp=.o)
MOBJ = $(MAIN:.cpp=.o)

OBJSx=$(SOBJ) $(MOBJx)
OBJS =$(SOBJ) $(MOBJ)

#-------------------------------------

$(MAINx) : $(OBJSx)
	$(CCMP) $(COPT) $(INC) -o $@ $(OBJSx) $(LIBS)

.cpp.o:
	$(CCMP) $(COPT) $(INC) -o $@ -c $<

#-------------------------------------

all: $(RUN)

$(RUN) : $(OBJS)
	$(CCMP) $(COPT) $(INT) -o $@ -c $<

clean:
	rm -f $(OBJS) $(EXE)
