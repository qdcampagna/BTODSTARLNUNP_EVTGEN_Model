TOP=/home/lopamudra/HepPackages/evtgen/EvtGen_Install
LIB=$(TOP)/lib
INC=$(TOP)/include

all: run_dstarlnu

run_dstarlnu: dstarlnu_test.cc
	g++ -I $(INC) -o $@ $? -L $(LIB) -lEvtGen -Wl,-rpath,$(LIB) `root-config --cflags --libs`
