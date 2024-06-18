all: Reweight

Reweight: reweightPPFX.cpp reweight_root_dict.o
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

reweight_root_dict.o:
	$(RM) reweight_root_dict*.*
	rootcling -f reweight_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o reweight_root_dict.o -c reweight_root_dict.cc
	$(RM) reweight_root_dict.cc

.PHONY: clean

.INTERMEDIATE: reweight_root_dict.o

clean:
	$(RM) Reweight reweight_root_dict.o reweight_root_dict_rdict.pcm