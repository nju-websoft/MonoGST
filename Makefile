CXX = g++
CXXFLAGS = -std=c++20 -O2 -g

# MonoGST paths
MONOGST_INC = src/MonoGST/include
MONOGST_SRC = src/MonoGST/src

.PHONY: all clean

all: MonoGST DST ImprovAPP PartialOPT PrunedDP

# MonoGST tests
MonoGST: src/MonoGST/tests/test_Improved2star.cpp $(MONOGST_SRC)/GlobalUtils.cpp $(MONOGST_SRC)/Log.cpp
	$(CXX) $(CXXFLAGS) -I$(MONOGST_INC) -DDEBUG -o $@ src/MonoGST/tests/test_Improved2star.cpp $(MONOGST_SRC)/GlobalUtils.cpp $(MONOGST_SRC)/Log.cpp

# DST

DST: 
	mkdir -p src/DST/build
	cmake -S src/DST -B src/DST/build
	cmake --build src/DST/build

# Other algorithms (single-file builds)
ImprovAPP: src/ImprovAPP/ImprovAPP.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/ImprovAPP/ImprovAPP.cpp

PartialOPT: src/PartialOPT/PartialOPT.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/PartialOPT/PartialOPT.cpp

PrunedDP: src/PrunedDP/PrunedDP.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/PrunedDP/PrunedDP.cpp

KKG_Index: src/kkg/test/create_index.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/kkg/test/create_index.cpp

KKG_Run: src/kkg/test/run_kkg.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/kkg/test/run_kkg.cpp

2starh: src/2starh/2starh.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/2starh/2starh.cpp

clean:
	rm -f MonoGST MonoGST_ab ImprovAPP PartialOPT PrunedDP KKG_Index KKG_Run 2starh *.o src/*.o