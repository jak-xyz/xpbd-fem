
SIM_H = $(filter-out %vectormathx8.h %ThreadPool.h,$(wildcard ../XPBDFEM/*.h))
SIM_CPP = $(filter-out %main.cpp %GlRenderer.cpp %ThreadPool.cpp,$(wildcard ../XPBDFEM/*.cpp))
SIM_SRC = $(SIM_H) $(SIM_CPP)
WASM_SRC = $(SIM_SRC:../XPBDFEM/%=build-artifacts/%)

BL_SRC = main.cpp printf.h printf.cpp wasm-compat.h wasm-compat.cpp
BR_SRC = $(BL_SRC:%=build-artifacts/%)

build : sim.wasm
.PHONY : build

build-artifacts/force:
	mkdir build-artifacts
	echo . > build-artifacts/force

$(WASM_SRC) : build-artifacts/force
$(WASM_SRC) : build-artifacts% : ../XPBDFEM%
	sed 's/^\s*#\s*include\s*<.*>/#include \"wasm-compat.h\"/' $< > $@

$(BR_SRC) : build-artifacts/force
$(BR_SRC) : build-artifacts/% : %
	cp $< $@

HDR = $(filter %.h,$(WASM_SRC) $(BR_SRC))
CPP = $(filter %.cpp,$(WASM_SRC) $(BR_SRC))
OBJ = $(CPP:.cpp=.o)

$(OBJ) : $(HDR) # Conservatively rebuild everything whenever a header changes
$(OBJ) : Makefile # Similarly for when the makefile changes

$(OBJ) : %.o : %.cpp
	clang-12 \
		--target=wasm32 -c -x c++ -std=c++17 -nostdlib -nostdinc -O3 \
		-fno-math-errno -fno-rtti -fno-exceptions -fno-vectorize -fno-slp-vectorize \
		-o $@ $<

sim.wasm : $(OBJ)
	wasm-ld-12 \
		--no-entry --export-dynamic --import-memory --initial-memory=50331648 -z stack-size=1048576 \
		--allow-undefined --export=__heap_base \
		-o $@ $^
