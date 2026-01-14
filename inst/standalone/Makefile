# Standalone Makefile

SRC := ../../src

# If MODEL is not provided, set blank as default
MODEL ?= $(SRC)

# Standard fields used for debugging
# exported by .__enclos_env__$private$.export_to_cpp_standalone
STD_IN := -A augm.csv -D data_wide.csv -L data_long.csv -G group_rates.csv -I infec.csv -P param.csv -S setup.csv -M multiplicity.csv

# Lost of source files with exclusions
EXCL_FILES := $(SRC)/export_bindings.cpp $(SRC)/export_wrappers.cpp $(SRC)/RcppExports.cpp $(SRC)/model.cpp
SOURCES := main.cpp $(MODEL)/model.cpp $(filter-out $(EXCL_FILES), $(wildcard $(SRC)/*.cpp))

# Object files to build (automatic)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(MODEL)/model.exe

# compilation flags
CXX := g++
CXXFLAGS := -std=c++17 -g -Wall -Wextra -pedantic -Wno-unknown-pragmas -I../template/include/ 
ifdef DEBUG
CXXFLAGS += -O0 -Werror -Wuninitialized
else
# standard Rcpp compilation flags
CXXFLAGS += -O2
# Change openMP flags according to OS
ifeq ($(shell uname -s),Linux)
CXXFLAGS += -fopenmp
else ifeq ($(shell uname -s),Darwin)
ifeq ($(shell command -v brew),)
$(info Homebrew is not installed. Falling back to default.)
else
ifeq ($(shell brew list --versions llvm),)
$(info LLVM is not installed via Homebrew. Falling back to default.)
else
CXX := $(shell brew --prefix llvm)/bin/clang++
CXXFLAGS += -fopenmp -I$(shell brew --prefix llvm)/include -I$(shell brew --prefix boost)/include
endif
endif
else
$(error Unsupported operating system.)
endif
endif

# Make rules
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

# For the model object
$(MODEL)/model.o: $(MODEL)/model.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# For the main object
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I$(MODEL) -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
clean:
	$(RM) main.o {$(MODEL),$(SRC),../models/*}/*.{o,so,exe}
	$(RM) -r ../template/src/

# Run make in debugging mode using standard housepi files
debug:
	$(MAKE) DEBUG=1

# Test log-likelihood evaluation 
test-log-lik:
	$(EXECUTABLE) $(STD_IN) --log-lik

# Test inference 
test-infer:
	$(EXECUTABLE) $(STD_IN) --infer -O out.csv

# Test simulation 
test-simulate:
	$(EXECUTABLE) $(STD_IN) --simulate -O out.csv

.PHONY: all clean debug test-log-lik test-infer test-simulate
