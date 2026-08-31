# Standalone Makefile

# If MODEL is not provided, set 'test' model as default
MODEL ?= ../../src/model_test.cpp

# Standard fields used for debugging
# exported by private$.export_to_cpp_standalone
STD_IN := -A latent.csv \
-C contact_risks.csv \
-D data_wide.csv \
-L data_long.csv \
-G group_risks.csv \
-I infec.csv \
-P param.csv \
-S setup.csv

# Lost of source files with exclusions
SRC := ../../src
EXCL_FILES := $(SRC)/RcppExports.cpp $(wildcard $(SRC)/model_*.cpp)
SOURCES := standalone.cpp $(MODEL) $(filter-out $(EXCL_FILES), $(wildcard $(SRC)/*.cpp))

# Object files to build (automatic)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
DIRMODEL := $(dir $(MODEL))
EXECUTABLE := $(DIRMODEL)/model.exe

# Boost headers come from the BH package in the R build; locate them for standalone
BH_INC := $(shell Rscript -e 'cat(system.file("include", package="BH"))')

# compilation flags
CXX := g++
CXXFLAGS := -std=c++17 -D_STANDALONE\
-I./ -I../include/ -I$(BH_INC) \
-Wall -Wextra -pedantic \
-Wno-unknown-pragmas
LDFLAGS :=

ifneq ($(filter debug,$(MAKECMDGOALS)),)
  DEBUG := 1
endif
ifneq ($(filter profile,$(MAKECMDGOALS)),)
  PROFILE := 1
endif

ifdef DEBUG
CXXFLAGS += -g3 -O0 \
-Werror -Wuninitialized -Wunused \
-fsanitize=address,undefined \
-fno-sanitize-recover=all \
-fno-omit-frame-pointer \
-fstack-protector-strong \
-D_GLIBCXX_DEBUG \
-D_GLIBCXX_DEBUG_PEDANTIC \
-fno-inline
else ifdef PROFILE
CXXFLAGS += -g -O3 -DNDEBUG\
-fno-omit-frame-pointer
else
CXXFLAGS += -g -O2 -DNDEBUG
endif

ifdef DISCRETE
CXXFLAGS += -D_HOUSEPI_DISCRETE_TIME
endif

# Standalone builds run single-threaded: parallelism is delegated to
# RcppParallel which only exists inside the R package. Pick up Boost
# headers from Homebrew on macOS when available.
ifeq ($(shell uname -s),Linux)
# nothing to add
else ifeq ($(shell uname -s),Darwin)
ifneq ($(shell command -v brew),)
ifneq ($(shell brew list --versions boost),)
CXXFLAGS += -I$(shell brew --prefix boost)/include
endif
endif
else
$(error Unsupported operating system.)
endif

# Make rules
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $(EXECUTABLE)

# For the model object
$(DIRMODEL)/model.o: $(MODEL)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# For the standalone object
standalone.o: standalone.cpp
	$(CXX) $(CXXFLAGS) -I$(DIRMODEL) -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
clean:
	$(RM) standalone.o
	$(RM) ../models/*/model.cpp
	$(RM) {$(DIRMODEL),$(SRC),../models/*}/*.{o,tmp,exe,so,dll}

# Run make in debugging mode
debug: all

# Run make in profiling mode
profile: all

# Test log-likelihood evaluation
test-log-lik: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --log-lik

# Test inference
test-infer: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --infer -O out.csv

# Test simulation
test-simulate: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --simulate -O out.csv

.PHONY: all clean debug profile test-log-lik test-infer test-simulate
