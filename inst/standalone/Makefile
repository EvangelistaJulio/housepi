# Standalone Makefile

# If MODEL is not provided, set a default value
MODEL ?= ../models/chain-binomial

# Lost of source files with exclusions
SRC := ../template/src
EXCL_FILES := ${SRC}/export_bindings.cpp ${SRC}/RcppExports.cpp ${SRC}/model.cpp
SOURCES := $(filter-out $(EXCL_FILES), $(wildcard ${SRC}/*.cpp) $(MODEL)/model.cpp main.cpp)

# Object files to build (automatic)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(MODEL)/model.exe

# Compilation flags
CXXFLAGS := -Wall -Wextra -std=gnu++14
ifdef DEBUG
CXXFLAGS += -g -O2 -Wno-unknown-pragmas
else
CXXFLAGS += -Ofast
# Change openMP flags according to OS
ifeq ($(shell uname -s),Linux)
CXX := g++
CXXFLAGS += -fopenmp
else ifeq ($(shell uname -s),Darwin)
ifeq ($(shell command -v brew),)
$(info Homebrew is not installed. Falling back to default.)
CXX := g++
else
ifeq ($(shell brew list --versions llvm),)
$(info LLVM is not installed via Homebrew. Falling back to default.)
CXX := g++
else
CXX := $(shell brew --prefix)/opt/llvm/bin/clang++
CXXFLAGS += -fopenmp -I$(shell brew --prefix)/opt/llvm/include
endif
endif	
else
$(error Unsupported operating system.)
endif
endif

# Make rules
.PHONY: all
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
	
# For the model object
$(MODEL)/model.o: $(MODEL)/model.cpp
	$(CXX) $(CXXFLAGS) -I${SRC} -c $< -o $@

# For the main object
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I${SRC} -I$(MODEL) -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(EXECUTABLE)

# To run make in a debug mode
.PHONY: debug
debug:
	$(MAKE) DEBUG=1