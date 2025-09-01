# Standalone Makefile

SRC := ../template/src
TMPL := ../models

# If MODEL is not provided, set blank as default
MODEL ?= ${SRC}

# Lost of source files with exclusions
EXCL_FILES := ${SRC}/export_bindings.cpp ${SRC}/RcppExports.cpp ${SRC}/model.cpp
SOURCES := main.cpp $(MODEL)/model.cpp $(filter-out $(EXCL_FILES), $(wildcard ${SRC}/*.cpp))

# Object files to build (automatic)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(MODEL)/model.exe

# Compilation flags
CXXFLAGS := -Wall -Wextra -std=gnu++14 -I/../include/
ifdef FAST
CXXFLAGS += -Ofast
else
CXXFLAGS += -g -O2 -Wno-unknown-pragmas

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
	$(RM) $(EXECUTABLE) main.o ${SRC}/*.{o,so,exe} ${TMPL}/*/*.{o,so,exe} ${TMPL}/*/.refer 

# To run make in a debug mode
.PHONY: fast
fast:
	$(MAKE) FAST=1