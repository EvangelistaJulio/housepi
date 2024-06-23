# Standalone Makefile

# temporay : added a DEBUG mode
DEBUG ?= 0

# If MODEL is not provided, set a default value
MODEL ?= ../models/chain_binomial/model.cpp

# Trim the file extension
NAME := $(notdir $(basename $(MODEL)))
FLDR := $(shell dirname ${MODEL})
MODEL_HEADER := $(NAME).h

# Lost of source files with exclusions
SRC := ../template/src
EXCL_FILES := ${SRC}/export_bindings.cpp ${SRC}/RcppExports.cpp ${SRC}/user_model.cpp
SOURCES := $(filter-out $(EXCL_FILES), $(wildcard ${SRC}/*.cpp) $(FLDR)/$(NAME).cpp main.cpp)

# Object files to build (automatic)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(FLDR)/$(NAME).exe

# Compilation flags
CXXFLAGS := -Wall -Wextra -std=gnu++11
ifeq ($(DEBUG), 1)
CXXFLAGS += -g -O2
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
$(FLDR)/$(NAME).o: $(FLDR)/$(NAME).cpp
	$(CXX) $(CXXFLAGS) -I${SRC} -c $< -o $@

# For the main object
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I${SRC} -I$(FLDR) -DMODEL_HEADER=\"$(MODEL_HEADER)\" -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
.PHONY: clean
clean:
	$(RM) $(OBJECTS)