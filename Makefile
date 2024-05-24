# If MODEL is not provided, set a default value
MODEL ?= src/models/chain_binomial/chain_binomial.cpp

# Trim the file extension
NAME := $(notdir $(basename $(MODEL)))
FLDR := $(shell dirname ${MODEL})
MODEL_HEADER := $(NAME).h

# List of source files as found in the software directory
SRCDIR := src
EXCL_FILES := Rcpp_export_bindings.cpp # FIXME (temporary)
SOURCES := $(filter-out $(addprefix $(SRCDIR)/,$(EXCL_FILES)), $(wildcard $(SRCDIR)/*.cpp) $(FLDR)/$(NAME).cpp)

# Object files to build (automatically generated from SRCS)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(FLDR)/$(NAME).exe

# Common compilation flags
CXXFLAGS := -Ofast -Wall -Wextra -std=gnu++11 -fopenmp -I$(shell pwd)/$(SRCDIR)

# Change options according to the operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CXX := g++

else ifeq ($(UNAME_S),Darwin)
	ifeq ($(shell command -v brew),)
		$(error Homebrew is not installed. Please install Homebrew.)
	endif		
	
	ifeq ($(shell brew list --versions llvm),)
		$(error LLVM is not installed via Homebrew. Please install LLVM.)
	endif
	
	CXX := $(shell brew --prefix)/opt/llvm/bin/clang++
	CXXFLAGS += -I$(shell brew --prefix)/opt/llvm/include
else
	$(error Unsupported operating system.)
endif

# Make rules
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(info Building $(EXECUTABLE)...)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
	
# For the model object
$(NAME).o: $(FLDR)/$(NAME).cpp
	$(CXX) $(CXXFLAGS) -I${FLDR} -c $< -o $@

# For the main object
$(SRCDIR)/main.o: $(SRCDIR)/main.cpp
	$(CXX) $(CXXFLAGS) -I$(FLDR) -DMODEL_HEADER=\"$(MODEL_HEADER)\" -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
clean:
	rm -rf ./*.o
	rm -rf src/*.o
	rm -rf ${FLDR}/*.o
