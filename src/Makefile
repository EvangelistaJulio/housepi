# If MODEL is not provided, set a default value
MODEL ?= models/chain_binomial/user_model.cpp

# Trim the file extension
NAME := $(notdir $(basename $(MODEL)))
FLDR := $(shell dirname ${MODEL})
MODEL_HEADER := $(NAME).h

# List of source files as found in the software directory
EXCL_FILES := Rcpp_export_bindings.cpp # FIXME (temporary)
SOURCES := $(filter-out $(EXCL_FILES), $(wildcard *.cpp) $(FLDR)/$(NAME).cpp standalone/main.cpp)

# Object files to build (automatically generated from SRCS)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(FLDR)/$(NAME).exe

# Common compilation flags
CXXFLAGS := -Ofast -Wall -Wextra -std=gnu++11

# Change openMP options according to the operating system
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

# Make rules
.PHONY: all
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
	
# For the model object
$(FLDR)/$(NAME).o: $(FLDR)/$(NAME).cpp
	$(CXX) $(CXXFLAGS) -I$(shell pwd) -c $< -o $@

# For the main object
standalone/main.o: standalone/main.cpp
	$(CXX) $(CXXFLAGS) -I$(shell pwd) -I$(FLDR) -DMODEL_HEADER=\"$(MODEL_HEADER)\" -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
.PHONY: clean
clean:
	$(RM) $(OBJECTS)