# If MODEL is not provided, set a default value
MODEL ?= example/default_model.cpp

# Trim the file extension
NAME := $(notdir $(basename $(MODEL)))
FLDR := $(shell dirname ${MODEL})

# List of source files
SRCDIR := src
SOURCES := $(wildcard $(SRCDIR)/*.cpp) $(FLDR)/$(NAME).cpp

# Object files to build (automatically generated from SRCS)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := $(FLDR)/$(NAME).exe

# compiler and rules
CXX := g++
CXXFLAGS := -Ofast -Wall -Wextra -std=gnu++11 -I$(shell pwd)/${SRCDIR}


# Make rules
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJECTS)
	$(info Building $(EXECUTABLE)...)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

# For the model object
$(NAME).o: $(FLDR)/$(NAME).cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
clean:
	rm -rf src/*.o
	rm -rf ${FLDR}/*.o
