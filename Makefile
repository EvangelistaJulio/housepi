# If MODEL is not provided, set a default value
MODEL ?= default_model

# Trim the file extension
NAME := $(notdir $(basename $(MODEL)))
FLDR := $(shell dirname ${MODEL})

CXX := g++
CXXFLAGS := -Ofast -Wall -Wextra -std=gnu++11 -I$(shell pwd)

# List of source files
SRCS := main.cpp mcmc.cpp $(FLDR)/$(NAME).cpp individual.cpp varia.cpp utils.cpp

# Object files to build (automatically generated from SRCS)
OBJS := $(patsubst %.cpp, %.o, $(notdir $(SRCS)))

# Executable name
EXECUTABLE := $(FLDR)/$(NAME).exe

# Make rules
all: $(EXECUTABLE)

# For the executable
$(EXECUTABLE): $(OBJS)
	$(info Building $(EXECUTABLE)...)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXECUTABLE)

# For the model object
$(NAME).o: $(FLDR)/$(NAME).cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# For all other objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# To clean objects
clean:
	rm -rf *.o
