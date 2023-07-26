CXX := g++
CXXFLAGS := -Wall -Wextra -std=gnu++14

# List of source files
SOURCES := main.cpp parameter.cpp household.cpp mcmc_object.cpp user_model.cpp

# Object files to build (automatically generated from SOURCES)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := program.exe

# Make rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.out
