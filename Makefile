CXX := g++
CXXFLAGS := -Wall -Wextra -std=gnu++17 -O3

# List of source files
SOURCES := main.cpp mcmc.cpp model.cpp individual.cpp varia.cpp utils.cpp

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
