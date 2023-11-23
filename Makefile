CXX := g++
CXXFLAGS := -Wall -Wextra -std=gnu++11 -Ofast

# List of source files
SOURCES := main.cpp mcmc.cpp ${MODEL}.cpp  individual.cpp varia.cpp utils.cpp

# Object files to build (automatically generated from SOURCES)
OBJECTS := $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE := run-$(basename ${MODEL}).exe

# Make rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.out
