# Standalone Makefile

# If MODEL is not provided, set 'test' model as default
MODEL ?= ../../src/model_test.cpp

# Standard fields used for debugging
# exported by private$.export_to_cpp_standalone
STD_IN := -A latent.csv \
-C contact_rates.csv \
-D data_wide.csv \
-L data_long.csv \
-G group_rates.csv \
-I infec.csv \
-P param.csv \
-S settings.csv

# Lost of source files with exclusions
SRC := ../../src
EXCL_FILES := $(SRC)/RcppExports.cpp $(wildcard $(SRC)/model_*.cpp)
SOURCES := standalone.cpp $(MODEL) $(filter-out $(EXCL_FILES), $(wildcard $(SRC)/*.cpp))
DIRMODEL := $(dir $(MODEL))

# Build flavour, picked up by the flags below
ifneq ($(filter debug,$(MAKECMDGOALS)),)
  DEBUG := 1
endif
ifneq ($(filter profile,$(MAKECMDGOALS)),)
  PROFILE := 1
endif
ifdef DEBUG
  FLAVOUR := debug
else ifdef PROFILE
  FLAVOUR := profile
else
  FLAVOUR := release
endif
ifdef DISCRETE
  FLAVOUR := $(FLAVOUR)-discrete
endif

# Objects and binary live under build/<flavour>, never in $(SRC): those object
# files belong to R CMD INSTALL. Sharing them makes the two builds trade object
# files -- the standalone link fails on R-only symbols (unif_rand, _main), and
# R's install silently pulls in objects compiled with -D_STANDALONE. One
# directory per flag set, so flavours never reuse each other's objects either.
BUILD := build/$(FLAVOUR)
OBJECTS := $(addprefix $(BUILD)/,$(notdir $(SOURCES:.cpp=.o)))
DEPENDS := $(OBJECTS:.o=.d)
EXECUTABLE := $(BUILD)/model.exe
MODEL_STAMP := $(BUILD)/model.stamp

# sources are gathered from several directories, objects land in one
vpath %.cpp . $(SRC) $(DIRMODEL)

# Boost headers come from the BH package in the R build; locate them for standalone
BH_INC := $(shell Rscript -e 'cat(system.file("include", package="BH"))')
# TBB headers and library come from the RcppParallel package, as in the R build
TBB_DIR := $(shell Rscript -e 'cat(system.file(package="RcppParallel"))')

# compilation flags
CXX := g++
CXXFLAGS := -std=c++17 -D_STANDALONE\
-I./ -I../include/ -I$(BH_INC) -I$(TBB_DIR)/include \
-Wall -Wextra -pedantic \
-Wno-unknown-pragmas \
-MMD -MP
LDFLAGS := -L$(TBB_DIR)/lib -Wl,-rpath,$(TBB_DIR)/lib -ltbb

ifdef DEBUG
CXXFLAGS += -g3 -O0 \
-Werror -Wuninitialized -Wunused \
-fsanitize=address,undefined \
-fno-sanitize-recover=all \
-fno-omit-frame-pointer \
-fstack-protector-strong \
-D_GLIBCXX_DEBUG \
-D_GLIBCXX_DEBUG_PEDANTIC \
-fno-inline
else ifdef PROFILE
CXXFLAGS += -g -O3 -DNDEBUG\
-fno-omit-frame-pointer
else
CXXFLAGS += -g -O2 -DNDEBUG
endif

ifdef DISCRETE
CXXFLAGS += -D_HOUSEPI_DISCRETE_TIME
endif

# Pick up Boost headers from Homebrew on macOS when available.
ifeq ($(shell uname -s),Linux)
# nothing to add
else ifeq ($(shell uname -s),Darwin)
ifneq ($(shell command -v brew),)
ifneq ($(shell brew list --versions boost),)
CXXFLAGS += -I$(shell brew --prefix boost)/include
endif
endif
else
$(error Unsupported operating system.)
endif

# Make rules
all: $(EXECUTABLE)

# For the executable, relinked when MODEL changes: object timestamps alone
# keep a binary linked against the previous model
$(EXECUTABLE): $(OBJECTS) $(MODEL_STAMP)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

# Stamp holding the model last linked, rewritten only when it differs
$(MODEL_STAMP): FORCE | $(BUILD)
	@echo '$(MODEL)' | cmp -s - $@ || echo '$(MODEL)' > $@

.PHONY: FORCE
FORCE:

# For the standalone object, the only one reading the model's header
$(BUILD)/standalone.o: standalone.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(DIRMODEL) -c $< -o $@

# For all other objects
$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD):
	@mkdir -p $@

# To clean objects
clean:
	$(RM) -r build
	$(RM) ../models/*/model.cpp

# Run make in debugging mode
debug: all

# Run make in profiling mode
profile: all

# Test log-likelihood evaluation
test-log-lik: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --log-lik

# Test inference
test-infer: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --infer -O out.csv

# Test simulation
test-simulate: $(EXECUTABLE)
	$(EXECUTABLE) $(STD_IN) --simulate -O out.csv

.PHONY: all clean debug profile test-log-lik test-infer test-simulate

# header dependencies recorded by -MMD, absent on the first build
-include $(DEPENDS)
