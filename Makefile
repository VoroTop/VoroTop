########################################################
####                                                ####
####   ******************************************   ####
####   *                                        *   ####
####   *     VoroTop: Voronoi Cell Topology     *   ####
####   *   Visualization and Analysis Toolkit   *   ####
####   *             (Version 1.0)              *   ####
####   *                                        *   ####
####   *           Emanuel A. Lazar             *   ####
####   *          Bar Ilan University           *   ####
####   *               June 2024                *   ####
####   *                                        *   ####
####   ******************************************   ####
####                                                ####
########################################################

### Compiler Configuration ###
CXX        := g++-mp-13
CXXSTD     := -std=c++11
CXXFLAGS   := -Wall -Wextra -O3 -fopenmp -MMD -MP
LDFLAGS    := -fopenmp
LDLIBS     := -lvoro++

### Project Structure ###
BUILD_DIR  := build
SRC_DIR    := .
TARGET     := VoroTop
PREFIX     := /usr/local# Default installation prefix
BIN_DIR    := $(PREFIX)/bin

# Source and object files
SOURCES    := $(wildcard $(SRC_DIR)/*.cc)
OBJECTS    := $(patsubst $(SRC_DIR)/%.cc,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS       := $(OBJECTS:.o=.d)

### Targets ###
.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(@D)
	$(CXX) $(CXXSTD) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

install: $(TARGET)
	@echo "Installing $(TARGET) to $(BIN_DIR)"
	@mkdir -p $(BIN_DIR)
	@install -m 755 $(TARGET) $(BIN_DIR)/$(TARGET)
	@echo "$(TARGET) installed successfully!"


