CXX      := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra

# ROOT flags
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# -----------------------------
# Directories
# -----------------------------
SRC_DIR := src
OBJ_DIR := obj
BIN     := build_pt2

# -----------------------------
# Sources & objects
# -----------------------------
SOURCES := \
	$(SRC_DIR)/main.cpp \
	$(SRC_DIR)/tools.cpp \
	$(SRC_DIR)/histograms.cpp \
	$(SRC_DIR)/plot_recipes.cpp \
	$(SRC_DIR)/plotting.cpp \
	$(SRC_DIR)/rootReader.cpp \
	$(SRC_DIR)/gator.cpp\
        $(SRC_DIR)/extrapolation.cpp

OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# -----------------------------
# Default target
# -----------------------------
all: $(BIN)

# -----------------------------
# Link
# -----------------------------
$(BIN): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(ROOTLIBS)

# -----------------------------
# Compile
# -----------------------------
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -I$(SRC_DIR) -c $< -o $@

# -----------------------------
# Create obj directory
# -----------------------------
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# -----------------------------
# Clean
# -----------------------------
clean:
	rm -rf $(OBJ_DIR) $(BIN)

.PHONY: all clean
