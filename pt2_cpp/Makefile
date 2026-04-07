# ===============================
# Compiler
# ===============================
CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -MMD -MP

# ===============================
# ROOT
# ===============================
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# ===============================
# ONNX Runtime
# ===============================
ONNX_DIR := /blue/avery/matthew.dittrich/PT2_Studies/lst_pt2/pt2_cpp/onnxruntime-linux-x64-1.17.0
ONNX_INC := -I$(ONNX_DIR)/include
ONNX_LIB := -L$(ONNX_DIR)/lib -lonnxruntime

# ===============================
# Directories
# ===============================
SRC_DIR   := src
TOOLS_DIR := $(SRC_DIR)/tools
INC_DIR   := include
BUILD_DIR := build
BIN_DIR   := bin
BIN       := $(BIN_DIR)/pt2

# ===============================
# Sources & Objects
# ===============================
SOURCES := $(SRC_DIR)/main.cpp $(wildcard $(TOOLS_DIR)/*.cpp)
OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(SOURCES)))

# ===============================
# Targets
# ===============================
all: $(BIN)

# Link
$(BIN): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(ONNX_LIB) $^ -o $@

# Compile main src files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(ONNX_INC) -I$(INC_DIR) -c $< -o $@

# Compile tools src files
$(BUILD_DIR)/%.o: $(TOOLS_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(ONNX_INC) -I$(INC_DIR) -c $< -o $@

# ===============================
# Create directories
# ===============================
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# ===============================
# Clean
# ===============================
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# ===============================
# Dependencies
# ===============================
-include $(OBJECTS:.o=.d)

.PHONY: all clean
