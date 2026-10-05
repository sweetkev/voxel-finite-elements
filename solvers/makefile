-include ../Make.user

ifndef MFEM_DIR
$(error MFEM_DIR is not set)
endif

ifndef SRC_DIR
$(error SRC_DIR is not set)
endif

include $(MFEM_DIR)/config/config.mk

BUILD_DIR := build

MESH_SRCS := \
	meshing/ppm.cpp \
	meshing/pixel_mesh.cpp

GRAPH_SRCS := \
	graphs/voxel_graph.cpp \
	graphs/voxel_graph_operator.cpp \
	graphs/voxel_graph_hierarchy.cpp

MESH_OBJS := $(MESH_SRCS:%.cpp=$(BUILD_DIR)/%.o)
GRAPH_OBJS := $(GRAPH_SRCS:%.cpp=$(BUILD_DIR)/%.o)

PROGRAMS := naive_solver graph_solver

CXXFLAGS := $(MFEM_FLAGS) -I$(SRC_DIR)/meshing -I$(SRC_DIR)/graphs -MMD -MP

.PHONY: all clean

all: $(PROGRAMS)

naive_solver: $(BUILD_DIR)/naive_solver.o $(MESH_OBJS) $(MFEM_LIB_FILE)

graph_solver: $(BUILD_DIR)/graph_solver.o $(MESH_OBJS) $(GRAPH_OBJS) $(MFEM_LIB_FILE)

$(PROGRAMS):
	$(MFEM_CXX) $(MFEM_FLAGS) $(filter %.o,$^) -o $@ $(MFEM_LIBS)

# library sources from SRC_DIR
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(MFEM_CXX) $(CXXFLAGS) -c $< -o $@

# solver drivers in this directory
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(MFEM_CXX) $(CXXFLAGS) -c $< -o $@

-include $(wildcard $(BUILD_DIR)/*.d $(BUILD_DIR)/*/*.d)

clean:
	rm -rf $(BUILD_DIR) $(PROGRAMS)
