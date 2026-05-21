-include ../Make.user

ifndef MFEM_DIR
$(error MFEM_DIR is not set)
endif

ifndef SRC_DIR
$(error SRC_DIR is not set)
endif

MFEM_BUILD_DIR = $(MFEM_DIR)
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
-include $(CONFIG_MK)

BUILD_DIR := build

SRCS := \
	graph_solver.cpp \
	$(SRC_DIR)/meshing/ppm.cpp \
	$(SRC_DIR)/meshing/pixel_mesh.cpp \
    $(SRC_DIR)/graphs/voxel_graph.cpp \
	$(SRC_DIR)/graphs/voxel_graph_operator.cpp \
	$(SRC_DIR)/graphs/voxel_graph_hierarchy.cpp

OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

DEPS := $(OBJS:.o=.d)

INCLUDES := $(MFEM_INCFLAGS) -I$(SRC_DIR)/meshing -I$(SRC_DIR)/graphs
CXXFLAGS := $(MFEM_CXXFLAGS) -fsanitize=address $(INCLUDES) -MMD -MP -g


all: dirichlet_2d graph_solver

dirichlet_2d: dirichlet_2d.cpp $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) -fsanitize=address $(MFEM_INCFLAGS) $< $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp -o $@ -L$(MFEM_DIR) $(MFEM_LIBS) -I$(SRC_DIR)/meshing

graph_solver: $(OBJS)
	$(MFEM_CXX) $(MFEM_CXXFLAGS) -fsanitize=address $^ -o $@ -L$(MFEM_DIR) $(MFEM_LIBS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(MFEM_CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

.PHONY: clean

clean:
	rm -rf $(BUILD_DIR) dirichlet_2d graph_solver
