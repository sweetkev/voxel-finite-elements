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

all: dirichlet_2d graph_solver_2d

dirichlet_2d: dirichlet_2d.cpp $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS) $< $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp -o $@ -L$(MFEM_DIR) $(MFEM_LIBS) -I$(SRC_DIR)/meshing

graph_solver: graph_solver.cpp $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp $(SRC_DIR)/graphs/graph.cpp $(SRC_DIR)/graphs/graph_operator.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS) $< $(SRC_DIR)/meshing/ppm.cpp $(SRC_DIR)/meshing/pixel_mesh.cpp $(SRC_DIR)/graphs/graph.cpp $(SRC_DIR)/graphs/graph_operator.cpp -o $@ -L$(MFEM_DIR) $(MFEM_LIBS) -I$(SRC_DIR)/meshing -I$(SRC_DIR)/graphs

.PHONY: clean

clean:
	rm -f dirichlet_2d graph_solver_2d
