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

.PHONY: clean

dirichlet_2d: dirichlet_2d.cpp $(SRC_DIR)/meshing_2d/ppm.cpp $(SRC_DIR)/meshing_2d/pixel_mesh.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS) $< $(SRC_DIR)/meshing_2d/ppm.cpp $(SRC_DIR)/meshing_2d/pixel_mesh.cpp -o $@ -L$(MFEM_DIR) $(MFEM_LIBS) -I$(SRC_DIR)/meshing_2d

clean:
	rm -f dirichlet_2d
