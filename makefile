MFEM_DIR = ${HOME}/mfem/mfem-4.8

ifndef MFEM_DIR
$(error MFEM_DIR is not set)
endif

MFEM_BUILD_DIR = $(MFEM_DIR)
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
INC = -I$(MFEM_DIR)

-include $(CONFIG_MK)

.PHONY: clean

pixels: pixels.cpp ppm.cpp make_mesh.cpp coarsen_mesh.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(INC) $< ppm.cpp make_mesh.cpp coarsen_mesh.cpp -o $@ -L$(MFEM_DIR) $(MFEM_LIBS)