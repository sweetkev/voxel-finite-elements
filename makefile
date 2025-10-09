MFEM_DIR = ${HOME}/mfem/mfem-4.8

ifndef MFEM_DIR
$(error MFEM_DIR is not set)
endif

MFEM_BUILD_DIR = $(MFEM_DIR)
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk

-include $(CONFIG_MK)

.PHONY: clean

make_mesh: make_mesh.cpp
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $< -o $@ -L$(MFEM_DIR) $(MFEM_LIBS) 