# -*- makefile-gmake -*-

# Pkpm link objects

PKPM_LINK_OBJS = ../$(BUILD_DIR)/pkpm/zero/*.o ../$(BUILD_DIR)/pkpm/apps/*.o ../$(BUILD_DIR)/pkpm/$(KERNELS_DIR)/*/*.o
# ... with GPUs
PKPM_LINK_CU_OBJS = 
