# -*- makefile-gmake -*-

# Moments include paths
MOMENTS_INCS = -I../moments/$(KERNELS_DIR)/fem_poisson -I../moments/apps -I../moments/zero

# Moments link objects

MOMENTS_LINK_OBJS = ../$(BUILD_DIR)/moments/zero/*.o ../$(BUILD_DIR)/moments/$(KERNELS_DIR)/*/*.o ../$(BUILD_DIR)/moments/apps/*.o ../$(BUILD_DIR)/moments/amr/*.o
# ... with GPUs
MOMENTS_LINK_CU_OBJS = ../$(BUILD_DIR)/moments/unit/*.cu.o
