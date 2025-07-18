# -*- makefile-gmake -*-

# Type "make help" to see help for this Makefile

# determine date of build
BUILD_DATE := $(shell date)
GIT_TIP := $(shell git describe --abbrev=12 --always --dirty=+)

# Build directory
BUILD_DIR ?= build

# Name of kernels directory
KERNELS_DIR := ker

ARCH_FLAGS ?= -march=native
CUDA_ARCH ?= 70
CFLAGS ?= -O3 -g -ffast-math -fPIC -MMD -MP -DGIT_COMMIT_ID=\"$(GIT_TIP)\" -DGKYL_BUILD_DATE="${BUILD_DATE}" -DGKYL_GIT_CHANGESET="${GIT_TIP}"
LDFLAGS = 
PREFIX ?= ${HOME}/gkylsoft

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = $(PREFIX)/OpenBLAS/include
LAPACK_LIB_DIR = $(PREFIX)/OpenBLAS/lib
LAPACK_LIB = -lopenblas

# Include config.mak file (if it exists) to overide defaults above
-include config.mak

INSTALL_PREFIX ?= ${PREFIX}
PROJ_NAME ?= gkeyll

# Determine OS we are running on
UNAME = $(shell uname)

# Read ADAS paths and flags if needed 
USING_ADAS =
ADAS_INC_DIR = zero # dummy
ADAS_LIB_DIR = .
ifeq (${USE_ADAS}, 1)
	USING_ADAS = yes
	CFLAGS += -DGKYL_HAVE_ADAS
endif

# Directory for storing shared data, like ADAS reaction rates and radiation fits
GKYL_SHARE_DIR ?= "${INSTALL_PREFIX}/${PROJ_NAME}/share"
CFLAGS += -DGKYL_SHARE_DIR=$(GKYL_SHARE_DIR)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_LIB_DIR = .
	LAPACK_INC = core # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
endif

# CUDA flags
USING_NVCC =
NVCC_FLAGS = 
CUDA_LIBS =
ifeq ($(CC), nvcc)
	USING_NVCC = yes
	CFLAGS = -O3 -g --forward-unknown-to-host-compiler --use_fast_math -ffast-math -MMD -MP -fPIC -DGIT_COMMIT_ID=\"$(GIT_TIP)\" -DGKYL_BUILD_DATE="${BUILD_DATE}" -DGKYL_GIT_CHANGESET="${GIT_TIP}"
	NVCC_FLAGS = -x cu -dc -arch=sm_${CUDA_ARCH} -rdc=true --compiler-options="-fPIC"
	LDFLAGS += -arch=sm_${CUDA_ARCH} -rdc=true
	ifdef CUDAMATH_LIBDIR
		CUDA_LIBS = -L${CUDAMATH_LIBDIR}
	else
		CUDA_LIBS =
	endif
	CUDA_LIBS += -lcublas -lcusparse -lcusolver
	SQL_CFLAGS = --forward-unknown-to-host-compiler -fPIC
endif

# MPI paths and flags
USING_MPI =
MPI_RPATH = 
MPI_INC_DIR = core # dummy
MPI_LIB_DIR = .
ifeq (${USE_MPI}, 1)
	USING_MPI = yes
	MPI_INC_DIR = ${CONF_MPI_INC_DIR}
	MPI_LIB_DIR = ${CONF_MPI_LIB_DIR}

ifdef USING_NVCC
	MPI_RPATH = -Xlinker "-rpath,${CONF_MPI_LIB_DIR}"
else
	MPI_RPATH = -Wl,-rpath,${CONF_MPI_LIB_DIR}
endif

	MPI_LIBS = -lmpi
	CFLAGS += -DGKYL_HAVE_MPI
endif

# Read NCCL paths and flags if needed (needs MPI and NVCC)
USING_NCCL =
NCCL_INC_DIR = zero # dummy
NCCL_LIB_DIR = .
ifeq (${USE_NCCL}, 1)
ifdef USING_MPI
ifdef USING_NVCC
	USING_NCCL = yes
	NCCL_INC_DIR = ${CONF_NCCL_INC_DIR}
	NCCL_LIB_DIR = ${CONF_NCCL_LIB_DIR}
	NCCL_LIBS = -lnccl
	CFLAGS += -DGKYL_HAVE_NCCL
endif
endif
endif

# Read CUDSS paths and flags if needed (needs MPI and NVCC)
USING_CUDSS =
CUDSS_INC_DIR = zero # dummy
CUDSS_LIB_DIR = .
CUDSS_RPATH =
ifeq (${USE_CUDSS}, 1)
ifdef USING_NVCC
	USING_CUDSS = yes
	CUDSS_INC_DIR = ${CONF_CUDSS_INC_DIR}
	CUDSS_LIB_DIR = ${CONF_CUDSS_LIB_DIR}
	CUDSS_RPATH = -Xlinker "-rpath,${CONF_CUDSS_LIB_DIR}"
	CUDSS_LIBS = -lcudss
	CFLAGS += -DGKYL_HAVE_CUDSS
endif
endif

# LUA paths and flags
USING_LUA =
LUA_RPATH = 
LUA_INC_DIR = core # dummy
LUA_LIB_DIR = .
ifeq (${USE_LUA}, 1)
	USING_LUA = yes
	LUA_INC_DIR = ${CONF_LUA_INC_DIR}
	LUA_LIB_DIR = ${CONF_LUA_LIB_DIR}

ifdef USING_NVCC
	LUA_RPATH = -Xlinker "-rpath,${CONF_LUA_LIB_DIR}"
else
	LUA_RPATH = -Wl,-rpath,${CONF_LUA_LIB_DIR}
endif

	LUA_LIBS = -l${CONF_LUA_LIB}
	CFLAGS += -DGKYL_HAVE_LUA
endif

# Build directory
ifdef USING_NVCC
	BUILD_DIR = cuda-build
endif

# Command to make dir
MKDIR_P ?= mkdir -p

# At this point, export all top-level variables to sub-makes and
# recurse downwards

.EXPORT_ALL_VARIABLES:

# Regression tests
${BUILD_DIR}/core/creg/%: core/creg/%.c ${BUILD_DIR}/core/libg0core.so
	cd core && $(MAKE) -f Makefile-core ../$@

${BUILD_DIR}/moments/creg/%: moments/creg/%.c ${BUILD_DIR}/moments/libg0moments.so
	cd moments && $(MAKE) -f Makefile-moments ../$@

${BUILD_DIR}/vlasov/creg/%: vlasov/creg/%.c ${BUILD_DIR}/vlasov/libg0vlasov.so
	cd vlasov && $(MAKE) -f Makefile-vlasov ../$@

${BUILD_DIR}/gyrokinetic/creg/%: gyrokinetic/creg/%.c ${BUILD_DIR}/gyrokinetic/libg0gyrokinetic.so
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic ../$@

${BUILD_DIR}/pkpm/creg/%: pkpm/creg/%.c ${BUILD_DIR}/pkpm/libg0pkpm.so
	cd pkpm && $(MAKE) -f Makefile-pkpm ../$@

# Unit tests
${BUILD_DIR}/core/unit/%: core/unit/%.c ${BUILD_DIR}/core/libg0core.so
	cd core && $(MAKE) -f Makefile-core ../$@

${BUILD_DIR}/moments/unit/%: core/moments/%.c ${BUILD_DIR}/moments/libg0moments.so
	cd moments && $(MAKE) -f Makefile-moments ../$@

${BUILD_DIR}/vlasov/unit/%: vlasov/unit/%.c ${BUILD_DIR}/vlasov/libg0vlasov.so
	cd vlasov && $(MAKE) -f Makefile-vlasov ../$@

${BUILD_DIR}/gyrokinetic/unit/%: gyrokinetic/unit/%.c ${BUILD_DIR}/gyrokinetic/libg0gyrokinetic.so
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic ../$@

${BUILD_DIR}/pkpm/unit/%: pkpm/unit/%.c ${BUILD_DIR}/pkpm/libg0pkpm.so
	cd pkpm && $(MAKE) -f Makefile-pkpm ../$@

# Declare sub-directories as phony targets
.PHONY: core moments vlasov gyrokinetic pkpm

all: gkeyll
	${MKDIR_P} ${INSTALL_PREFIX}/${PROJ_NAME}/share/adas

everything: regression unit gkeyll ## Build everything, including unit, regression and gkeyll exectuable

## Core infrastructure targets
core:  ## Build core infrastructure code
	cd core && $(MAKE) -f Makefile-core

core-unit: ## Build core unit tests
	cd core && $(MAKE) -f Makefile-core unit

core-regression: ## Build core regression tests
	cd core && $(MAKE) -f Makefile-core regression

core-install: ## Install core infrastructure code
	cd core && $(MAKE) -f Makefile-core install
	test -e config.mak && cp -f config.mak ${INSTALL_PREFIX}/${PROJ_NAME}/share/config.mak || echo "No config.mak"

core-clean: ## Clean core infrastructure code
	cd core && $(MAKE) -f Makefile-core clean

core-check: core ## Run unit tests in core
	cd core && $(MAKE) -f Makefile-core check

core-valcheck: core ## Run valgrind on unit tests in core
	cd core && $(MAKE) -f Makefile-core valcheck

## Moments infrastructure targets
moments: core  ## Build moments infrastructure code
	cd moments && $(MAKE) -f Makefile-moments

moments-unit: moments ## Build moments unit tests
	cd moments && $(MAKE) -f Makefile-moments unit

moments-regression: moments ## Build moments regression tests
	cd moments && $(MAKE) -f Makefile-moments regression

moments-amr-regression: moments ## Build moments AMR regression tests
	cd moments && $(MAKE) -f Makefile-moments amr_regression

moments-install: core-install ## Install moments infrastructure code
	cd moments && $(MAKE) -f Makefile-moments install

moments-clean: ## Clean moments infrastructure code
	cd moments && $(MAKE) -f Makefile-moments clean

moments-check: moments ## Run unit tests in moments
	cd moments && $(MAKE) -f Makefile-moments check

moments-valcheck: moments ## Run valgrind on unit tests in moments
	cd moments && $(MAKE) -f Makefile-moments valcheck

## Vlasov infrastructure targets
vlasov: moments  ## Build Vlasov infrastructure code
	cd vlasov && $(MAKE) -f Makefile-vlasov

vlasov-unit: vlasov ## Build Vlasov unit tests
	cd vlasov && $(MAKE) -f Makefile-vlasov unit

vlasov-regression: vlasov ## Build Vlasov regression tests
	cd vlasov && $(MAKE) -f Makefile-vlasov regression

vlasov-install: moments-install ## Install Vlasov infrastructure code
	cd vlasov && $(MAKE) -f Makefile-vlasov install

vlasov-clean: ## Clean Vlasov infrastructure code
	cd vlasov && $(MAKE) -f Makefile-vlasov clean

vlasov-check: vlasov ## Run unit tests in Vlasov
	cd vlasov && $(MAKE) -f Makefile-vlasov check

vlasov-valcheck: vlasov ## Run valgrind on unit tests in Vlasov
	cd vlasov && $(MAKE) -f Makefile-vlasov valcheck

## Gyrokinetic infrastructure targets
gyrokinetic: vlasov  ## Build Gyrokinetic infrastructure code
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic

gyrokinetic-unit: gyrokinetic ## Build Gyrokinetic unit tests
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic unit

gyrokinetic-regression: gyrokinetic ## Build Gyrokinetic regression tests
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic regression

gyrokinetic-install: vlasov-install ## Install Gyrokinetic infrastructure code
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic install

gyrokinetic-clean: ## Clean Gyrokinetic infrastructure code
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic clean

gyrokinetic-check: gyrokinetic ## Run unit tests in Gyrokinetics
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic check

gyrokinetic-valcheck: gyrokinetic ## Run valgrind on unit tests in Gyrokinetics
	cd gyrokinetic && $(MAKE) -f Makefile-gyrokinetic valcheck

## PKPM infrastructure targets
pkpm: gyrokinetic  ## Build PKPM infrastructure code
	cd pkpm && $(MAKE) -f Makefile-pkpm

pkpm-unit: pkpm ## Build PKPM unit tests
	cd pkpm && $(MAKE) -f Makefile-pkpm unit

pkpm-regression: pkpm ## Build PKM regression tests
	cd pkpm && $(MAKE) -f Makefile-pkpm regression

pkpm-install: gyrokinetic-install ## Install PKPM infrastructure code
	cd pkpm && $(MAKE) -f Makefile-pkpm install

pkpm-clean: ## Clean PKPM infrastructure code
	cd pkpm && $(MAKE) -f Makefile-pkpm clean

pkpm-check: pkpm ## Run unit tests in PKPM
	cd pkpm && $(MAKE) -f Makefile-pkpm check

pkpm-valcheck: pkpm ## Run valgrind on unit tests in PKPM
	cd pkpm && $(MAKE) -f Makefile-pkpm valcheck

## Top-level Gkeyll target
gkeyll: pkpm
	cd gkeyll && ${MAKE} -f Makefile-gkeyll gkeyll

gkeyll-install: pkpm-install gkeyll
	cd gkeyll && ${MAKE} -f Makefile-gkeyll install

## Targets to build things all parts of the code

# build all unit tests 
unit: pkpm-unit gyrokinetic-unit vlasov-unit moments-unit core-unit ## Build all unit tests

# build all regression tests 
regression: pkpm-regression gyrokinetic-regression vlasov-regression moments-regression core-regression ## Build all regression tests

# Install everything
install: gkeyll-install  ## Install all code

# Clean everything
clean:
	rm -rf ${BUILD_DIR}

# Check everything
check: core-check moments-check vlasov-check gyrokinetic-check pkpm-check ## Run all unit tests

# From: https://www.client9.com/self-documenting-makefiles/
.PHONY: help
help: ## Show help
	@echo "Gkeyll Makefile help. You can set parameters on the command line:"
	@echo ""
	@echo "make CC=cc -j"
	@echo ""
	@echo "Or run the configure script to set various parameters. Usually"
	@echo "defaults are all you need, specially if the dependencies are in"
	@echo "${HOME}/gkylsoft and you are using standard compilers (not building on GPUs)."
	@echo ""
	@echo "See ./configure --help for usage of configure script."
	@echo ""
	@echo "You can build only portions of the code using the specific targers below."
	@echo "Typing \"make all\" will build the complete code"
	@echo ""
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {\
        printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF \
        }' $(MAKEFILE_LIST)
