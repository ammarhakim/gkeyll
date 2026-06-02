# -*- makefile-gmake -*-
# Generated automatically! Do not edit

.PHONY:
all: core core-unit core-regression ## Build only specified Apps (core)

.PHONY: install
install: core-install ## Install gkeyll executable
	cd gkeyll && ${MAKE} -f Makefile-gkeyll install

