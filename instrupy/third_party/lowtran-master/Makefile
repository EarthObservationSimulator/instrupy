###############################################################################
# Project:        TATC for instrument module using Python
# Created by:     Joey Gurganus
# Date:           2019.04.17
###############################################################################
#
PROJECT := lowtran

.DEFAULT_GOAL := all

GFORTRAN_FLAG := $(shell which gfortran 2> /dev/null)
ifeq ($(GFORTRAN_FLAG),)
    $(error "No gfortran installed.   Please install it.")
endif

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  all        to perform clean-up and installation"
	@echo "  install    to set up the python package (pip install -e .)"
	@echo "  runtest    to perform unit testing"
	@echo "  testlog    to perform unit testing with no log capture"
	@echo "  fulltest   to perform unit testing with verbose"
	@echo "  clean      to remove *.pyc files and __pycache__ directories"
	@echo "  bare       to uninstall the package and remove *egg*"

all: bare install

install:
	pip install -e .

runtest:
	@nosetests

testlog:
	@nosetests --nologcapture

fulltest:
	@nosetests --nologcapture -v

clean:
	@echo "Cleaning up..."
	@find . -name "*.pyc" -delete
	@find . -type d -name __pycache__ -print0 | xargs -0 rm -rf

bare: clean
	pip uninstall -y $(PROJECT) 
	rm -rf $(PROJECT).egg-info .eggs
