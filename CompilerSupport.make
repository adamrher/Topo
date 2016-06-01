# Find a usable python command
ifndef python_command
  python_command := $(shell command -v  $(HOME)/anaconda/bin/python 2> /dev/null)
endif
ifndef python_command
  python_command := $(shell command -v  python 2> /dev/null)
endif

# Attempt to identify machine for special support
machname := $(HOSTNAME)
ifndef machname
  machname := $(shell uname -n)
endif
# Special case for Yellowstone
machine := $(patsubst yslogin%,yellowstone,$(machname))

# Find a suitable compiler -- f77 is not suitable
ifeq ($(origin FC),default)
  undefine FC
endif

ifndef FC
  FC := gfortran
  ifeq ($(machine),yellowstone)
    FC := ifort
  endif
endif

all:
ifndef python_command
	$(error "python not found, please ensure python is in your PATH")
endif
	echo "python = $(python_command)"
	echo "machine = '$(machine)'"
	echo "FC = $(FC)"
	echo $(origin FC)
