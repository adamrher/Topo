# Find a usable python command
ifndef python_command
  python_command := $(shell command -v  $(HOME)/anaconda/bin/python 2> /dev/null)
endif
ifndef python_command
  python_command := $(shell command -v  python 2> /dev/null)
endif

# Attempt to identify machine for special support
MACNAME = generic_mac
ifeq ($(shell uname -s),Darwin)
  # Special case for Macs -- ignore machine name
  machname := $(MACNAME)
else
  machname := $(HOSTNAME)
  ifndef machname
    machname := $(shell uname -n)
  endif
  # Special case for Yellowstone
  machine := $(patsubst yslogin%,yellowstone,$(machname))
endif

# Find a suitable compiler -- f77 is not suitable
ifeq ($(origin FC),default)
  undefine FC
endif

ifndef FC
  FC := gfortran
  ifeq ($(machine),yellowstone)
    FC := ifort
  endif
  export FC
endif

#
# compiler settings
#

#
#------------------------------------------------------------------------
# ifort
#------------------------------------------------------------------------
#
ifeq ($(FC),ifort)
  INC_NETCDF :=${NETCDF}/include
  LIB_NETCDF :=${NETCDF}/lib

  FFLAGS = -c -g -r8 -O1 -I$(INC_NETCDF)
  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
#  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff

  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -Mbounds -traceback -Mchkfpstk
  else
    FFLAGS += -O
  endif
endif


#------------------------------------------------------------------------
# GFORTRAN
#------------------------------------------------------------------------
#
ifeq ($(FC),gfortran)
  ifeq ($(machine),$(MACNAME))
    INC_NETCDF := /opt/local/include
    LIB_NETCDF := /opt/local/lib
  endif
  ifeq ($(machine),harmon)
    INC_NETCDF :=/usr/local/netcdf-gcc-g++-gfortran/include
    LIB_NETCDF :=/usr/local/netcdf-gcc-g++-gfortran/lib
  endif

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff 
  FFLAGS   := -c  -fdollar-ok  -I$(INC_NETCDF)

  ifeq ($(DEBUG),TRUE)
#   FFLAGS += --chk aesu  -Cpp --trace
    FFLAGS += -Wall -fbacktrace -fbounds-check -fno-range-check
  else
    FFLAGS += -O
  endif

endif

#------------------------------------------------------------------------
# NAG
#------------------------------------------------------------------------
ifeq ($(FC),nagfor)

#  INC_NETCDF :=/usr/local/netcdf-gcc-nag/include
#  LIB_NETCDF :=/usr/local/netcdf-pgi/lib

  INC_NETCDF :=/usr/local/netcdf-gcc-nag/include
  LIB_NETCDF :=/usr/local/netcdf-gcc-nag/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff
  FFLAGS   := -c  -I$(INC_NETCDF)

  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -C
  else
    FFLAGS += -O
  endif

endif

#------------------------------------------------------------------------
# PGF95
#------------------------------------------------------------------------
#
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib
#

ifeq ($(FC),pgf95)
  INC_NETCDF :=/opt/local/include
  LIB_NETCDF :=/opt/local/lib

  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff
  FFLAGS   := -c -Mlarge_arrays -I$(INC_NETCDF)


  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -Mbounds -traceback -Mchkfpstk
  else
    FFLAGS += -O
  endif

endif

export INC_NETCDF
export LIB_NETCDF
export LDFLAGS
export FFLAGS

comp_check:
ifndef python_command
	$(error "python not found, please ensure python is in your PATH")
endif
	echo "python = $(python_command)"
	echo "machine = '$(machine)'"
	echo "FC = $(FC)"
	echo $(origin FC)
