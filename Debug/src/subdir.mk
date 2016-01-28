################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/1global.f90 \
../src/1tools.f90 \
../src/3qff.f90 \
../src/3states.f90 \
../src/4formPES.f90 \
../src/4integral.f90 \
../src/5modal.f90 \
../src/6oldvscf.f90 \
../src/6vscf.f90 \
../src/6xvscf.f90 \
../src/7vci.f90 \
../src/7vpt.f90 \
../src/main.f90 

OBJS += \
./src/1global.o \
./src/1tools.o \
./src/3qff.o \
./src/3states.o \
./src/4formPES.o \
./src/4integral.o \
./src/5modal.o \
./src/6oldvscf.o \
./src/6vscf.o \
./src/6xvscf.o \
./src/7vci.o \
./src/7vpt.o \
./src/main.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	gfortran -c -std=legacy -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace  -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


