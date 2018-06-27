################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../Header.f90 \
../Main.f90 \
../routines.f90 \
../routines_generic.f90 

OBJS += \
./Header.o \
./Main.o \
./routines.o \
./routines_generic.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -cpp  -fopenmp -c -fmessage-length=0  -ffree-line-length-none -std=f2008 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

Header.o: ../Header.f90

Main.o: ../Main.f90 Header.o routines.o

routines.o: ../routines.f90 Header.o routines_generic.o

routines_generic.o: ../routines_generic.f90 Header.o


