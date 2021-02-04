#!/bin/bash

# Select input file:
# ===================================================
INPUT_FILE="Bfield_b.txt"
NZ_STR=501

# Set number of threads:
# ===================================================
OMP_NUM_THREADS=4

# Set processor binding for openMP:
# ===================================================
OMP_PROC_BIND=true

# Get repo directory:
# ===================================================
REPO_DIR=$(echo $PWD | sed "s|/src||g")

# Get input file directory:
# ===================================================
INPUT_FILE_DIR=$REPO_DIR/BfieldData/$INPUT_FILE

# Make environment variables available to shell:
# ===================================================
export REPO_DIR
export INPUT_FILE_DIR
export INPUT_FILE
export NZ_STR
export OMP_NUM_THREADS
export OMP_PROC_BIND

# Compile source code:
# ===================================================
make -f $REPO_DIR/src/make_test_spline.f
mv $REPO_DIR/src/test_spline $REPO_DIR/bin/

# Run code:
# ===================================================
if [ $? -eq 0 ] ; then
	echo ''
	echo '**********************'
	echo 'Succesful compilation!'
	echo ''
	echo 'Code is running...'
	echo '**********************'
	echo ''
	$REPO_DIR/bin/test_spline
	echo ''
	echo 'Calculation complete!'
	echo ''
else
	echo ''
	echo 'Compilation error'
	echo ''
fi
