#!/bin/bash

# Select input file:
# ===================================================
INPUT_FILE="xp_TestCollisions.in"
#INPUT_FILE="xp_IonSlowingDown.in"
#INPUT_FILE="xp_Test.in"

# Set number of threads:
# ===================================================
export OMP_NUM_THREADS=1

# Set processor binding for openMP:
# ===================================================
export OMP_PROC_BIND=true

# Get repo directory:
# ===================================================
REPO_DIR=$(echo $PWD | sed "s|/src||g")

# Get input file directory:
# ===================================================
INPUT_FILE_DIR=$REPO_DIR/InputFiles/$INPUT_FILE

# Populate directory information in the input file:
# ===================================================
sed -i "s|REPO_DIR|"$REPO_DIR"|g" $INPUT_FILE_DIR

# Make environment variables available to shell:
# ===================================================
export REPO_DIR
export INPUT_FILE_DIR
export INPUT_FILE

# Compile source code:
# ===================================================
make -f $REPO_DIR/src/make_linFP.f
mv $REPO_DIR/src/linFP $REPO_DIR/bin/

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
	$REPO_DIR/bin/linFP
	echo ''
	echo 'Calculation complete!'
	echo ''
else
	echo ''
	echo 'Compilation error'
	echo ''
fi
