#!/bin/bash

# Select input file
INPUT_FILE="xp_Xray.in"

# Populate the correct directory information in the input file:
REPO_DIR=$(echo $PWD | sed "s|/src||g")
INPUT_FILE=$REPO_DIR/InputFiles/$INPUT_FILE
sed -i "s|REPO_DIR|"$REPO_DIR"|g" $INPUT_FILE

# Compile:
make -f $REPO_DIR/src/Makefile_openMP.f

# Run code:
if [ $? -eq 0 ] ; then
echo 'Succesfull compilation'
echo 'Code is running...'
./MPEX
echo 'Calculation complete!'
else
echo 'Compilation error'
fi
