#!/bin/bash

# This is a script to write a string as the prefix in a csv-file.

# $1 = prefix string
# $2 = input file
# $3 = output file

prefix=$1

sed "s/^/${prefix}/g" $2 >> $3
