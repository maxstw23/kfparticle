#!/bin/bash

cd output
mkdir -p temp1
mkdir -p temp2
cp ../hadd2.pl ./temp1/
cp ../hadd2.pl ./temp2/

cp EPDShiftOutput_* ./temp1/
cd temp1
perl hadd2.pl
rm EPDShiftOutput_*
perl hadd2.pl
cd ../


