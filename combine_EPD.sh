#!/bin/bash

cd output
mkdir -p temp1
mkdir -p temp2
cp ../hadd2.pl ./temp1/
cp ../hadd2.pl ./temp2/

cp cent_EPD_CorrectionOutput_* ./temp1/
cd temp1
perl hadd2.pl
rm cent_EPD_CorrectionOutput_*
perl hadd2.pl
cd ../


