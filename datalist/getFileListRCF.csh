#!/bin/csh

#input begin
set energy=27
set nRun=18

#run18List27
set outListRaw=run${nRun}List${energy}_raw.list
if(-e $outListRaw)rm -r $outListRaw
set outList=run${nRun}List${energy}.list
if(-e $outList)rm -r $outList

get_file_list.pl -keys path,filename -delim '/' -cond "production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber[]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics" -limit 0 > $outListRaw
awk '{printf "root://xrdstar.rcf.bnl.gov:1095/%s", $0}' $outListRaw > $outList