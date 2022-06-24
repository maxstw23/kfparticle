#!/bin/csh

#input begin
set energy=27.0
set nRun=18
set mRun=Run$nRun
set production=P19ib
set library=SL19b
set trgsetupname=27GeV_production_2018
#input end

#run18List27.0
#set outList=./$energy.list
set outList=run${nRun}List${energy}.list
if(-e $outList)rm -r $outList

#get_file_list.pl -keys 'runnumber' -delim '/' -cond 'production=P19ib,filetype=daq_reco_PicoDst,trgsetupname=27GeV_production_2018,tpx=1,filename~st_physics,sanity=1,storage!=hpss' -limit 0 > $outList
get_file_list.pl -keys 'runnumber' -delim '/' -cond "production=$production,filetype=daq_reco_PicoDst,trgsetupname~$trgsetupname,runnumber[]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics" -limit 0 > $outList
