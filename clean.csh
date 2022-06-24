#!/bin/csh

echo "Do you want to delete all the output/log/xml? [Y/n]"
if( !($< =~ [Yy]*) ) then
	echo "Nothing happened."
	exit
else
	rm ./output/*
	rm ./log/*
	rm ./xml/*
endif


