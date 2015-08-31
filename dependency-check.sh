#!/bin/sh

myboola=`command -v intersectBed`
if [ ! -z "$myboola" ]
then
	echo "Verified that BEDTools is installed"
else
	echo "ERROR: BEDTools is NOT installed. Install BEDTools to proceed"
	exit 1
fi
myboolb=`command -v ./bigWigAverageOverBed`
if [ ! -z "$myboolb" ]
then
	echo "Verified that bigWigAverageOverBed"
else
	echo "ERROR: bigWigAverageOverBed is NOT installed. Install bigWigAverageOverBed to proceed"
	exit 1
fi
myboolc=`stat annotations/replication_timing.bw`
if [ ! -z "$myboolc" ]
then
	echo "Verified that DNA replication timing file is accounted for"
else
	echo "ERROR: Cannot find DNA replication timing file. Place a DNA replication timing bigWig file in the annotations folder to proceed. Expected name: \"replication_timing.bw\""
	exit 1
fi
myboold=`stat annotations/blacklist_regions.bed`
if [ ! -z "$myboold" ]
then
	echo "Verified that blacklist regions file is accounted for"
else
	echo "ERROR: Cannot find blacklist regions file. Place a blacklist regions file in the annotations folder to proceed. Expected name: \"blacklist_regions.bed\""
	exit 1
fi
myboole=`stat annotations/genes.bed`
if [ ! -z "$myboole" ]
then
	echo "Verified that gene annotation file is accounted for"
else
	echo "ERROR: Cannot find gene annotation file. Place a gene annotation file in the annotations folder to proceed. Expected name: \"genes.bed\""
	exit 1
fi
myboolf=`stat annotations/pgenes.bed`
if [ ! -z "$myboolf" ]
then
	echo "Verified that pseudogene annotation file is accounted for"
else
	echo "ERROR: Cannot find pseudogene annotation file. Place a pseudogene annotation file in the annotations folder to proceed. Expected name: \"pgenes.bed\""
	exit 1
fi

exit 0
