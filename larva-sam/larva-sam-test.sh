#!/bin/sh

# This script runs all LARVA rand pipeline test cases (future version will validate
# that the case produced the correct value)

testdir="../../larva-sam-test-cases"

for i in {1..2}
do
	echo "<== TEST CASE $i ==>"
	perl larva-sam.pl $testdir/t$i/vfiles.txt $testdir/t$i/afiles.txt 8 8 e $testdir/t$i/test_db.db
done

for i in {3..5} # Normal range: 3..5
do
# 	echo "<== TEST CASE $i (serial) ==>"
# 	perl larva-sam.pl $testdir/t$i/vfiles.txt $testdir/t$i/afiles.txt 8 1 e $testdir/t$i/test_db.db
	echo "<== TEST CASE $i (parallel) ==>"
	perl larva-sam.pl $testdir/t$i/vfiles.txt $testdir/t$i/afiles.txt 8 8 e $testdir/t$i/test_db.db
done
