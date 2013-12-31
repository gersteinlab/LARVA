#!/bin/sh

# This script runs all LARVA test cases (future version will validate that the
# case produced the correct value)

testdir="../larva-test-cases"

for i in {1..14}
do
	echo "<== TEST CASE $i ==>"
	# echo "$testdir/t$i/qfiles.txt $testdir/t$i/dfiles.txt $testdir/t$i/test.db"
	perl larva-core.pl $testdir/t$i/qfiles.txt $testdir/t$i/dfiles.txt $testdir/t$i/test.db
done
