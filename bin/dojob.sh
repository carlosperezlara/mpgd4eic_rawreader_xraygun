#!/bin/sh
export N=$(($1 + 1))
export A=`head -$N runs.dat| tail -1`
echo $A
B=-1
echo "./datamonitor test0001/dream-0000${A}-0000.evt"
./datamonitor test0001/dream-0000${A}-0000.evt
