#!/bin/bash

for i in {2..4}
do
  ./admixture --cv ../../poecileADMIXTURE.bed $i > admixture${i}.out
done
