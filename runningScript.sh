#/bin/bash
nEvents=5000
count=0
#for pmin in 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5  ; do 
#root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$pmin', '$pmin'+1,'$etamin', '$etamin'+0.25, 4, "job'$count'")'
#let count=count+1
#done
#for pmin in 11 13 15 17 19  ; do 
#root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$pmin', '$pmin'+2,'$etamin', '$etamin'+0.25, 4, "job'$count'")'
#let count=count+1
#done

pmin=1
pmax=20
for etamin in 0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 ; do 
root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$pmin', '$pmax','$etamin', '$etamin'+0.25, 4, "job'$count'_NoGEM")'
let count=count+1
done
