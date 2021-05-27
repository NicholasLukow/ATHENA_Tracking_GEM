#/bin/bash
nEvents=5000
count=0
for pmin in 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5  ; do 
root -b -l -q './Fun4All_G4_ModifiedSi.C('$nEvents' , '$pmin', '$pmin'+1, 4, "job'$count'")'
let count=count+1
done
for pmin in 11 13 15 17 19  ; do 
root -b -l -q './Fun4All_G4_ModifiedSi.C('$nEvents' , '$pmin', '$pmin'+2, 4, "job'$count'")'
let count=count+1
done
