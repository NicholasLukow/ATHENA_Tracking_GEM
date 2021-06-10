#/bin/bash
nEvents=20000

etamin=0
etamax=3.5
		
pmin=1
pmax=10

BField=5

#detversion="noGem"
#detversion="1ffg"
detversion="2ffg"


root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$pmin', '$pmax','$etamin', '$etamax', '$BField', "'$detversion'")'
