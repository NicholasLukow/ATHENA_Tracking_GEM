#/bin/bash

etamin=1.0
etamax=3.5
etainc=0.5
eta=$etamin


while [ "$(bc <<< "$eta < $etamax")" == "1" ]; do
pmin=1
pmax=20
pinc=1
p=$pmin

	while [ "$(bc <<< "$p < $pmax")" == "1" ] ; do
		ptemp=$(bc <<< "$p+$pinc")
		etatemp=$(bc <<< "$eta+$etainc")
		
		#Output Root file name identifies the parameters of the generator
		#root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$p', '$ptemp','$eta', '$etatemp', 4, "P_'$p'_'$ptemp'_Eta_'$eta'_'$etatemp'_")'
		root -b -l -q './RunAnalysis.C('$p', '$ptemp','$eta', '$etatemp')'
		#echo Eta: $eta P: $p 


		p=$(bc <<< "$p+$pinc")
		
		#incrementing by 2 above 10
		if [ "$(bc <<< "$p == 10.0")" == "1" ]; then
			pinc=2
		fi

	done
	eta=$(bc <<< "$eta+$etainc")

done
