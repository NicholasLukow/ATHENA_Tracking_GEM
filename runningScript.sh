#/bin/bash
nEvents=5000


#detversion="noGem"
#detversion="1ffg"
detversion="2ffg"

for BField in 4 5 ; do 

	etamin=-3.5
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
			root -b -l -q './Fun4All_G4_HybridGEM.C('$nEvents' , '$p', '$ptemp','$eta', '$etatemp', '$BField', "P_'$p'_'$ptemp'_Eta_'$eta'_'$etatemp'_'$detversion'")'
			p=$(bc <<< "$p+$pinc")
		
			#incrementing by 2 above 10
			if [ "$(bc <<< "$p == 10.0")" == "1" ]; then
				pinc=2
			fi

		done
		eta=$(bc <<< "$eta+$etainc")

	done
done
