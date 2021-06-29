#/bin/bash

declare -a eta_RangeValues
declare -a mom_RangeValues

etabins=1
#eta_RangeValues=( -3.5 -2.5 -1.0 1.0 2.5 3.5)
#eta_RangeValues=( -1.0 1.0)
eta_RangeValues=( 1.0 2.5 3.5)

mombins=29
mom_RangeValues=(1 30)

#0 to use momentum binning, 1 to use pt binning
PtOption=0

#Variable to be plotted on Y-axis
#Yaxis: 0=momentum resolution  1=transverse pointing resolution  2=tracking efficiency  3=projection resolutions (requires further choices below)  4=longitudinal pointing resolution
Yaxis=4
#if Yaxis=3, the following must be specified
#ProjectionSurface:  0=DIRC  1=RICH  2=EMCal
ProjectionSurface=1
#ProjectionValue: 0=Pointing resolution on surface  1: Theta resolution  2: Phi resolution
ProjectionValue=2

#Fit the plot with the functional form specified by the PWG requirement: 0 is off, 1 is on
DrawFits=0
#Draw the PWG requirements if the specifed range and plotting variables match one from table 11.2 in Yellow Report: 0 is off, 1 is on
DrawPWGReq=0
	
#for etaindex in {0..4..1}; do
for etaindex in {0..1..1}; do
	for momindex in {0..0..1}; do
	
		root -b -q -l './Analysis.C('${etabins}','${eta_RangeValues[$etaindex]}','${eta_RangeValues[$etaindex+1]}','${mombins}','${mom_RangeValues[$momindex]}','${mom_RangeValues[$momindex+1]}','${PtOption}','${Yaxis}','${ProjectionSurface}','${ProjectionValue}','${DrawFits}','${DrawPWGReq}')'

	done
done	
