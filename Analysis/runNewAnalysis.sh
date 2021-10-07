#/bin/bash

declare -a eta_RangeValues
declare -a mom_RangeValues

etabins=1
#eta_RangeValues=( -3.5 -2.5 -1.0 1.0 2.5 3.5)
#eta_RangeValues=( 2.5 3.5)
#eta_RangeValues=( -3.5 3.5)
#eta_RangeValues=( -2.5 -1.0)
#eta_RangeValues=( -2.25 -1.25)
#eta_RangeValues=( -3.5 -2.5)
#eta_RangeValues=( -3.5 -2.5 -1.0)
#eta_RangeValues=( 1.0 2.5 3.5 )
#eta_RangeValues=( 1.0 2.5)
#eta_RangeValues=( -3.5 -3.0 -2.5 -1.0 1.0 2.5 3.0 3.5)

#eta_RangeValues=(1.0 1.5 2.0 2.5 3.0 3.5)
#eta_RangeValues=(-1.0 -0.5 0 0.5 1.0)
#eta_RangeValues=( -3.5 -3.0 -2.5 -2.0 -1.5 -1.0)

#eta_RangeValues=(-3.5 -2.0 -1.0 1.0 2.0 3.5)

#eta_RangeValues=( -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0)
eta_RangeValues=( -3.5 -2.5 -1.0 0.0)

#eta_RangeValues=( -3.45 -3.35 )
#eta_RangeValues=(  -1.95 -1.85 -1.75)

#eta_RangeValues=( -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5)
#eta_RangeValues=( -2.1 -2.05)
#eta_RangeValues=( -1.5 -1.4 -1.3 -1.2 -1.1 -1.0 -0.9)
#eta_RangeValues=( -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5)
#eta_RangeValues=( -2.3 -2.2 -2.1 -2.0 -1.9 -1.8 -1.7)
#eta_RangeValues=( -2.3 -2.2 -2.1 -2.0 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1 -1.0 -0.9 -0.8 -0.7)
#eta_RangeValues=( -2.3 -2.2 -2.1 -2.0 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1)

#eta_RangeValues=( -1.775 -1.725)

#eta_RangeValues=( -3.5 -3.0)
#eta_RangeValues=( -3.0 -2.5)
#eta_RangeValues=( -2.5 -1.0)
#eta_RangeValues=( -1.0 1.0)
#eta_RangeValues=( 1.0 2.5)
#eta_RangeValues=( 2.5 3.0)
#eta_RangeValues=( 3.0 3.5)

mombins=15
mom_RangeValues=( 0.0 30.0)

#0 to use momentum binning, 1 to use pt binning
PtOption=0

#Variable to be plotted on Y-axis
#Yaxis: 0=momentum resolution  1=transverse pointing resolution  2=tracking efficiency  3=projection resolutions (requires further choices below)  4=longitudinal pointing resolution
Yaxis=0
#if Yaxis=3, the following must be specified
#ProjectionSurface:  0=DIRC  1=dRICH Entrance  2=mRICH 3=dRICH Exit 4=Vertex (not from a dummy layer) 
ProjectionSurface=0
#ProjectionValue: 0=Pointing resolution on surface  1: Theta resolution  2: Phi resolution
ProjectionValue=1

#Fit the plot with the functional form specified by the PWG requirement: 0 is off, 1 is on
DrawFits=0
#Draw the PWG requirements if the specifed range and plotting variables match one from table 11.2 in Yellow Report: 0 is off, 1 is on
DrawPWGReq=1


for plotVal in  1; do
	Yaxis=$plotVal
	#ProjectionValue=$plotVal
	
	#for etaindex in {0..15..1}; do
	#for etaindex in {0..13..1}; do
	#for etaindex in {0..11..1}; do
	#for etaindex in {0..6..1}; do
	#for etaindex in {0..5..1}; do
	#for etaindex in {0..4..1}; do
	#for etaindex in {0..3..1}; do
	for etaindex in {0..2..1}; do
	#for etaindex in {0..1..1}; do
	#for etaindex in {0..0..1}; do
		for momindex in {0..0..1}; do
		#for momindex in {0..3..1}; do
	
	#		root -b -q -l './NewFitAnalysis.C('${etabins}','${eta_RangeValues[$etaindex]}','${eta_RangeValues[$etaindex+1]}','${mombins}','${mom_RangeValues[$momindex]}','${mom_RangeValues[$momindex+1]}','${PtOption}','${Yaxis}','${ProjectionSurface}','${ProjectionValue}','${DrawFits}','${DrawPWGReq}')'
			root -b -q -l './NewAnalysis.C('${etabins}','${eta_RangeValues[$etaindex]}','${eta_RangeValues[$etaindex+1]}','${mombins}','${mom_RangeValues[$momindex]}','${mom_RangeValues[$momindex+1]}','${PtOption}','${Yaxis}','${ProjectionSurface}','${ProjectionValue}','${DrawFits}','${DrawPWGReq}')'

		done
	done	
done
