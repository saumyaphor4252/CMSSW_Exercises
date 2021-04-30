#!/bin/sh

for mu in PF Glob Loose Soft Tight Iso TrkIso; do

file=treeTMPCOUNT

# for reco Vtx , use following histograms
#histo1d=hFillNVtx
#histo1n=hFillNVtx${mu} 

histo1d=hFillProbePt
histo1n=hFillProbePt${mu}

if [ "$mu" == "Iso" ]; then
histo1d=hFillProbePtTight
histo1n=hFillProbePt${mu}
#histo1d=hFillNVtxTight
#histo1n=hFillNVtx${mu}
fi

xmin=20
xmax=200
rebin=8
#xmin=-3
#xmax=3
#rebin=4

Xtitle="probe p_{T}"
#Xtitle="# good reconstructed vertices"
#Xtitle="probe #eta"
Ytitle="Efficiency of ${mu} muons"
echo $Ytitle
root -l -b -q  "plotter.C(\"${file}\", \"${histo1d}\",\"${histo1n}\" ,\"${Xtitle}\",\"${Ytitle}\",${xmin},${xmax},${rebin})"

done
