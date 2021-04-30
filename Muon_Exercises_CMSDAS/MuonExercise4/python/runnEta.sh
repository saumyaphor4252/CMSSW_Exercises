#!/bin/sh

for mu in PF Glob Loose Soft Tight Iso; do

file=treeTMPCOUNT

histo1d=hFillProbeEta
histo1n=hFillProbeEta${mu} 
if [ "$mu" == "Iso" ]; then
histo1d=hFillProbeEtaTight
histo1n=hFillProbeEta${mu}

fi

xmin=-2.4
xmax=2.4
rebin=4

Xtitle="probe #eta"
Ytitle="Efficiency of ${mu} muons"
echo $Ytitle
root -l -b -q  "plotter.C(\"${file}\", \"${histo1d}\",\"${histo1n}\" ,\"${Xtitle}\",\"${Ytitle}\",${xmin},${xmax},${rebin})"

done
