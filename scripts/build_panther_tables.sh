#!/bin/sh

# Run pdflatex twice, in case table widths change

for X in 1 2; do
    cd results/panther_results/parsed_results/windowed/enrich/
    for tex in $( ls panther_enrichment_outliers_*.tex ); do
	pdflatex $tex
    done
    cd ../../../../..

    cd results/panther_results/parsed_results/windowed/overrep/
    for tex in $( ls panther_overrep_outliers_*.tex ); do
	pdflatex $tex
    done
    cd ../../../../..
    
    cd results/panther_results/parsed_results/roi/enrich/
    for tex in $( ls panther_enrichment_outliers_*.tex ); do
        pdflatex $tex
    done
    cd ../../../../..

done
