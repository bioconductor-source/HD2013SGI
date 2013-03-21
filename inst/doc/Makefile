all:
	cd src
	$(R_HOME)/bin/R CMD Stangle da01_dataAccess.Rnw
	$(R_HOME)/bin/R CMD Stangle da02_figures.Rnw
	$(R_HOME)/bin/R CMD Stangle da11_convertData.Rnw
	$(R_HOME)/bin/R CMD Stangle da12_rawDataProcessing.Rnw
	$(R_HOME)/bin/R CMD Stangle da13_exampleInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle da14_siRNAcorrelation.Rnw
	$(R_HOME)/bin/R CMD Stangle da15_heatmaps.Rnw
	$(R_HOME)/bin/R CMD Stangle da21_smallCellNumber.Rnw
	$(R_HOME)/bin/R CMD Stangle da22_distributionOfInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle da23_scatterplots.Rnw
	$(R_HOME)/bin/R CMD Stangle da24_screenplots.Rnw
	$(R_HOME)/bin/R CMD Stangle da25_table.Rnw
	cd ..
	$(R_HOME)/bin/R CMD Sweave src/da01_dataAccess.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da02_figures.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da11_convertData.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da12_rawDataProcessing.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da13_exampleInteractions.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da14_siRNAcorrelation.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da15_heatmaps.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da21_smallCellNumber.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da22_distributionOfInteractions.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da23_scatterplots.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da24_screenplots.Rnw
	$(R_HOME)/bin/R CMD Sweave src/da25_table.Rnw
	$(R_HOME)/bin/R CMD Sweave HD2013SGI.Rnw
	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex
	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex

clean:
	make texclean
	rm -rf result

texclean:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f *.out
	rm -f *.toc

pdfclean:
	make clean
	rm -f Figures.R
	rm -f *.pdf
