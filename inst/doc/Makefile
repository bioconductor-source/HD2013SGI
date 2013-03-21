all:
	$(R_HOME)/bin/R CMD Stangle src/sec02_dataAccess.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec03_figures.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec04_convertData.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec05_rawDataProcessing.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec06_exampleInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec07_siRNAcorrelation.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec08_heatmaps.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec09_smallCellNumber.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec10_distributionOfInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec11_scatterplots.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec12_screenplots.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec13_table.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec02_dataAccess.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec03_figures.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec04_convertData.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec05_rawDataProcessing.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec06_exampleInteractions.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec07_siRNAcorrelation.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec08_heatmaps.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec09_smallCellNumber.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec10_distributionOfInteractions.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec11_scatterplots.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec12_screenplots.Rnw
	$(R_HOME)/bin/R CMD Sweave src/sec13_table.Rnw
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
