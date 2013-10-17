all:
	$(R_HOME)/bin/R CMD Stangle src/sec02_dataAccess.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec03_figures.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec04_imageProcessing.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec05_convertData.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec06_rawDataProcessing.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec07_exampleInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec08_siRNAcorrelation.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec09_heatmaps.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec10_smallCellNumber.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec11_distributionOfInteractions.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec12_scatterplots.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec13_screenplots.Rnw
	$(R_HOME)/bin/R CMD Stangle src/sec14_table.Rnw
	echo 'Sweave("src/sec02_dataAccess.Rmw")' | R --slave --vanilla
	echo 'Sweave("src/sec03_figures.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec04_imageProcessing.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec05_convertData.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec06_rawDataProcessing.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec07_exampleInteractions.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec08_siRNAcorrelation.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec09_heatmaps.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec10_smallCellNumber.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec11_distributionOfInteractions.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec12_scatterplots.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec13_screenplots.Rnw")' | R --slave --vanilla
	echo 'Sweave("src/sec14_table.Rnw")' | R --slave --vanilla
#	$(R_HOME)/bin/R CMD Sweave src/sec02_dataAccess.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec03_figures.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec04_imageProcessing.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec05_convertData.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec06_rawDataProcessing.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec07_exampleInteractions.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec08_siRNAcorrelation.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec09_heatmaps.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec10_smallCellNumber.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec11_distributionOfInteractions.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec12_scatterplots.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec13_screenplots.Rnw
#	$(R_HOME)/bin/R CMD Sweave src/sec14_table.Rnw
	$(R_HOME)/bin/R CMD Sweave --pdf HD2013SGI.Rnw
#	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex
#	$(R_HOME)/bin/R CMD bibtex HD2013SGI
#	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex
#	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex

clean:
	make texclean
	rm -rf result

texclean:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.bbl
	rm -f *.blg

pdfclean:
	make clean
	rm -f Figures.R
	rm -f sec*.R
	rm -f *.pdf
