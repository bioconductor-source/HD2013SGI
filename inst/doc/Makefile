all:
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec02_dataAccess.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec03_figures.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec04_imageProcessing.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec05_convertData.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec06_rawDataProcessing.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec07_exampleInteractions.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec08_siRNAcorrelation.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec09_heatmaps.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec10_smallCellNumber.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec11_distributionOfInteractions.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec12_scatterplots.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec13_screenplots.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("src/sec14_table.Rnw",latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)'
	$(R_HOME)/bin/R -e 'tools:buildVignette("HD2013SGI.Rnw",latex=TRUE,tangle=TRUE,quite=FALSE,clean=FALSE)'
#	$(R_HOME)/bin/R CMD Sweave --pdf HD2013SGI.Rnw
	make clean

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
