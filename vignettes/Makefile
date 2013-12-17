all:
	echo 'Stangle(file.path("src","sec02_dataAccess.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec03_figures.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec04_imageProcessing.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec05_convertData.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec06_rawDataProcessing.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec07_exampleInteractions.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec08_siRNAcorrelation.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec09_heatmaps.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec10_smallCellNumber.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec11_distributionOfInteractions.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec12_scatterplots.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec13_screenplots.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle(file.path("src","sec14_table.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Stangle("HD2013SGI.Rnw")' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec02_dataAccess.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec03_figures.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec04_imageProcessing.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec05_convertData.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec06_rawDataProcessing.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec07_exampleInteractions.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec08_siRNAcorrelation.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec09_heatmaps.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec10_smallCellNumber.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec11_distributionOfInteractions.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec12_scatterplots.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec13_screenplots.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave(file.path("src","sec14_table.Rnw"))' | $(R_HOME)/bin/R --slave --vanilla
	echo 'Sweave("HD2013SGI.Rnw")' | $(R_HOME)/bin/R --slave --vanilla
	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex
	$(R_HOME)/bin/R CMD bibtex HD2013SGI
	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex
	$(R_HOME)/bin/R CMD pdflatex HD2013SGI.tex

#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec02_dataAccess.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec03_figures.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec04_imageProcessing.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec05_convertData.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec06_rawDataProcessing.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec07_exampleInteractions.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec08_siRNAcorrelation.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec09_heatmaps.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec10_smallCellNumber.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec11_distributionOfInteractions.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec12_scatterplots.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec13_screenplots.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette(file.path("src", "sec14_table.Rnw"),latex=FALSE,tangle=TRUE,quite=FALSE,clean=FALSE)"
#	$(R_HOME)/bin/R -e "tools::buildVignette("HD2013SGI.Rnw",latex=TRUE,tangle=TRUE,quite=FALSE,clean=FALSE)"
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
