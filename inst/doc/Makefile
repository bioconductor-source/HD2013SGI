all:
#	$(R_HOME)/bin/R CMD Stangle Figures.Rnw
	$(R_HOME)/bin/R CMD Sweave Figures.Rnw
	$(R_HOME)/bin/R CMD pdflatex Figures.tex
	$(R_HOME)/bin/R CMD pdflatex Figures.tex

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
