PROJECT = docs
TEX = pdflatex 
BIB = bibtex
BUILDTEX = $(TEX) $(PROJECT).tex

all:
	    $(BUILDTEX)
		 $(BIB) $(PROJECT).aux
	    $(BUILDTEX)
	    evince $(PROJECT).pdf & 

clean-all:
	    rm -f *.dvi *.log *.bak *.aux *.bbl *.blg *.idx *.ps *.eps $(PROJECT).pdf *.toc *.out 

clean:
	    rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out 
