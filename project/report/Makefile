MAIN=casari_project

.PHONY: all clean

${MAIN}.pdf: ${MAIN}.tex ${MAIN}.bbl
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex
	touch ${MAIN}.bbl
	touch ${MAIN}.pdf

${MAIN}.bbl: ${MAIN}.bcf bibliography.bib
	biber ${MAIN}

${MAIN}.bcf: ${MAIN}.tex
	pdflatex ${MAIN}.tex

all: ${MAIN}.tex bibliography.bib
	pdflatex ${MAIN}.tex
	biber ${MAIN}
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex

clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.bcf
	rm -f *.blg
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.xml
