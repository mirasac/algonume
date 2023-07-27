MAIN=main

all: main bibliography
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex

main: ${MAIN}.tex
	pdflatex ${MAIN}.tex

bibliography: bibliography.bib
	biber ${MAIN}

clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.bcf
	rm -f *.blg
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.xml
