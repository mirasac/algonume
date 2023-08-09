MAIN=casari_project

all: ${MAIN}.pdf ${MAIN}.bbl

${MAIN}.pdf: ${MAIN}.tex
	pdflatex ${MAIN}.tex

${MAIN}.bbl: bibliography.bib ${MAIN}.bcf
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
