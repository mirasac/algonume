MAIN=casari_project
DIR_FIGURES=./figures
FIGURES=${DIR_FIGURES}/spectral_irradiance.pdf ${DIR_FIGURES}/spectral_irradiance_logy.pdf

.PHONY: all all_figures clean

${MAIN}.pdf: ${MAIN}.tex ${MAIN}.bbl all_figures
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex
	touch ${MAIN}.bbl
	touch ${MAIN}.pdf

${MAIN}.bbl: ${MAIN}.bcf bibliography.bib
	biber ${MAIN}

${MAIN}.bcf: ${MAIN}.tex
	pdflatex ${MAIN}.tex

all: ${MAIN}.tex bibliography.bib all_figures
	pdflatex ${MAIN}.tex
	biber ${MAIN}
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex

all_figures: ${FIGURES}

${DIR_FIGURES}/spectral_irradiance.pdf: ${DIR_FIGURES}/spectral_irradiance_input.tex ${DIR_FIGURES}/spectral_irradiance_input.pdf
	sh standalone.sh spectral_irradiance ${DIR_FIGURES}

${DIR_FIGURES}/spectral_irradiance_logy.pdf: ${DIR_FIGURES}/spectral_irradiance_logy_input.tex ${DIR_FIGURES}/spectral_irradiance_logy_input.pdf
	sh standalone.sh spectral_irradiance_logy ${DIR_FIGURES}

clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.bcf
	rm -f *.blg
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.xml
