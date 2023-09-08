MAIN=casari_project
CLASS_OPTIONS=draft
MAIN_CLASS_OPTIONS="\providecommand\classoptions{${CLASS_OPTIONS}}\input{${MAIN}.tex}"
DIR_FIGURES=./figures
FIGURES=${DIR_FIGURES}/spectral_irradiance.pdf ${DIR_FIGURES}/spectral_irradiance_logy.pdf ${DIR_FIGURES}/spectral_irradiance_bands.pdf ${DIR_FIGURES}/spectral_irradiance_wavelength.pdf ${DIR_FIGURES}/shortwave_radiation.pdf ${DIR_FIGURES}/temperature_radiative_equilibrium.pdf ${DIR_FIGURES}/irradiance_radiative_equilibrium.pdf

.PHONY: all class_options all_figures clean

${MAIN}.pdf: ${MAIN}.tex ${MAIN}.bbl ${FIGURES}
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex
	touch --no-create ${MAIN}.pdf ${MAIN}.bbl ${FIGURES}

${MAIN}.bbl: ${MAIN}.bcf bibliography.bib
	biber ${MAIN}

${MAIN}.bcf: ${MAIN}.tex
	pdflatex ${MAIN}.tex

all: ${MAIN}.tex bibliography.bib all_figures
	pdflatex ${MAIN}.tex
	biber ${MAIN}
	pdflatex ${MAIN}.tex
	pdflatex ${MAIN}.tex

class_options: ${MAIN}.tex bibliography.bib
	pdflatex ${MAIN_CLASS_OPTIONS}
	biber ${MAIN}
	pdflatex ${MAIN_CLASS_OPTIONS}
	pdflatex ${MAIN_CLASS_OPTIONS}

all_figures: ${FIGURES}

${DIR_FIGURES}/spectral_irradiance.pdf: ${DIR_FIGURES}/spectral_irradiance_input.tex ${DIR_FIGURES}/spectral_irradiance_input.pdf
	sh standalone.sh spectral_irradiance ${DIR_FIGURES}

${DIR_FIGURES}/spectral_irradiance_logy.pdf: ${DIR_FIGURES}/spectral_irradiance_logy_input.tex ${DIR_FIGURES}/spectral_irradiance_logy_input.pdf
	sh standalone.sh spectral_irradiance_logy ${DIR_FIGURES}

${DIR_FIGURES}/spectral_irradiance_bands.pdf: ${DIR_FIGURES}/spectral_irradiance_bands_input.tex ${DIR_FIGURES}/spectral_irradiance_bands_input.pdf
	sh standalone.sh spectral_irradiance_bands ${DIR_FIGURES}

${DIR_FIGURES}/spectral_irradiance_wavelength.pdf: ${DIR_FIGURES}/spectral_irradiance_wavelength_input.tex ${DIR_FIGURES}/spectral_irradiance_wavelength_input.pdf
	sh standalone.sh spectral_irradiance_wavelength ${DIR_FIGURES}

${DIR_FIGURES}/shortwave_radiation.pdf: ${DIR_FIGURES}/shortwave_radiation.tex
	pdflatex -jobname=shortwave_radiation -output-directory=${DIR_FIGURES} ${DIR_FIGURES}/shortwave_radiation
	rm ${DIR_FIGURES}/shortwave_radiation.aux
	rm ${DIR_FIGURES}/shortwave_radiation.log

${DIR_FIGURES}/temperature_radiative_equilibrium.pdf: ${DIR_FIGURES}/temperature_radiative_equilibrium_input.tex ${DIR_FIGURES}/temperature_radiative_equilibrium_input.pdf
	sh standalone.sh temperature_radiative_equilibrium ${DIR_FIGURES}

${DIR_FIGURES}/irradiance_radiative_equilibrium.pdf: ${DIR_FIGURES}/irradiance_radiative_equilibrium_input.tex ${DIR_FIGURES}/irradiance_radiative_equilibrium_input.pdf
	sh standalone.sh irradiance_radiative_equilibrium ${DIR_FIGURES}

clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.bcf
	rm -f *.blg
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.xml

