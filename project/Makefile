DIR_CODE=code
DIR_DATA=data
DIR_REPORT=report
MAIN=main
FIGURE=figure

.PHONY: all clean

all:
	make --directory=${DIR_CODE} MAIN=${MAIN}
	(cd ${DIR_CODE}; ./${MAIN}.out)
	(cd ${DIR_DATA}; sh plot.sh ${FIGURE} exportonly)
	make --directory=${DIR_REPORT} MAIN=casari_project all

clean:
	make --directory=${DIR_CODE} clean
	make --directory=${DIR_REPORT} clean
