DIR_CODE=code
DIR_DATA=data
DIR_REPORT=report
MAIN=main
FIGURE=

.PHONY: single_figure clean

single_figure:
	@test -n "${FIGURE}" || (echo 'Variable FIGURE must be set to run this rule' && exit 1)
	$(MAKE) --directory=${DIR_CODE} MAIN=${MAIN}
	cd ${DIR_CODE} && ./${MAIN}.out
	cd ${DIR_DATA} && sh plot.sh ${FIGURE} exportonly
	$(MAKE) --directory=${DIR_REPORT} MAIN=casari_project all

clean:
	$(MAKE) --directory=${DIR_CODE} clean
	$(MAKE) --directory=${DIR_REPORT} clean
