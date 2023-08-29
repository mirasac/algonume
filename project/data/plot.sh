#!/bin/sh

# Global variables and functions definition.
FILENAME_SCRIPT='plot.sh'
DEF_DIR_EXPORT='../report/figures'
DEF_NULL='NULL'

print_info() {
	echo "${FILENAME_SCRIPT} INFO: $1"
}

print_error() {
	echo "${FILENAME_SCRIPT} ERROR: $1"
}

print_usage() {
	cat <<-EOF
		${FILENAME_SCRIPT}
		Create plot using a specific gnuplot script.
		
		The first argument is mandatory and its value is used as basename for the gnuplot script.
		If the second argument is specified, if its value is 'export' then the plot is also exported to a file with same basename of the gnuplot script, if its value is 'onlyexport' then the plot is not displayed in a separate window and only the export is performed (e.g. useful in systems which lack of GUI support), if its value is 'exportpdf' then the same exporting of value 'export' is performed but with terminal pdfcairo, if its value is 'onlypdf' then the same exporting of value 'onlyexport' is performed but with terminal pdfcairo.
		If the third argument is specified, its value is used as directory for the exported file, else the default directory is '${DEF_DIR_EXPORT}'.
	EOF
}

# Check mandatory arguments.
basename="$1"
if [ -z "${basename}" ]
then
	print_error 'the first argument is mandatory'
	print_usage
	return 1
fi

option="${2:-${DEF_NULL}}"

# Main plot command.
command="load '${basename}.gp'"

# Add export commands.
case "${option}" in
export*|only*)
	dir_export="${3:-${DEF_DIR_EXPORT}}"
	basename_export="${basename}"
	if [ "${option}" = 'export' ] || [ "${option}" = 'onlyexport' ]
	then
		command_export=$(
		cat <<-EOF
			set terminal cairolatex pdf input colourtext noheader color size 14.85cm,10.5cm
			set output '${basename_export}_input.tex'
		EOF
		)
	elif [ "${option}" = 'exportpdf' ] || [ "${option}" = 'onlypdf' ]
	then
		command_export=$(
		cat <<-EOF
			set terminal pdfcairo noenhanced color font ",10" size 14.85cm,10.5cm
			set output '${basename_export}.pdf'
		EOF
		)
	fi
	command=$(
	cat <<-EOF
		${command}
		${command_export}
		replot
	EOF
	)
	case "${option}" in
	only*)
		command=$(
		cat <<-EOF
			set terminal unknown
			${command}
		EOF
		)
	;;
	esac
;;
"${DEF_NULL}")
	:
;;
*)
	print_info "invalid value for the second argument, it is ignored"
;;
esac

# Launch gnuplot.
$(
gnuplot --persist <<-EOF
	${command}
EOF
)

# Move exported files to the specified directory.
if [ "${option}" = 'export' ] || [ "${option}" = 'onlyexport' ]
then
	mv "${basename_export}_input.tex" "${dir_export}/${basename_export}_input.tex"
	mv "${basename_export}_input.pdf" "${dir_export}/${basename_export}_input.pdf"
elif [ "${option}" = 'onlypdf' ]
then
	mv "${basename_export}.pdf" "${dir_export}/${basename_export}.pdf"
fi
