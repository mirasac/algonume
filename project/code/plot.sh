#!/bin/sh

# Global variables and functions definition.
FILENAME_SCRIPT='plot.sh'
DEF_DIR_EXPORT='../report/figures'

print_info() {
	echo "${FILENAME_SCRIPT} INFO: $1"
}

print_error() {
	echo "${FILENAME_SCRIPT} ERROR: $1"
}

print_usage() {
	cat <<-EOF
		${FILENAME_SCRIPT}
		Create plot in a separate window using a specific gnuplot script.
		
		The first argument is mandatory and its value is used as basename for the gnuplot script.
		If the second argument is specified, the plot is also exported to a file. If the second argument has value 'export', then the exported file has same basename of the gnuplot script, else it has basename equal to the specified value.
		If the third argument is specified, the directory specified as value is used to save the exported file, else the default path '${DEF_DIR_EXPORT}' is used.
	EOF
}

# Get arguments.
basename="$1"
if [ -z "${basename}" ]
then
	print_error 'the first argument is mandatory'
	print_usage
else
	basename_export="$2"
	dir_export="${3:-${DEF_DIR_EXPORT}}"

	# Main plot command.
	command="load '${basename}.gp'"

	# Add export commands.
	if [ -n "${basename_export}" ]
	then
		if [ "${basename_export}" = 'export' ]
		then
			basename_export="${basename}"
		fi
		command=$(
		cat <<-EOF
			${command}
			set terminal pdfcairo
			set output "${dir_export}/${basename_export}.pdf"
			replot
		EOF
		)
	fi

	# Launch gnuplot.
	$(
	gnuplot --persist <<-EOF
		${command}
	EOF
	)
fi
