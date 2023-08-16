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
		Create plot using a specific gnuplot script.
		
		The first argument is mandatory and its value is used as basename for the gnuplot script.
		If the second argument is specified with value 'export', then the plot is also exported to a file with same basename of the gnuplot script. The default directory for the exported file is '${DEF_DIR_EXPORT}'.
		If the third argument is specified, its value is used as directory for the exported file.
	EOF
}

# Script body.
basename="$1"
if [ -z "${basename}" ]
then
	print_error 'the first argument is mandatory'
	print_usage
else
	option="$2"
	
	# Main plot command.
	command="load '${basename}.gp'"
	
	# Add export commands.
	if [ "${option}" = 'export' ]
	then
		dir_export="${3:-${DEF_DIR_EXPORT}}"
		basename_export="${basename}"
		command=$(
		cat <<-EOF
			${command}
			set terminal cairolatex pdf input colourtext noheader color size 14.85cm,10.5cm
			set output '${dir_export}/${basename_export}.tex'
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
