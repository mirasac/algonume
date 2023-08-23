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
		If the second argument is specified, if its value is 'export' then the plot is also exported to a file with same basename of the gnuplot script, if its value is 'exportonly', the plot is not displayed in a separate window and only the export is performed (e.g. useful in systems which lack of GUI support).
		If the third argument is specified, its value is used as directory for the exported file, else the default directory is '${DEF_DIR_EXPORT}'.
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
	if [ "${option}" = 'export' ] || [ "${option}" = 'exportonly' ]
	then
		dir_export="${3:-${DEF_DIR_EXPORT}}"
		basename_export="${basename}"
		command_export=$(
		cat <<-EOF
			set terminal cairolatex pdf input colourtext noheader color size 14.85cm,10.5cm
			set output '${basename_export}_input.tex'
		EOF
		)
		command=$(
		cat <<-EOF
			${command}
			${command_export}
			replot
		EOF
		)
		if [ "${option}" = 'exportonly' ]
		then
			command=$(
			cat <<-EOF
				set terminal unknown
				${command}
			EOF
			)
		fi
	fi
	
	# Launch gnuplot.
	$(
	gnuplot --persist <<-EOF
		${command}
	EOF
	)
	
	# Move exported files to the specified directory.
	if [ "${option}" = 'export' ] || [ "${option}" = 'exportonly' ]
	then
		mv "${basename_export}_input.tex" "${dir_export}/${basename_export}_input.tex"
		mv "${basename_export}_input.pdf" "${dir_export}/${basename_export}_input.pdf"
	fi
fi
