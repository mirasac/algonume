#!/bin/sh

# Global variables and functions definition.
FILENAME_SCRIPT='plot.sh'
DEF_DIR_EXPORT='../report/figures'
BASENAME_LUA_TIKZ='gnuplot-lua-tikz'

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
			set terminal lua tikz latex color gparrows gppoints picenvironment latex createstyle
			set output "${dir_export}/${basename_export}.tex"
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
	
	# Convert TikZ source to PDF figure.
	if [ "${option}" = 'export' ]
	then
		filename_tmp_source="${dir_export}/${basename_export}.tmp"
		# Create LaTeX source.
		$(
		cat <<-EOF > "${filename_tmp_source}"
			\documentclass[tikz, border=0mm]{standalone}
			\usepackage{${BASENAME_LUA_TIKZ}}
			\begin{document}
			\include{${basename_export}.tex}
			\end{document}
		EOF
		)
		# Convert LaTeX source to PDF figure.
		pdflatex -cnf-line='main_memory = 10000000' -cnf-line='main_memory.pdflatex = 10000000' -output-directory="${dir_export}" "${filename_tmp_source}"
		# Delete temporary files.
		rm -f "${filename_tmp_source}*"
		rm "${BASENAME_LUA_TIKZ}.tex"
		rm "t-${BASENAME_LUA_TIKZ}.tex"
	fi
fi
