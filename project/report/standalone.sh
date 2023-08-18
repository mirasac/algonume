#!/bin/sh

# Global variables and functions definition.
FILENAME_SCRIPT='standalone.sh'
DEF_DIR_FIGURE='./figures'

print_error() {
	echo "${FILENAME_SCRIPT} ERROR: $1"
}

print_usage() {
	cat <<-EOF
		${FILENAME_SCRIPT}
		Create standalone PDF figure starting from its LaTeX source.
		
		The first argument is mandatory and it is the basename of the output PDF figure. The LaTeX source file must be in the same directory and have the same basename with suffix '_input'.
		If the second argument is specified, its value is the directory of the LaTeX source file, else the default value '${DEF_DIR_FIGURE}' is used.
	EOF
}

# Check mandatory arguments.
basename_figure="$1"
if [ -z "${basename_figure}" ]
then
	print_error 'the first argument is mandatory'
	print_usage
	return 1
fi

# Create standalone LaTeX source file.
dir_figure="${2:-${DEF_DIR_FIGURE}}"
filename_tmp="${dir_figure}/${basename_figure}.tmp"
if [ -f "$(pwd)/standalone_packages.tex" ]
then
	preamble="\input{$(pwd)/standalone_packages.tex}"
else
	preamble=$(
	cat <<-EOF
		\usepackage{graphicx}
		\usepackage{xcolor}
	EOF
	)
fi
$(
cat <<-EOF > "${filename_tmp}"
	\documentclass[class=article,border=0mm,10pt]{standalone}
	${preamble}
	\begin{document}
		\input{${basename_figure}_input.tex}
	\end{document}
EOF
)

# Compose standalone LaTeX source file.
pdflatex -jobname="${basename_figure}" -output-directory="${dir_figure}" "${filename_tmp}"

# Remove temporary files.
rm "${filename_tmp}"
rm "${dir_figure}/${basename_figure}.aux"
rm "${dir_figure}/${basename_figure}.log"
