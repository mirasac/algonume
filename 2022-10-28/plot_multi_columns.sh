#!/bin/sh

# Default value for the data file.
DEF_FILENAME='plotfile.dat'

# Get name of data file.
filename=${1:-${DEF_FILENAME}}

# Get the number of columns from the header of the data file.
read -r header < ${filename}
set -- ${header}
shift 1
header=$*

# Write command to launch in gnuplot.
#plot 'plotfile.dat' using 0:5, '' using 0:2, '' using 0:3, '' using 0:4
command=''
i=2
for item in ${header}
do
	command="${command}${command:+,} '${filename}' using 1:${i}"
	filename=''
	i=$(( i + 1 ))
done
command="plot${command}"

echo ${command}

# Launch command.
$(
gnuplot --persist <<-EOF
	set log x
	set log y
	${command}
EOF
)

