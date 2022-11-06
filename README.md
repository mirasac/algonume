# algonume
## Introduction
This repository contains source codes of exercises performed during the course of [Numerical Algorithms for Physics](https://fisica-sc.campusnet.unito.it/do/corsi.pl/Show?_id=3881), taken during the first semester of the academic year 2022 for the master's degree in Physics of Complex Systems at UniTo.

Professor Mignone stores slides for each lesson on [his site](http://personalpages.to.infn.it/%7emignone/Numerical_Algorithms/).

### Disclaimer
This repository exists just as an effective way to save the code while the exercises are being done, it is not meant for presenting them in a complete and nice way. Just the information sufficient for the course are included.

## Repository structure
Files and exercises used in each lesson are contained in folders named with the day of the relative lesson in [ISO format](https://en.wikipedia.org/wiki/ISO_8601#Dates) for calendar dates, with dashes.

## Local environment setup
I have a PC with Windows 10 but a Unix environment is mandatory for the course. Since I already have a VM with an Ubuntu image, I use it for the purpose.

I use VS Code as visual text editor, which I installed on the VM.

### Failed attempts
I tried to use WSL2 and install all the request software there, including [gnuplot](http://www.gnuplot.info/), which is used for the graphical part of programs we write. To display the plots created by gnuplot on the WSL2 side, a [X11 server](https://sourceforge.net/projects/xming/) on Windows was needed and I tried to install it following [this guide](https://blog.karatos.in/a?ID=01700-6d257862-8225-4d2a-b4cd-140b2fba8020), but it did not work.
