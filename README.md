# algonume
## Introduction
This repository contains source codes of exercises performed during the course of [Numerical Algorithms for Physics](https://fisica-sc.campusnet.unito.it/do/corsi.pl/Show?_id=3881), taken during the first semester of the academic year 2022 for the master's degree in Physics of Complex Systems at UniTo.

The regular teacher of this course is Professor Mignone. He stores slides for each lesson on [his site](http://personalpages.to.infn.it/%7emignone/Numerical_Algorithms/).

### Disclaimer
This repository exists just as an effective way to save the code while the exercises are being done, it is not meant for presenting them in a complete and nice way. Just the information sufficient for the course are included.

### Chapters
Course topics are grouped in chapters as shown on the course web page. A chapter can span more lessons, so the following is the agenda:
1. 2022-09-29 - Chapter 1
2. 2022-09-30 - Chapter 1
3. 2022-10-06 - Chapter 2
4. 2022-10-07 - Chapter 2, chapter 3
5. 2022-10-13 - Chapter 3
6. 2022-10-14 - Chapter 3
7. 2022-10-20 - Chapter 5
8. 2022-10-21 - Chapter 5
9. 2022-10-27 - Chapter 5
10. 2022-10-28 - Chapter 6
11. 2022-11-03 - Chapter 4
12. 2022-11-04 - Chapter 4, chapter 7
13. 2022-11-10 - Chapter 7
14. 2022-11-11 - Chapter 7
15. 2022-11-17 - Chapter 8
16. 2022-11-18 - Chapter 8
17. 2022-11-24 - Chapter 9
18. 2022-11-25 - Chapter 9
19. 2022-12-01 - Chapter 10
20. 2022-12-02 - Chapter 10

## Repository structure
Files and exercises used in each lesson are contained in folders named with the day of the relative lesson in [ISO format](https://en.wikipedia.org/wiki/ISO_8601#Dates) for calendar dates, with dashes.

## Local environment setup
I have a PC with Windows 10 but a Unix environment is mandatory for the course. Since I already have a VM with an Ubuntu image, I use it for the purpose.

I use VS Code as visual text editor, which I installed on the VM.

### Failed attempts
I tried to use WSL2 and install all the request software there, including [gnuplot](http://www.gnuplot.info/), which is used for the graphical part of programs we write. To display the plots created by gnuplot on the WSL2 side, a [X11 server](https://sourceforge.net/projects/xming/) on Windows was needed and I tried to install it following [this guide](https://blog.karatos.in/a?ID=01700-6d257862-8225-4d2a-b4cd-140b2fba8020), but it did not work.
