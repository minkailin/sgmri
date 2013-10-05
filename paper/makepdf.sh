#!/bin/bash
latex paper.tex
dvips -o paper.ps paper.dvi
ps2pdf -sPAPERSIZE=a4 paper.ps
echo 'done'
