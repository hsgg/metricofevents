#!/bin/bash
# HSGG 20080721

INFILE="$1"
TEXBASE="${2%.pdf}"

PATH="$PATH:`dirname "$0"`"
echo PATH=$PATH

if [[ "$#" != "2" ]]; then
	echo "Usage: $0 <infile.metric> <outfile.pdf>"
	exit 1
fi

metric2energy "$INFILE" "$TEXBASE".tex || exit $?
pdflatex -halt-on-error -interaction=nonstopmode "$TEXBASE".tex || exit $?
evince "$TEXBASE".pdf || exit $?
