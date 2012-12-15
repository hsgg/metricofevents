#!/bin/sh

#mencoder mf://p-*.png -mf w=800:h=600:fps=200 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o movie.avi

mencoder mf://tmpdir/png/p-*.png -mf w=800:h=600:fps=100:type=png -ovc lavc -oac copy -o movie.avi
