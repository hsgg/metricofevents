#!/bin/sh

if test ! -e tmpdir; then
	tmp=`mktemp -d --tmpdir binary.XXXXXXXXXX`
	ln -s "$tmp" tmpdir
fi
tmp=tmpdir
ln -sf "`pwd`" "$tmp/srcdir"

f=$tmp/plot.gplt

mkdir -p $tmp/png

# need noly first two columns
echo Creating $tmp/plot0-first.dat...
cut -f1-2 -d' ' plot0.dat >$tmp/plot0-first.dat
echo Creating $tmp/plot1-first.dat...
cut -f1-2 -d' ' plot1.dat >$tmp/plot1-first.dat

echo Initializing iteration...
max=`wc -l <$tmp/plot0-first.dat`
for i in `seq -w 4 500 $max`
do
echo $i/$max

echo "#set title \"Metric\"
set polar
set grid polar
set nokey

set size square
set xrange[-6e11:6e11]
set yrange[-6e11:6e11]

set terminal png size 800,600

set output \"$tmp/png/p-$i.png\"

plot \"sphere.dat\", \"sphere2.dat\", \"$tmp/plot0-tmp.dat\" w l lw 2, \"$tmp/plot1-tmp.dat\" w l lw 2" >$f

head -n $i $tmp/plot0-first.dat >$tmp/plot0-tmp.dat
head -n $i $tmp/plot1-first.dat >$tmp/plot1-tmp.dat

gnuplot $f
done
