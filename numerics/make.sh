#!/bin/sh

GINACDIR=/home/gebhardt/public_bin/ginac-1.5.6-bin
CLNDIR=/home/gebhardt/public_bin/cln-1.3.1-bin

CFLAGS="-I$GINACDIR/include -I$CLNDIR/include"
LDFLAGS="-static $GINACDIR/lib/libginac.a $CLNDIR/lib/libcln.a -lgmp"


make CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS"
