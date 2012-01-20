mkdir -p $1
redo-ifchange $1/$1.tex

( cd $1; latexmk -pdf $1.tex >&2 )

cp $1/$1.pdf $3
