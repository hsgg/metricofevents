mkdir -p $2
redo-ifchange $2/$2.tex

( cd $2; latexmk -pdf $2.tex >&2 )

cp $2/$1 $3
