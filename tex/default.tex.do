pwd >&2
name="${1##*/}"
redo-ifchange ../metrics/$name.metric ../calculus/metric2energy

../calculus/metric2energy ../metrics/$name.metric $3 >&2
