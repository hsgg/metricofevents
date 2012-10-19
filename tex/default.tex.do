name="${2##*/}"
redo-ifchange ../metrics/$name.metric

../calculus/metric2energy ../metrics/$name.metric $3 >&2
