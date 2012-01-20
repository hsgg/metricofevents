metrics=""
for i in ../metrics/*.metric; do
	name="${i##../metrics/}"
	name="${name%.metric}"
	case "$name" in
		"schwarzschild-cyl") ;;
		"robertson-walker-schwarzschild") ;;
		"kerr-newman") ;;
		*) metrics="$metrics ${name}.pdf" ;;
	esac
done
echo $metrics >&2
redo-ifchange $metrics
