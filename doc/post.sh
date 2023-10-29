for i in force.*.dat
do  set -- `basename $i .dat | sed 's/\./ /g'`
    lv=$2
    re=$3
    case "$re" in
	'') re=1000
    esac
    Cd=`awk '$2 > 60 { s += $3; n++} END {if (n > 0) printf "%.16e\n", s / n / 75}' "$i"`
    case "$Cd" in
	'') ;;
	*) echo $lv $re $Cd
    esac
done
