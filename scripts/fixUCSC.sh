## color ucsc tracks by string matching
awk ' { if ( $3 ~ /hs/ ) print $0,"color=0,0,255"; else if ( $3 ~ /pt/ ) print $0,"color=0,255,0"; else print $0,"color=255,0,0;" } ' UCSC.bw.txt | sort > UCSC.bw.fix.txt
