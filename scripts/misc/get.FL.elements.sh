## get nearly full length elements from hg38.fa.out
## Allow a truncation from start site of 2% of full length of consensus, and 0.1% of end.

cat hg38.fa.out.sortedPercDiv.ERE | tr -d '()' | awk ' {
  fl=$12+$13+$14;

  if ( $9 == "+" ) {
    if ( $12/fl < 0.02 && $14/fl < 0.001 ) {
      print $0;
    }
  } else {
    if ( $14/fl < 0.02 && $12/fl < 0.001 ) {
      print $0;
    }
  }
} ' > hg38.fa.out.sortedPercDiv.ERE.fl_0.02start_0.001end
