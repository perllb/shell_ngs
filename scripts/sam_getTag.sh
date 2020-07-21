#!/bin/sh

## 
##
## Script to extract values of tag in sam file
## sh sam_getTag.sh -s <samfile> -t <tag to extract> -h <help>
##

unset samfile
unset tag
seq=0

while getopts ":i:t:sh" OPTION; do
    case "$OPTION" in
	i) samfile=$OPTARG ;;
	t) tag=$OPTARG ;;
	s) seq=1 ;;
	h)     
	    echo "Usage: "
	    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"
	    exit 0 ;;
	\?) 
	    echo "Invalid option: $OPTARG" 
	    echo "Usage: "
	    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"
	    exit 1;;
	: ) 
	    echo "Invalid option_ $OPTARG requires an argument" 
	    echo "Usage: "
	    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"
	    exit 1;;
	    
    esac
done

if [ $OPTIND -eq 1 ]
then
    echo "Usage: "
    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"   
    exit 1;
fi

if [ -z "$samfile" ]
then
    echo "> Error: No samfile argument"
    echo "Usage: "
    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"    
    exit 1;
fi

if [ -z "$tag" ]
then
    echo "> Error: No sam TAG argument"
    echo "Usage: "
    echo "   sam_getTag.sh -i <samfile> -t <tag to extract> -s <set if mapped sequence should be included> -h <help>"   
    exit 1;
fi

# Check if sequence should be printed
if [ $seq == 0 ]
then
    
    awk -v tagvar="$tag" '{for (I=1;I<NF;I++) if ($I ~ tagvar) {print $I};}' $samfile | sed 's/'"${tag}"'//g'

else

    awk -v tagvar="$tag" '{for (I=1;I<NF;I++) if ($I ~ tagvar) {print $I "\t" $10 };}' $samfile | sed 's/'"${tag}"'//g'

fi
