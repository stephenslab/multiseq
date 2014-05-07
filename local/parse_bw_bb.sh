chrom_file=$1
if [ -z $chrom_file ]; then
    echo "error: provide file with chormosome sizes"
    exit 1
fi
for i in `find . -name "*.bw"`; do
    echo $i
    l=`bigWigInfo -chroms $i | grep chr | wc -l`
    echo $l
    if [ "$l" == 1 ]; then
        bigWigToWig $i $i".wig"
        cat $i".wig" | sed 's/chrom=/chrom=chr/g' > $i".cp.wig"
        wigToBigWig $i".cp.wig" $chrom_file $i".cp.bw"
        mv $i".cp.bw" $i
        rm $i".wig"
        rm $i".cp.wig"
    fi
done
for i in `find . -name *.bb`; do
    echo $i
    l=`bigBedInfo -chroms $i | grep chr | wc -l`
    echo $l
    if [ "$l" == 1 ]; then
        bigBedToBed $i $i".bed"
        cat $i".bed" | awk '{OFS="\t"}{print "chr"$0}' > $i".cp.bed"
        bedToBigBed $i".cp.bed" $chrom_file $i".cp.bb"
        mv $i".cp.bb" $i
        rm $i".bed"
        rm $i".cp.bed"
    fi
done

