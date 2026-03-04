
# input alphaSat suprachromosomal (SF) file and the name of the assembly 
SF_file=$1
ASM_Name=$2

# split into each superfamily 
# dealing with the active SF1 and the subtypes 
grep -w -e "S3-1" -e "S4-1" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-1.merged.bed

grep -w -e "S3-2" -e "S4-2" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-2.merged.bed

grep -w -e "S3-3" -e "S4-3" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-3.merged.bed

grep -w -e "S3-4" -e "S4-4" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-4.merged.bed

grep -w -e "S3-5" -e "S4-5" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-5.merged.bed

grep -w -e "S3-6" -e "S4-6" $SF_file > active_dimer.bed
bedtools merge -d 350 -c 4 -o distinct -i   active_dimer.bed | awk '($3-$2) >= 350'  > active-6.merged.bed



#inactive dimeric
grep -w -e S3b -e S4b $SF_file >  S3bS4b.bed
bedtools merge -d 350 -i  S3bS4b.bed | awk '($3-$2) >= 700' >  S3bS4b.merged.bed


grep -w -e S3c -e S4c $SF_file >  S3cS4c.bed
bedtools merge -d 350 -i  S3cS4c.bed | awk '($3-$2) >= 700'  > S3cS4c.merged.bed

grep -w -e S3d -e S4d $SF_file >  S3dS4d.bed
bedtools merge -d 350 -i  S3dS4d.bed | awk '($3-$2) >= 700'  > S3dS4d.merged.bed

#dHOR

grep -w -e a -e b -e c -e d -e e -e f -e g -e h -e i -e j -e k $SF_file > dhor.bed
bedtools merge -d 350 -i  dhor.bed | awk '($3-$2) >= 700'  > dhor.merged.bed


rm active_dimer.bed
# pave over lines and merge them into a single file for the possibly active SFs 
for f in active*merged.bed ; do
    # I only want the dimeric so I'm pulling those out 
    name=$(grep "," $f | awk '{print $4}' | sort | uniq) 
    bedtools sort -i $f | bedtools merge -d 9500 -c 4 -o distinct -i stdin | awk '($3-$2) >= 900' > tmp.active.bed
    awk 'BEGIN{OFS="\t"} {
        split($4,a,","); 
        delete seen; 
        out="";
        for(i in a){
            if(!(a[i] in seen)){
                seen[a[i]]=1; 
                out=(out==""?a[i]:out","a[i]);
            }
        }
        $4=out; 
        print
    }' tmp.active.bed |  awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, "0", ".", $2, $3, "255,102,0"}' | bedtools sort -i stdin >> active_dimer.bed
done


# extract a list of chromosome 

chromosomes=( $(cut -f1 active_dimer.bed | uniq) )

#incase of rerun 
rm dimer.bed

# for each chromosome determine the largest array of active dimers and mark as active 
for chrom in ${chromosomes[@]}; do 
    grep $chrom active_dimer.bed > tmp.bed \
    && awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, ($3-$2)}' tmp.bed | sort -nk 10,10 -r > tmp.sorted.bed \
    && head -n 1 tmp.sorted.bed | awk '{gsub("mon/", "");print}' | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "active(NSF1-"$4")", $5, $6, $7, $8, "153,0,0"}' >> dimer.bed \
    && tail -n+2 tmp.sorted.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, "dimeric(NSF1-"$4")", $5, $6, $7, $8, $9}' >> dimer.bed
done

# now add in the other dimeric regions 
rm dimer.inactive.bed
bedtools merge -d 1000 -i S3bS4b.merged.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dimeric(NSF2)", "0", ".", $2, $3, "255,102,0"} '>> dimer.inactive.bed
bedtools merge -d 1000 -i S3cS4c.merged.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dimeric(NSF3)", "0", ".", $2, $3, "255,102,0"}'  >> dimer.inactive.bed
bedtools merge -d 1000 -i S3dS4d.merged.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dimeric(NSF4)", "0", ".", $2, $3, "255,102,0"}' >> dimer.inactive.bed

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dhor(NSF6)", "0", ".", $2, $3, "255,146,0"}' dhor.merged.bed >> dimer.inactive.bed

cat dimer.bed dimer.inactive.bed | bedtools sort -i stdin > dimer.final.bed


# identify the monomeric regions 
bedtools merge -d 500 -i $SF_file >  mon_merged.bed

# in case of rerun 
rm mon.final.bed
# subtract and pave over LINES 
bedtools subtract -A -a  mon_merged.bed -b  dimer.final.bed | awk '($3-$2) >= 350' | bedtools merge -d 6500  \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "mon", "0", ".", $2, $3, "255,204,153"}' >>  mon.final.bed

# format for UCSC genome browser and censat workflow 
echo 'track name="Simple_'$ASM_Name'_Alpha_Summary" visibility=2 itemRgb="On"' > ${ASM_Name}_alphaSummary.sorted.bed
cat mon.final.bed dimer.final.bed | bedtools sort -i stdin | uniq  >>  ${ASM_Name}_alphaSummary.sorted.bed

rm tmp*