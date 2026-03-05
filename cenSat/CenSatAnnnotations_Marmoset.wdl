version 1.0

workflow cenSatAnnotation{
    input {
         File RMOut
         File aSatBed
         File aSatStrand
         File rDNABed
         File gapBed
         File deNovoSat
         String fName=basename(sub(sub(sub(sub(RMOut, "\\.bed$", ""), "\\_rm$", ""), "\\.fa$", ""), "\\.fasta$", ""))
        
         Int threadCount = 20
         Int preemptible = 1
         Int diskSize = 32
         Int memSizeGB = 32
     }
    

    call createAnnotations {
        input:
            RMOut=RMOut,
            aSatBed=aSatBed,
            aSatStrand=aSatStrand,
            rDNABed=rDNABed,
            gapBed=gapBed,
            deNovoSat=deNovoSat,
            fName=fName,

            preemptible=preemptible,
            threadCount=threadCount,
            diskSize=diskSize,
            memSizeGB=memSizeGB
    }


    output {
        File cenSatAnnotations = createAnnotations.cenSatAnnotations
        File cenSatStrand = createAnnotations.cenSatStrand
        File centromeres = createAnnotations.centromeres
    }
    

    parameter_meta {
         RMOut: "RepeatMasker annotation converted to bed file by python script https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py"
         aSatBed: "Bed file containing alpha Satellite annotations from HMMR"
         rDNABed: "Bed file with rDNA annotations"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

task createAnnotations {
    input {
        File RMOut
        File aSatBed
        File aSatStrand
        File rDNABed
        File gapBed
        File deNovoSat
        String fName

        Int memSizeGB
        Int preemptible
        Int threadCount
        Int diskSize
    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # HSAT1A - SAR - 5kb merge color code 0,222,96
        grep SAR ~{RMOut} > HSAT1A.bed || true
        if [ -s HSAT1A.bed ]; then
            bedtools merge -s -c 6 -o distinct -s -d 50 -i HSAT1A.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "hsat1A", "0", $4, $2, $3, "."}' > strandInfo.bed  # store for strand information later 
            bedtools merge -d 1000 -i HSAT1A.bed > HSAT1A.merged.bed
            sed 's/$/\tHSat1A\t0\t.\t.\t.\t0,222,96/' HSAT1A.merged.bed > HSAT1A.merged.named.bed
            awk '$7=$2' OFS='\t' HSAT1A.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1A.diff.bed
        fi

        # HSAT1B - HSATI - 5kb merge color code 27,153,139
        grep -w "HSATI" ~{RMOut} > HSAT1B.bed || true
        if [ -s HSAT1B.bed ]; then
            bedtools sort -i HSAT1B.bed | bedtools merge -d 5000 -i stdin > HSAT1B.merged.bed
            bedtools sort -i HSAT1B.bed | bedtools merge -s -c 6 -o distinct -s -d 50 -i stdin | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "hsat1B", "0", $4, $2, $3, "."}' >> strandInfo.bed  # store for strand information later 
            sed 's/$/\tHSat1B\t0\t.\t.\t.\t27,153,139/' HSAT1B.merged.bed > HSAT1B.merged.named.bed
            awk '$7=$2' OFS='\t' HSAT1B.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1B.part.bed
        fi

        # BetaSats - BSAT, LSAU, BSR 250,153,255
        grep -e BSAT -e LSAU -e BSR ~{RMOut} > BSAT.bed || true
        if [ -s BSAT.bed ]; then
            bedtools sort -i BSAT.bed > BSAT.sorted.bed
            bedtools merge -s -c 4,6 -o distinct -s -d 50 -i BSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed # store for strand information later 
            bedtools merge -d 5000 -c 4 -o distinct -i BSAT.sorted.bed > BSAT.merged.bed 
            sed 's/$/\t0\t.\t.\t.\t250,153,255/' BSAT.merged.bed > BSAT.merged.named.bed
            awk '$7=$2' OFS='\t' BSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="bSat("$4")"' OFS='\t' > BSAT.part.bed
        fi


        # P-Censat - CER, SATR, SST1, ACRO - many more - 1 kb merge more fine tuned 0,204,204
        grep -e CER -e SATR -e SST1 -e ACRO -e Hsap -e rnd -e HSAT5 -e 5SRNA -e TAF11 -e HSAT4  ~{RMOut} > cenSAT.bed || true
        grep -v BSAT cenSAT.bed > cenSAT.filtered.bed || true



        # iterate through unique elements and merge 
        touch censat.merged.bed
        satellites=( $(awk '{print $4}' cenSAT.filtered.bed | sort | uniq) )
        echo $satellites
        bedtools sort -i cenSAT.filtered.bed > cenSAT.sorted.bed
        
        for sat in ${satellites[@]}; do 
            echo $sat
            grep -w $sat cenSAT.filtered.bed > tmp.bed || true 
            if [ -s tmp.bed ]; then
                bedtools sort -i tmp.bed > tmp.sorted.bed 
                bedtools merge -c 4 -o distinct -d 10000 -i tmp.sorted.bed >> censat.merged.bed
            fi
        done


        if [ -s censat.merged.bed ]; then
            bedtools merge -s -c 4,6 -o distinct -s -d 150 -i cenSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed || true # store for strand information later 
            bedtools sort -i censat.merged.bed > censat.sorted.bed
            bedtools merge -d 6000 -c 4 -o distinct -i cenSAT.sorted.bed > cenSAT.merged.bed || true
            sed 's/$/\t0\t.\t.\t.\t0,204,204/' cenSAT.merged.bed > cenSAT.merged.named.bed
            awk '$7=$2' OFS='\t' cenSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="cenSat("$4")"' OFS='\t' | bedtools sort -i stdin > cenSAT.part.bed
        fi 

        # add the new satellite elements 

        cat ~{deNovoSat} >> deNovo.part.bed
        bedtools merge -s -c 4,6 -o distinct -s -d 50 -i ~{deNovoSat} | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed

        # GammaSats - GSAT, TAR1 
        grep -e GSAT -e TAR1 ~{RMOut} > GSAT.bed || true
        if [ -s GSAT.bed ]; then
            bedtools sort -i GSAT.bed > GSAT.sorted.bed
            bedtools merge -s -c 4,6 -o distinct -s -d 50 -i GSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed || true # store for strand information later 
            bedtools merge -d 2000 -c 4 -o distinct -i GSAT.sorted.bed | bedtools sort -i stdin > GSAT.merged.bed
            # there are some compound arrays in primates with gamma and other censat - so to resolve we subtract the censat from the gamma 
            bedtools subtract -a GSAT.merged.bed -b cenSAT.part.bed > tmp.bed && mv tmp.bed GSAT.merged.bed
            sed 's/$/\t0\t.\t.\t.\t172,51,199/' GSAT.merged.bed > GSAT.merged.named.bed
            awk '$7=$2' OFS='\t' GSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="gSat("$4")"' OFS='\t' > GSAT.part.bed
        fi

        # pull the simple repeats - and subtracting with the repeatModeler models since there is overlap 
        # pentamer GGCA + strand 
        grep -e "(GGCAA)" -e "(GCAAG)" -e "(CAAGG)" -e "(AAGGC)" -e "(AGGCA)" ~{RMOut} | bedtools sort | bedtools merge -d 3000 | bedtools subtract -a stdin -b ~{deNovoSat} | awk '{print $1, $2, $3, "pentameric(GGCAA)", "0", "+", $2, $3, "51,51,102"}' OFS='\t' > simple.repeats.part.bed
        #pentamer GGCA - strand 
        grep -e "(TTGCC)" -e "(TGCCT)" -e "(GCCTT)" -e "(CCTTG)" -e "(CTTGC)"  ~{RMOut} | bedtools sort | bedtools merge -d 3000 | bedtools subtract -a stdin -b ~{deNovoSat} | awk '{print $1, $2, $3, "pentameric(GGCAA)", "0", "-", $2, $3, "51,51,102"}' OFS='\t' >> simple.repeats.part.bed


        # add strand info from simple repeats 
        cat simple.repeats.part.bed >> strandInfo.bed

        


                # AlphaSat - resolve overlaps 
        # merge any overlaps smaller than 400 bp to upstream annotation - 2+ alpha monomers 
        bedtools closest -D a -iu -t last -a ~{aSatBed} -b ~{aSatBed} | awk ' BEGIN {OFS="\t"} {if ($19<1&&(($3-$11)<400)&&($2!=$11||$3!=$12)) {$3=$8=($11-1)}} {print $1,$2,$3,$4,$5,$6,$2,$3,$9}' > smallMerged.bed

        # merge all alpha annotations for coverage calculation
        bedtools merge -d 2000000 -i ~{aSatBed} > wholealpha.bed

        # make a file of the regions with overlaps
        bedtools coverage -a wholealpha.bed -b smallMerged.bed -d | awk '{if ($5 >= 2) {print $1, ($2+$4-1), ($2+$4)}}' OFS='\t' | bedtools sort -i stdin | bedtools merge > overlapping.smallMerged.bed
        # add them to the final overlaps file 
        cat overlapping.smallMerged.bed > ~{fName}.overlaps_resolved.bed
        cat ~{fName}.overlaps_resolved.bed

        # intersect the overlaps bed file 
        bedtools intersect -wb -a smallMerged.bed -b overlapping.smallMerged.bed | awk '{print $1,$2,$3,$4,$5,$6,$2,$3,$9}' OFS='\t' | bedtools sort -i > mixedOverlaps.bed

        # extract the names of overlapping alpha annotations
        grep "(" mixedOverlaps.bed > mixedOverlaps.toRename.bed || true 
        sed 's/.*(\(.*\))/\1/' mixedOverlaps.toRename.bed | awk '{print $1}' > alpha_names.txt

        # cleaning up the alpha HOR names
        awk 'FNR==NR{a[NR]=$1;next}{$4=a[FNR]}1' alpha_names.txt mixedOverlaps.toRename.bed | awk '{print $1, $2, $3, $4, $5,$6,$2,$3,$9}' OFS='\t' > mixedOverlaps.named.bed
        grep -v "(" mixedOverlaps.bed >> mixedOverlaps.named.bed || true
        bedtools sort -i mixedOverlaps.named.bed | uniq > mixedOverlaps.named.sorted.bed

        # merge the overlapping regions
        bedtools merge -c 4 -o distinct -i mixedOverlaps.named.sorted.bed > mixedOverlaps.merged.bed

        # format bed entries 
        sed 's/$/\t0\t.\t.\t.\t204,0,0/' mixedOverlaps.merged.bed > mixedOverlaps.merged.named.bed
        awk '$7=$2' OFS='\t' mixedOverlaps.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="mixedAlpha("$4")"' OFS='\t' > mixedOverlaps.part.bed
        bedtools subtract -a smallMerged.bed -b mixedOverlaps.part.bed | awk '{print $1,$2,$3,$4,$5,$6,$2,$3,$9}' OFS='\t'  > overlapping.filtered.part.bed

        # create the final alpha file 
        cat mixedOverlaps.part.bed  > overlapsResolved.alpha.bed
        cat overlapping.filtered.part.bed >> overlapsResolved.alpha.bed
        bedtools sort -i overlapsResolved.alpha.bed > overlapsResolved.alpha.sorted.bed

        # merge all annotations into one file and sort 
        for f in *part.bed ; do cat $f >> ~{fName}.cenSat.bed ; done 
        cat overlapsResolved.alpha.sorted.bed >> ~{fName}.cenSat.bed
        cat ~{rDNABed} >> ~{fName}.cenSat.bed
        bedtools sort -i ~{fName}.cenSat.bed | uniq | awk '($3-$2) >= 800' > ~{fName}.cenSat.sorted.bed

        

        # sort out entries smaller than 800 bp - removes single monomers etc
        # also fix that bedtools subtract only alters columns 2 and 3 and not 7 and 8 & subtract one from end because bedtools subtract leaves overlaps of 1 bp 
        # this will be fixed in the next steps 
        bedtools sort -i ~{fName}.cenSat.sorted.bed | awk '{print $1, $2, $3, $4, $5,$6,$2, $3,$9}' OFS='\t' > ~{fName}.sorted.bed

        # close gaps smaller than 2000 bp - avoid tiny CT annotations
        # this closes gaps by expanding the annotation upstream
        bedtools closest -io -D a -iu -a ~{fName}.sorted.bed -b ~{fName}.sorted.bed | awk ' BEGIN {OFS="\t"} {if ($19 > 0 && $19 < 2000) ($3=$8=($8+$19-1))} {print $1,$2,$3,$4,$5,$6,$2,$3,$9 }' > tmp.txt && mv tmp.txt ~{fName}.sorted.bed
        
        # now add the gap annotations - these override any existing annotation 
        cat ~{gapBed} | awk '($3-$2) >= 1'  > ~{gapBed}.filtered
        bedtools subtract -a ~{fName}.sorted.bed -b ~{gapBed}.filtered | awk '{print $1, $2, $3, $4, $5,$6, $2, $3 ,$9}' OFS='\t'> ~{fName}.gap.merged.bed
        cat ~{gapBed} | awk '{print $1, $2, $3, $4, $5,$6, $2, $3 ,$9}' OFS='\t' >> ~{fName}.gap.merged.bed
        bedtools sort -i ~{fName}.gap.merged.bed > ~{fName}.sorted.bed

        # create the CT annotation and define centromere intervals
        bedtools merge -d 2000000 -i ~{fName}.sorted.bed > centromeres.bed
        grep active ~{aSatBed} > active_arrays.bed || true
        bedtools intersect -wa -u -a centromeres.bed -b active_arrays.bed > ~{fName}.active.centromeres.bed
        awk '{print $1, $2, $3}'  OFS='\t' ~{fName}.active.centromeres.bed > active_centromeres.bed
        bedtools coverage -a active_centromeres.bed -b ~{fName}.sorted.bed -d | awk '{ if ($5 == 0) { print $1, ($2+$4-1), ($2+$4)} }' OFS='\t' | bedtools merge | awk '{print $1, $2, $3}' OFS='\t' > CT.bed

        sed 's/$/\tct\t0\t.\t.\t.\t224,224,224/' CT.bed > CT.named.bed
        awk '$7=$2' OFS='\t' CT.named.bed | awk '$8=$3' OFS='\t' >> ~{fName}.sorted.bed

        # finalize bed file 
        echo 'track name="'~{fName}'" visibility=2 itemRgb="On"' > ~{fName}.cenSat.bed
        bedtools sort -i ~{fName}.sorted.bed >> ~{fName}.cenSat.bed

        # Finalize the strand track 
        # first let's add the strand information into our strand file 
        bedtools merge -s -c 4,6 -o distinct -s -d 500 -i ~{aSatStrand} | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "AS_strand", "0", $5, $2, $3, "."}' >> strandInfo.bed 
        # Now sort and rename the entries 
        bedtools sort -i strandInfo.bed | awk '($3-$2) >= 1000' > strandInfo.sorted.bed
        echo 'track name="'~{fName}'_Satellite_Strand" visibility=2 itemRgb="On"' > ~{fName}.cenSatStrand.bed
        cat strandInfo.sorted.bed | awk 'BEGIN{OFS="\t"} {if ($6 == "+") {($4=($4"_Plus_Strand")) && ($9="0,0,255")} else {($4=($4"_Minus_Strand")) && ($9="255,0,0")} print}' >> ~{fName}.cenSatStrand.bed
        

        #clean up the directory 
        # rm CT*bed
        # rm active_arrays.bed
        # rm *SAT*bed
        # rm *merged*bed
        # rm *sorted.bed
        # rm *filtered.bed
        # rm *part.bed
        # rm mixed*bed

    >>> 
    output {
        File cenSatAnnotations = "~{fName}.cenSat.bed"
        File cenSatStrand = "~{fName}.cenSatStrand.bed"
        File centromeres = "~{fName}.active.centromeres.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}
