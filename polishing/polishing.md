#### Polishing marmoset verkko hi-c assemblies with DV

Alignment and DeepVariant submitted in batch here:

Alignment and PHARAOH:
https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows/hprc_DeepPolisher

DeepVariant:
https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows/deepvariant

Subset DV vcf to pass only variants
```
bcftools view -Oz -f "PASS" calJac240_verkko_hic_deepvariant.vcf.gz > calJac240_verkko_hic_deepvariant.PASS.vcf.gz
```

Separate unpolished assemblies by haplotype
```
docker run -it -v /private/groups:/private/groups nanozoo/seqkit /bin/bash

  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220/calJac220_final.oriented.fa | seqkit grep -r -p '.*_hap1' > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220/calJac220_final.hap1.fasta
  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220/calJac220_final.oriented.fa | seqkit grep -v -r -p '.*_hap1' > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220/calJac220_final_hap2.fasta

  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final.oriented.fa | seqkit grep -r -p '.*_hap1' > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final_hap1.fasta

  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final.oriented.fa | seqkit grep -r -p '.*chrX' >> /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final_hap1.fasta

  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final.oriented.fa | seqkit grep -r -p '.*_hap2' > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final_hap2.fasta

  cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final.oriented.fa | seqkit grep -r -p '.*chrY' >> /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final_hap2.fasta
```

Polish assemblies with DV pass vcf
```
bcftools consensus -H1 -f /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac240/calJac240_final.oriented.fa /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/calJac240_verkko_final/analysis/deepvariant_outputs/calJac240_verkko_final_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/calJac240_verkko_final/analysis/deepvariant_outputs/calJac240_final.oriented.DV_pass.fa

bcftools consensus -H1 -f /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220/calJac220_final.oriented.fa /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/calJac220_verkko_final/analysis/deepvariant_outputs/calJac220_verkko_final_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/calJac220_verkko_final/analysis/deepvariant_outputs/calJac220_final.oriented.DV_pass.fa
```

Rename final assemblies to avoid duplicates
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    pegi3s/seqkit seq --rename "s/(.*)/\1_primary/" /private/groups/cgl/pnhebbar/marmoset/patched_calJac220_v0305/final/calJac220_primary_v0305.fa \
    > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/calJac220_primary_v0305.fa
```
