version 1.0

workflow marmoset_ubam_to_methylbam {
	input {
		Array[File] ubams
        File ref_fasta
		String sample_id

		# typically either "map-ont" or "map-pb"
		String minimap2_preset = "map-ont"
	}

    call index_ref_fasta {
        input:
            ref_fasta = ref_fasta
    }

    scatter (ubam in ubams) {
        call samtools_fastq {
            input:
                ubam = ubam
        }

        call minimap2_align {
            input:
                ref_fasta = index_ref_fasta.fa,
                fq = samtools_fastq.fq,
                minimap2_preset = minimap2_preset
        }
    }

	call samtools_merge_index {
        input:
			sams = minimap2_align.sam,
			sample_id = sample_id
    }

    call modkit_pileup {
        input:
            ref_fa = index_ref_fasta.fa,
            ref_fa_fai = index_ref_fasta.fai,
			bam = samtools_merge_index.bam,
            bai = samtools_merge_index.bam_index,
			sample_id = sample_id
    }

    call bedgraph_to_bigwig {
        input:
            pileup_bedgraph = modkit_pileup.modkit_pileup_bedgraph,
            fa_genome = index_ref_fasta.fa_genome,
    }

	output {
        File bam = samtools_merge_index.bam
        File bai = samtools_merge_index.bam_index
        File bedMethyl = modkit_pileup.modkit_pileup_bed
        File bedgraph_5mC = modkit_pileup.modkit_pileup_bedgraph
        File bigwig_5mC = bedgraph_to_bigwig.pileup_bigwig
	}

	meta {
		author: "Julian Menendez"
		email: "jmmenend@ucsc.edu"
		description: "Maroset Workflow (align, pileup, make bigwig...)"
	}
}

task index_ref_fasta {
	input {
		File ref_fasta

		Int threadCount = 16
		Int memSizeGB = 128
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(ref_fasta, "GB"))       
	Int final_disk_size = input_file_size * 12

    String fai_output_prefix = basename(ref_fasta)

	command <<<
		set -eux -o pipefail

        samtools faidx --fai-idx "~{fai_output_prefix}.fai" ~{ref_fasta}
		cut -f 1,2 "~{fai_output_prefix}.fai" > "~{fai_output_prefix}.genome"
	>>>

	output {
        File fa = ref_fasta
		File fai = "~{fai_output_prefix}.fai"
        File fa_genome = "~{fai_output_prefix}.genome"
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_size + " SSD"
		docker: "staphb/samtools:1.21" 
		preemptible: 1
	}
}

task samtools_fastq {
	input {
		File ubam

		Int threadCount = 16
		Int memSizeGB = 128
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(ubam, "GB"))       
	Int final_disk_size = input_file_size * 12

    String ubam_basename = basename(ubam)
    String fq_output = sub(ubam_basename, "\\.bam$", ".fq")

	command <<<
		set -eux -o pipefail

		samtools fastq \
            -T '*' \
            -@ ~{threadCount} \
            ~{ubam} > ~{fq_output}
	>>>

	output {
		File fq = fq_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_size + " SSD"
		docker: "staphb/samtools:1.21" 
		preemptible: 1
	}
}

task minimap2_align {
	input {
        File fq
		File ref_fasta

		String minimap2_preset

		Int threadCount = 32
		Int memSizeGB = 256
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(fq, "GB")) + ceil(size(ref_fasta, "GB"))   
	Int final_disk_dize = input_file_size * 6

	# set outputs as wdl variables
    String fq_basename = basename(fq)
    String aligned_sam_output_prefix = sub(fq_basename, "\\.fq$", "")

	command <<<
		set -eux -o pipefail

		minimap2 \
			-I 16G \
			--eqx \
            --cs \
			-Y \
			-L \
            -y \
			-p 0.5 \
            -k 17 \
            -ax ~{minimap2_preset} \
            -t ~{threadCount} \
            ~{ref_fasta} ~{fq} > ~{aligned_sam_output_prefix}.minimap2.sam
	>>>

	output {
		File sam = "~{aligned_sam_output_prefix}.minimap2.sam"
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/minimap2:2.28" 
		preemptible: 1
	}
}

task samtools_merge_index {
	input {
		Array[File] sams
		String sample_id

		Int threadCount = 32
		Int memSizeGB = 256
	}

	# Estimate disk size required
	Int input_sam_size = ceil(size(sams, "GB")) 
	Int final_disk_dize = input_sam_size * 12

    File sam_header = sams[0]
    String merged_sam = "~{sample_id}.sam"
    String merged_bam_output = "~{sample_id}.bam"
    String merged_bai_output = "~{sample_id}.bam.bai"

	command <<<
		set -eux -o pipefail

        samtools merge \
            -@ ~{threadCount} \
            -o ~{merged_sam} \
            -h ~{sam_header} \
            ~{sep=" " sams}

		samtools sort \
			-@ ~{threadCount} \
			~{merged_sam} > ~{merged_sam}.tmp && mv ~{merged_sam}.tmp ~{merged_sam}

        samtools view \
            -@ ~{threadCount} \
            -bh \
            ~{merged_sam} > ~{merged_bam_output}

        samtools index -@ ~{threadCount} ~{merged_bam_output}

	>>>

	output {
		File bam = merged_bam_output
		File bam_index = merged_bai_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/samtools:1.21" 
		preemptible: 1
	}
}

task modkit_pileup {
	input {
		File ref_fa 
		File ref_fa_fai

		File bam
		File bai

		String sample_id

		String modkit_thresholds = "--mod-thresholds m:0.8 --filter-threshold C:0.5"

		Int threadCount = 32
		Int memSizeGB = 256
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(bam, "GB"))       
	Int final_disk_size = input_file_size * 3

	String modkit_pileup_bed_output = "~{sample_id}.CpG.pileup.bed"
    String modkit_pileup_bedgraph_output = "~{sample_id}.5mC.bedgraph"

	String ref_basename = basename(ref_fa)
	String ref_fai_basename = basename(ref_fa_fai)

	command <<<
		set -eux -o pipefail

		ln -s ~{ref_fa} ./~{ref_basename}
        ln -s ~{ref_fa_fai} ./~{ref_fai_basename}

		modkit pileup \
			-t ~{threadCount} \
			--force-allow-implicit \
			--cpg \
			--ref ./~{ref_basename} \
			--combine-strands \
			~{modkit_thresholds} \
            ~{bam} ~{modkit_pileup_bed_output}

        # use awk to create bedgraph
        awk -v OFS='\t' '{if ($4=="m") print $1, $2, $3, $11}' ~{modkit_pileup_bed_output} > ~{modkit_pileup_bedgraph_output}
	>>>

	output {
		File modkit_pileup_bed = modkit_pileup_bed_output
        File modkit_pileup_bedgraph = modkit_pileup_bedgraph_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_size + " SSD"
		docker: "jmmenend/modkit:0.4.2"
		preemptible: 1
	}
}

task bedgraph_to_bigwig {
	input {
		File pileup_bedgraph
		File fa_genome

		Int threadCount = 32
		Int memSizeGB = 256
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(pileup_bedgraph, "GB"))     
	Int final_disk_size = input_file_size * 6

	String bedgraph_basename = basename(pileup_bedgraph, ".bedgraph")
	String pileup_bigwig_output = "~{bedgraph_basename}.bigwig"

	command <<<
		set -eux -o pipefail

		bedGraphToBigWig \
			~{pileup_bedgraph} \
			~{fa_genome} \
			~{pileup_bigwig_output}
	>>>

	output {
		File pileup_bigwig = pileup_bigwig_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_size + " SSD"
		docker: "jmmenend/bedgraphtobigwig:latest"
		preemptible: 1
	}
}