#!/bin/bash


i=$SLURM_ARRAY_TASK_ID

if [[ -s read_${i}-BC1.fastq.gz ]] && [[ -s read_${i}-BC2.fastq.gz ]]; then
	echo "read_${i}-BC1.fastq.gz and read_${i}-BC2.fastq.gz already exist."
	echo "Skip processing."
	exit 0
fi

echo "
$tools/scaff10x/Scaff10X-4.1/src/scaff_reads -nodes $SLURM_CPUS_PER_TASK input$i.dat read_${i}-BC1.fastq.gz read_${i}-BC2.fastq.gz"
$tools/scaff10x/Scaff10X-4.1/src/scaff_reads -nodes $SLURM_CPUS_PER_TASK input$i.dat read_${i}-BC1.fastq.gz read_${i}-BC2.fastq.gz

