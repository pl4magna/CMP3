#!/bin/bash

display_usage() {

	printf "%s\n" "usage example: go samtools view [bamfile].. "

	exit
}


if [[ $# -eq 0 ]] ; then

    display_usage

    exit 1

fi


if [[ ( $# == "--help") ||  $# == "-h" ]] ; then

    display_usage

    exit 0

fi



if [[ $# -eq 0 ]] ; then

    display_usage

    exit 1

fi


if [[ ( $# == "--help") ||  $# == "-h" ]] ; then

    display_usage

    exit 0

fi


case $1 in
	samtools) singularity run --bind $PWD,$3 $BIN/samtools_latest.sif samtools ${@:2};;
	bcftools) singularity run --bind $PWD,$ARCHIVE $BIN/bcftools-1.16.sif bcftools ${@:2};;
	fastqc) singularity exec --bind $PWD $BIN/fastqc_latest.sif fastqc ${@:2};;
	multiqc) singularity exec --bind $PWD /hpcshare/genomics/bioinfo_parabricks/SamBcfBedMultiqcTools.sif multiqc ${@:2};;
	plink) /hpcshare/genomics/bin/plink ${@:2};;
	bwa) singularity exec --bind $PWD $BIN/VGII_SingularityCont_v2.sif bwa ${@:2};;
	gatk) singularity exec --bind $PWD $BIN/gatk_4.1.9.0.sif ${@:2};;
	parabriks) exit 0;;
	picard) ;;
	manta) ;;
	melt) singularity exec --bind $PWD $BIN/MELT-2.2.2.sif ${@:2};;
	bedtools) singularity exec --bind $PWD $BIN/VGII_SingularityCont_v2.sif bedtools ${@:2} ;;
	vcftools) singularity exec --bind $PWD $BIN/vcftools.sif vcftools ${@:2} ;;
	--help | -h) display_usage;;
	*) printf "%s\n" "Please check the tool spell"; display_usage;;
esac

