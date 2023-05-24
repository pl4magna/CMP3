#!/bin/bash

variants=$1

_WORKING_DIR=""
BIN=""

mkdir -p ${_WORKING_DIR}


# display usage function
display_usage(){

echo
echo "Retrieves information from vargenius postgres database given SVs"
echo "Usage example:"
echo "sv_backup/code/which_cluster.sh 2_212258815_212314121_DUP_1,6_162478270_162548897_DEL_1"
echo

}


# display usage in cases
if [[ $1 == "-h" || $1 == "--help" ]];then

	display_usage
	exit

elif [[ -z $1 ]]; then
	
	display_usage
	exit

fi


# run connect_to_database.py
singularity exec python3 connect_to_database.py -a ask -w $variants -o $_WORKING_DIR


# check if more than one variant, in case use "cut"
if [[ ! -z "$(grep "," <<< $variants)" ]]; then
	variants=($(cut -d"," -f1- --output-delimiter " " <<< $variants))
fi


# loop through variants array
for var in ${variants[@]};do

	for f in $(seq 0.70 0.01 1.00);do

		while true;do

			cluster=$(singularity exec --bind ${_WORKING_DIR},$PWD $BIN/VGII_SingularityCont_v2.sif bedtools intersect -a ${_WORKING_DIR}${var}.bed -b ${_WORKING_DIR}clusters_$(cut -d"_" -f4 <<< $var).bed -f $f -r -wa -wb | wc -l)

			if ((${cluster} == 0));then
				
				printf "\n***\n%s not clusterable (no overlaps of 70%% or greater with any cluster) :( \n\n***\n\n" "$var" 
				break 2


			elif ((${cluster} > 1));then

				break 


			elif ((${cluster} == 1));then

				cluster=$(singularity exec --bind ${_WORKING_DIR},$PWD $BIN/VGII_SingularityCont_v2.sif bedtools intersect -a ${_WORKING_DIR}${var}.bed -b ${_WORKING_DIR}clusters_$(cut -d"_" -f4 <<< $var).bed -f $f -r -wa -wb | cut -f4,5,6 --output-delimiter=",")

				break 2

			fi

		done

	done

	# if no overlaps, continue to the next variant
	if ((${cluster} == 0)); then
		continue
	fi

	# retrieve information from postgres tables
	variants_sv_entry=$(grep $cluster "${_WORKING_DIR}variants_sv_$(cut -d"_" -f4 <<< $var)_pulled")

	varid=$(cut -d"," -f1 <<< $variants_sv_entry)

	genotypes=$(grep "^,${varid}," "${_WORKING_DIR}genotypes_sv_$(cut -d"_" -f4 <<< $var)_pulled")

	cluster_id=$(cut -d"," -f2 --output-delimiter=" " <<< $variants_sv_entry)


	# print cluster_id, position and cluster frequency output
	printf "\nVariant \"$var\" clusters in:\n"
	printf "\n%-20s %s" "cluster_id:" $cluster_id
	printf "\n%-20s %s-%s-%s" "position:" $(cut -d"," -f3,4,5 --output-delimiter=" " <<< $variants_sv_entry)
	printf "\n%-20s %f\n\n" "cluster frequency:" $(cut -d"," -f8 --output-delimiter=" " <<< $variants_sv_entry)


	# print genotypes and phenotypes output into a fancy table
	seperator=--------------------------------------------
	seperator=$seperator$seperator

	printf "%-$((${#var}+1))s| %-18s %s %-10s %s %s\n" cluster_id samplename "|" genotype "|" phenotype
	printf "%.70s\n" "$seperator"

	rows="%s | %-18s %s %-10s %s %s\n"

	while read out;do

		sampleid=$(cut -d"," -f3 <<< $out)
		geno=$(cut -d"," -f4 <<< $out)
		samplename=$(grep "^${sampleid}," samples_df | cut -d"," -f6)

		if [ ${samplename:0:4} == "PPMI" ];then
			pheno=$(grep ${samplename##*I} pheno_samples | cut -d" " -f1 )

		else
			pheno=$(cat CN EMCI LMCI AD | grep $samplename | cut -d" " -f1)
		fi

		printf "$rows" "$var" "$samplename" "|" "$geno" "|" "$pheno"
	
	done <<< $genotypes

	echo

done


