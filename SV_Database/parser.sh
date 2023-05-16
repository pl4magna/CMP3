#!/bin/bash

# read config
cat << EOF > /hpcshare/genomics/plamagna/sv_backup/config_clustering
$(printf "%s\n" $(grep 'clust' -A3 sv_backup/code/config.ini | tail -n+2))
EOF

source /hpcshare/genomics/plamagna/sv_backup/config_clustering

_WORKING_DIR=$path_out


cd ${_WORKING_DIR}"/clusters"

len_clusters=$(find . -iname "cluster_*" | xargs ls -lrt | wc -l)

cd ${_WORKING_DIR}


launch(){

	list_jobs=()

	for n in ${subset[@]};do

		JOB=$(echo "singularity exec --bind $HOUSE,/hpcshare/genomics/afant/,/archive/s2/genomics/plamagna/adni/ /hpcshare/genomics/plamagna/sv_backup/container/clustering3.sif python /hpcshare/genomics/plamagna/sv_backup/code/retrieve_GT.py -c /hpcshare/genomics/plamagna/sv_backup/code/config.ini -k ${n},$((n+${range}))" | qsub -q fatnodes -e /hpcshare/genomics/plamagna/sv_backup/code/retrieve_GT.err -l walltime=100:00:00 -N retrieve_gt_${n}_$((n+${range}))) 

		#read -r "JOB" < <(echo "${subset[$((n+=1))]}")
		#JOB=$(echo "${subset[$((n+=1))]}")
		#echo $JOB
		list_jobs+=(${JOB})

	done

	echo ${list_jobs[@]}
}


if (( len_clusters % 10 != 0 ));then

	range=$(( (len_clusters - len_clusters % 10) / 10 ))
	subset=$(seq 0 ${range} ${len_clusters})
	subset=($(for x in $subset ${len_clusters};do echo $x;done))

	launch

else

	range=$(( len_clusters / 10 ))
	subset=($(seq 0 ${range} ${len_clusters}))

	launch

fi


echo "awk 'FNR>1' ${_WORKING_DIR}/genotypes_df_* >> ${_WORKING_DIR}/genotypes_df && sed -i '1 i id\tsamplename\tGT' ${_WORKING_DIR}/genotypes_df" | qsub -W depend=afterok:${list_jobs[@]} -N parse_genotypes_df
echo "awk 'FNR>1' ${_WORKING_DIR}/variants_df_* >> ${_WORKING_DIR}/variants_df && sed -i '1 i id\tchrom\tstart\tend\talt\tcompid\ttargeted_allele_freq' ${_WORKING_DIR}/variants_df" | qsub -W depend=afterok:${list_jobs[@]} -N parse_variants_df
echo "awk 'FNR>1' ${_WORKING_DIR}/clusters_df_* >> ${_WORKING_DIR}/clusters_df && sed -i '1 i clusters\tvariants' ${_WORKING_DIR}/clusters_df" | qsub -W depend=afterok:${list_jobs[@]} -N parse_clusters_df







