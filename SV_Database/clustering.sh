#!/bin/bash

# qsub -q fatnodes -e sv_backup/output/clustering.err sv_backup/code/clustering.sh


# read config
cat << EOF > /hpcshare/genomics/plamagna/sv_backup/config_clustering
$(printf "%s\n" $(grep 'clust' -A3 sv_backup/code/config.ini | tail -n+2))
EOF

source /hpcshare/genomics/plamagna/sv_backup/config_clustering

_WORKING_DIR=$path_out
_TMPFILE_NCLS="tmp_cls"
_TMPFILE_NCHR="tmp_chr"

echo 0 > ${_WORKING_DIR}/${_TMPFILE_NCLS}

mkdir -p ${_WORKING_DIR}/clusters/
mkdir -p ${_WORKING_DIR}/chr_intersections/

rm /hpcshare/genomics/plamagna/sv_backup/config_clustering


find_inners(){ ## returns mean point from an array of coordinates

	arr=("$@")
	tot=0

	for i in ${arr[@]}; do
	  let tot+=$i
	done

	local mean_point=$(( $tot / ${#arr[@]} ))

	echo $mean_point
}


give_inners(){ ## returns inner starts and ends

	st_arr=$(cut -d " " -f3 $1)
	end_arr=$(cut -d " " -f4 $1)


	inner_start="$(find_inners ${st_arr[@]})"
	inner_end="$(find_inners ${end_arr[@]})"

	echo ${1##*/} ${c} $inner_start $inner_end >> ${_WORKING_DIR}/inners_clusters
}


### --------- MAIN --------- ###

for c in {1..22..1} X Y;do 

	# splitting group dataframe per chr
	grep  "	${c}	"  ${_WORKING_DIR}/groups_df_${var_type} | cut -f2,3,4 > ${_WORKING_DIR}/chr_intersections/chr$c

	echo ${c} > ${_WORKING_DIR}/${_TMPFILE_NCHR}

	while true;do

		c=$( cat ${_WORKING_DIR}/${_TMPFILE_NCHR} )

		sed -n 1p ${_WORKING_DIR}/chr_intersections/chr$c > ${_WORKING_DIR}/bed
		
		go bedtools intersect -a ${_WORKING_DIR}/chr_intersections/chr$c -b ${_WORKING_DIR}/bed -f 0.7 -r -wa -wb > ${_WORKING_DIR}/intersections

		n=$[$(cat ${_WORKING_DIR}/${_TMPFILE_NCLS}) + 1]
		echo $n > ${_WORKING_DIR}/${_TMPFILE_NCLS}

		cluster=cluster_${n}_${var_type}
		
		if [ ${#cluster} -eq 13 ];then
			cluster=${cluster:0:7}_0${cluster:8}
		fi

		# while read chr start1 end1 n start2 end2;do

		# 	echo $chr $start1 $end1 >> ${_WORKING_DIR}/outputs/clusters/${cluster}

		# 	sed -i "/${chr}\t${start1}\t${end1}/d" ${_WORKING_DIR}/outputs/chr_intersections/chr$c 

		# done < ${_WORKING_DIR}/outputs/intersections

		# writing clusters
		while read chr start1 end1 n start2 end2;do

			echo ${cluster} $chr $start1 $end1 >> ${_WORKING_DIR}/clusters/${cluster}

			grep -P --invert-match "${chr}\t${start1}\t${end1}" ${_WORKING_DIR}/chr_intersections/chr$c > ${_WORKING_DIR}/chr_intersections/chr$c_tmp

			mv ${_WORKING_DIR}/chr_intersections/chr$c_tmp ${_WORKING_DIR}/chr_intersections/chr$c

		done < ${_WORKING_DIR}/intersections


		give_inners ${_WORKING_DIR}/clusters/${cluster}

		#n=$((n + 1))
		#(( n+=1 ))
		#n=$[$n +1]

		if [ $( wc -l ${_WORKING_DIR}/chr_intersections/chr$c | cut -d " " -f1 ) -eq 0 ];then
			break
		fi

	done

done
