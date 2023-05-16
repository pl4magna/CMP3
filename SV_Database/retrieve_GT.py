#!/usr/bin/python3.6

""" LAUNCH
echo "singularity exec --bind $HOUSE,/archive/s2/genomics/plamagna/adni/,/hpcshare/genomics/afant/PPMI/ /hpcshare/genomics/plamagna/sv_backup/container/clustering3.sif python /hpcshare/genomics/plamagna/sv_backup/code/retrieve_GT.py -c /hpcshare/genomics/plamagna/sv_backup/code/config.ini -k" | qsub -e /hpcshare/genomics/plamagna/sv_backup/code/retrieve_GT.err -l walltime=100:00:00
"""


import os, re, sys, mmap
import pandas as pd
import argparse as ap
import configparser
import ast


parser = ap.ArgumentParser()
parser.add_argument("-c", "--configuration", type=str, default='config.ini', help='Parameters for configuration')
parser.add_argument("-k", "--coordinates", type=str, default='config.ini', help='SV coordinates , delimited')   
args = vars(parser.parse_args())  
config = configparser.ConfigParser()
config.read(args['configuration'])

_WORKING_DIR = str(ast.literal_eval(config['clust']['path_out']))
VAR_TYPE = str(ast.literal_eval(config['clust']['var_type']))

start = args["coordinates"].split(",")[0]
end = args["coordinates"].split(",")[1]

os.chdir(_WORKING_DIR)


# get list of full path clusters
def list_full_paths(directory):
	return [os.path.join(directory, file) for file in os.listdir(directory)]

clusters = list_full_paths(_WORKING_DIR+"clusters/")


# create dictionaries of genotypes per single variants (genotypes) and cluster info (variants)
genotypes = {
	"cluster_id": [],
	"samplename": [],
	"GT": []
}

variants = {
	"cluster_id": [],
	"chrom": [],
	"start": [],
	"end": [],
	"SV_type": [],
	"compid": [],
	"WGS_allele_freq": []
}

clusters_variants = {
	"cluster": [],
	"variants": []
}


def main(start, end):

	# loop over clusters
	for cluster in clusters[int(start):int(end)]:

		# getting supergroups info in current cluster
		supergroups = []

		with open(cluster, "r") as clst:

		  for line in clst.readlines():
		    supergroups.append("{},{},{}".format(line.split(" ")[1], line.split(" ")[2], line.split(" ")[3])[:-1])
		  clst.close()


		# retrieving variant ids within supergroups
		varid = []

		with open(_WORKING_DIR+"intervals_super_group_"+VAR_TYPE, "r", encoding="utf-8") as intervals_super_group:
			with mmap.mmap(intervals_super_group.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
				for supergroup in supergroups:
					sup = bytes(supergroup, encoding='utf-8')
					o = mmap_obj.find(sup)
					if o != -1:
						mmap_obj.seek(o)
						l = mmap_obj.readline()
						w = str(l).split(",")
						for i in w:
							if i[:14] == " \\'sv_id\\': \\'":
								if re.search("DNA", i[14:]):
									post = i[re.search("DNA", i).span()[0]-10:-2]
									pre = i[14:re.search("DNA", i).span()[0]-10]
									varid.append(pre + "-" + post)
								elif re.search("PPMI", i[14:]):
									post = i[re.search("PPMI", i).span()[0]:-2]
									pre = i[14:re.search("PPMI", i).span()[0]]
									varid.append(pre + "-" + post)
								elif re.search("ASL_ONC_", i[14:]) or re.search("ASL_NEU_", i[14:]) or re.search("ASL_NED_", i[14:]) or re.search("ASL_TRP_", i[14:]) or re.search("ASL_VDA_", i[14:]):
									post = i[re.search("ASL_", i).span()[0]:-2]
									pre = i[14:re.search("ASL_", i).span()[0]]
									varid.append(pre + "-" + post)
		    		# else:
		    		# 	print(o)
		    		# 	print(sup)
		    		# 	# raise Exception("no match")


		# updating variants dataframe with SVs within the cluster
		clusters_variants["variants"] = clusters_variants["variants"] + [v.split("-")[0] for v in varid]
		clusters_variants["cluster"] = clusters_variants["cluster"] + [cluster.split("/")[-1] for v in varid]


		# retrieving genotype info for each variant
		geno = []

		for var in varid:
			with open("/hpcshare/genomics/plamagna/sv_backup/path_list", "r") as path_list:
				for n, l in enumerate(path_list):
					if len(var.split('-')) == 3:
						if var.split('-')[2][:3] == "DNA":
							if var.split('-')[1] + "-" + var.split('-')[2] in l:
								if l[-1:] == "\n":
									l = l[:-1]
								with open(l, "r") as vcf:
									GT = []
									for line_n, line in enumerate(vcf):
										if var.split('-')[0] in line:
											s = line.split("\t")
											gt = s[14].split(":")[0]
											GT.append(gt)
									GT = list(set(GT))
									if len(GT) == 1:
										geno.append(GT)
									else:
										print(len(GT))
										print(var)
										print(l)
										raise Exception()
								vcf.close()
					elif len(var.split('-')) == 2:
						if var.split('-')[1] + "/" in l:
							if l[-1:] == "\n":
								l = l[:-1]
							with open(l, "r") as vcf:
								GT = []
								for line_n, line in enumerate(vcf):
									if var.split('-')[0] in line:
										s = line.split("\t")
										gt = s[14].split(":")[0]
										GT.append(gt)
								GT = list(set(GT))
								if len(GT) == 1:
									geno.append(GT)
								else:
									print(len(GT))
									print(var)
									print(l)
									raise Exception()
							vcf.close()


		# appending varid, sampleid and genotypes to genotypes dictionary
		genotypes["cluster_id"] = genotypes["cluster_id"] + [cluster.split("/")[-1] for v in varid]
		genotypes["samplename"] = genotypes["samplename"] + [s.split("-")[1] if len(s.split("-")) == 2 else s.split("-")[1] + "-" + s.split("-")[2] for s in varid]
		genotypes["GT"] = genotypes["GT"] + [s[x] for s in geno for x in range(0, len(s))]


		# filling variants dictionary as well 
		n_cls = cluster.split("/")[-1].split("_")[1]


		with open(_WORKING_DIR+"inners_clusters") as inners_clusters:
			coordinates = [[x] for x in [line.split(" ") for line_n, line in enumerate(inners_clusters) if "cluster_" + n_cls in line][0]]


		variants["cluster_id"] = variants["cluster_id"] + [cluster.split("/")[-1]]
		variants["chrom"] = variants["chrom"] + ["chr" + str(coordinates[1][0])]
		variants["start"] = variants["start"] + coordinates[2]
		variants["end"] = variants["end"] + [coordinates[3][0][:-1]]
		variants["SV_type"] = variants["SV_type"] + [VAR_TYPE]
		variants["compid"] = variants["compid"] + ["chr" + str(coordinates[1][0]) + "_" + str(coordinates[2][0]) + "_" + str(coordinates[3][0][:-1]) + "_" + VAR_TYPE]
		variants["WGS_allele_freq"] = variants["WGS_allele_freq"] + list("*")


	# create pandas dataframe from previous dictionaries
	genotypes_df = pd.DataFrame(genotypes)
	variants_df = pd.DataFrame(variants)
	clusters_df = pd.DataFrame(clusters_variants)

	genotypes_df.to_csv(_WORKING_DIR + "genotypes_df" + "_" + str(start) + "_" + str(end), sep = "\t", index = False)
	variants_df.to_csv(_WORKING_DIR + "variants_df" + "_" + str(start) + "_" + str(end), sep = "\t", index = False)
	clusters_df.to_csv(_WORKING_DIR + "clusters_df" + "_" + str(start) + "_" + str(end), sep = "\t", index = False)


main(start=start, end=end)
