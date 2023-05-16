#!/usr/bin/python3.6

""" LAUNCH
echo "singularity exec --bind $HOUSE $BIN/pyodbc.sif python3 /hpcshare/genomics/plamagna/sv_backup/code/connect_to_database.py -a new -o /hpcshare/genomics/plamagna/sv_backup/output/PPMI_ADNI/" | qsub -N connect_to_db -e /hpcshare/genomics/plamagna/sv_backup/code/connect_to_database.err -l walltime=100:00:00
"""

# psql -h 10.244.0.91 -U vargenius vargenius_db

import os, sys, time
import argparse as ap
from getpass import getpass
import pyodbc
import pandas as pd


# get arguments
parser = ap.ArgumentParser()
parser.add_argument("-s", "--server", default="10.244.0.91", type=str, help='server name or ip..')
parser.add_argument("-d", "--database", default="vargenius_db", type=str, help='name of SQL database..')
parser.add_argument("-u", "--username", default="vargenius", type=str, help='postgres user ID..')

parser.add_argument("-a", "--action", choices=["new", "pull", "upgrade", "ask", "annotate"], required=True, type=str, help='[new | pull | upgrade | cluster_on_tsv] to push new clusters, pull variants and genotypes tables, just upgrade them or cluster SVs within a TSV..')

parser.add_argument("-ty", "--type", type=str, help='')
parser.add_argument("-w", "--which", type=str, help='')
parser.add_argument("-t", "--tsv_path", type=str, help='path of tsv to annotate')
parser.add_argument("-c", "--clusters_path", type=str, help='path to the file that contains SVs for each cluster')
parser.add_argument("-tb", "--table", type=str, help='table to download when -a pull')
parser.add_argument("-o", "--output", type=str, help='output path')
args = parser.parse_args()

# password = getpass("password for postgres server: ")
password = "vargenius_pwd"


# set global variables
_WORKING_DIR="/hpcshare/genomics/plamagna/sv_backup/output/PPMI_ADNI/"
VAR_TYPE = args.type
VARIANTS_TO_APPEND=_WORKING_DIR+"variants_df"
# VARIANTS_TO_ADD="/hpcshare/genomics/plamagna/to_add"
GENOTYPES_TO_APPEND=_WORKING_DIR+"genotypes_df"
CLUSTERS_TO_APPEND=_WORKING_DIR+"clusters_df"
CLUSTERED_SAMPLES_PATH="/hpcshare/genomics/plamagna/sv_backup/path_list"


os.chdir(_WORKING_DIR)


class Connect_To_Database:

	def __init__(self, args):

		# set attributes
		self.server = args.server
		self.database = args.database 
		self.username = args.username


	# connect to postgres
	def connect(self):

		cnxn = pyodbc.connect('DRIVER={PostgreSQL ANSI};SERVER='+self.server+';DATABASE='+self.database+';UID='+self.username+';PWD='+ password)
		cursor = cnxn.cursor()
		cnxn.setdecoding(pyodbc.SQL_CHAR, encoding='utf8')
		cnxn.setdecoding(pyodbc.SQL_WCHAR, encoding='utf8')
		cnxn.setencoding(encoding='utf8')
		setattr(self, 'cnxn', cnxn)
		setattr(self, 'cursor', cursor)


	# insert dataframes into SQL Server (all_col=T) or upgrade table col (all_col=F):
	def push(self, table, df, all_col=True):

		if all_col:

			if table == f"variants_sv":

				# cols = ", ".join([str(i) for i in df.columns])
				cols = ["\"end\"" if col == "end" else col for col in df.columns]
				cols = ", ".join([i for i in cols])
				len_col = '?,' * len([str(i) for i in df.columns])

				for index, row in df.iterrows():

					values = ", ".join([str(value) for value in list(row.values)])

					statement = f"INSERT INTO {table} ({cols}) VALUES ({len_col[:-1]});"

					self.cursor.execute(statement, row.cluster_id, row.chrom, row.start, row.end, row.SV_type, row.compid, row.WGS_allele_freq, row.WGS_freq_factors)		
					self.cnxn.commit()

				self.cursor.close()

			if table == f"genotypes_sv":

				# cols = ", ".join([str(i) for i in df.columns])
				cols = ["\"end\"" if col == "end" else col for col in df.columns]
				cols = ["\"GT\"" if col == "GT" else col for col in df.columns]
				cols = ", ".join([i for i in cols])
				len_col = '?,' * len([str(i) for i in df.columns])

				for index, row in df.iterrows():

					values = ", ".join([str(value) for value in list(row.values)])

					statement = f"INSERT INTO {table} ({cols}) VALUES ({len_col[:-1]});"

					self.cursor.execute(statement, row.cluster_id, row.sampleid, row.GT)
					self.cnxn.commit()

				self.cursor.close()

			if table == f"clusters":

				cols = ", ".join([str(i) for i in df.columns])
				len_col = '?,' * len([str(i) for i in df.columns])

				for index, row in df.iterrows():

					values = ", ".join([str(value) for value in list(row.values)])

					statement = f"INSERT INTO {table} ({cols}) VALUES ({len_col[:-1]});"

					self.cursor.execute(statement, row.cluster, row.variants)
					self.cnxn.commit()

				self.cursor.close()

		else:

			for index, row in df.iterrows():

				values = [str(value) for value in list(row.values)]

				statement = f"UPDATE {table} SET wgs_allele_freq = (?) WHERE id = (?);"

				self.cursor.execute(statement, values[1], values[0])		
				self.cnxn.commit()

			self.cursor.close()


	# pull SQL table into a dataframe
	def pull(self, table, df_name):

		df = pd.read_sql(f"SELECT * FROM {table};", self.cnxn)
		df.to_csv(path_or_buf=args.output+df_name, index=False)


	# intersect genotypes table
	def join_tables(self, samples, variants, genotypes):
		
		samples = pd.read_csv(samples, sep=',')
		variants = pd.read_csv(variants, sep=',')

		genotypes.columns = ["cluster_id", "samplename", "GT"]

		merged1 = pd.merge(genotypes, variants, on=["cluster_id"])
		merged1 = merged1[["varid", "samplename", "GT"]]
		merged2 = pd.merge(merged1, samples, on=["samplename"], how="left")

		merged2 = merged2[["varid", "sampleid", "GT"]]
		merged2.columns = ["cluster_id", "sampleid", "GT"]
		merged2["sampleid"] = merged2["sampleid"].fillna("0")
		merged2["sampleid"] = merged2["sampleid"].astype(int)
		# merged2.info()


		global geno_df
		geno_df = merged2

		return geno_df


	def add_freq_to_tsv(self, tsv, clusters):

		tsv_df = pd.read_csv(tsv)
		clusters_df = pd.read_csv(clusters)		


	# create new SQL database and tables with clusters info
	def create_tables(self):

  		# create variants_sv table
		self.cursor.execute(f'''
			CREATE TABLE variants_sv (
			varid serial PRIMARY KEY,
			cluster_id VARCHAR(10000),
			chrom VARCHAR(50),
			start int,
			"end" int,
			SV_type VARCHAR(10000),
			compid VARCHAR(1000),
			WGS_allele_freq decimal,
			WGS_freq_factors VARCHAR(1000)
			);
			''')

		self.cnxn.commit()

		# create genotypes_sv table
		self.cursor.execute(f'''
			CREATE TABLE genotypes_sv (
			analysisid int,
			cluster_id int,
			sampleid int,			
			"GT" VARCHAR(3),
			FOREIGN KEY(analysisid) REFERENCES analyses(analysisid)
			);
	        ''')

		self.cnxn.commit()
		
		# create clusters_sv table
		# self.cursor.execute(f'''
		# 	CREATE TABLE clusters_sv_{VAR_TYPE} (
		# 	cluster VARCHAR(1000),
		# 	variants VARCHAR(1000)
		# 	);
		# 	''')

		# self.cnxn.commit()


def merge_var_types():

	global var_df
	global geno_df
	
	var_df_del = pd.read_csv(_WORKING_DIR+"output_del/variants_df", sep="\t")
	geno_df_del = pd.read_csv(_WORKING_DIR+"output_del/genotypes_df", sep="\t")

	var_df_dup = pd.read_csv(_WORKING_DIR+"output_dup/variants_df", sep="\t")
	geno_df_dup = pd.read_csv(_WORKING_DIR+"output_dup/genotypes_df", sep="\t")

	var_df_ins = pd.read_csv(_WORKING_DIR+"output_ins/variants_df", sep="\t")
	geno_df_ins = pd.read_csv(_WORKING_DIR+"output_ins/genotypes_df", sep="\t")

	var_df = pd.concat([var_df_del, var_df_dup, var_df_ins], ignore_index = True, sort = False)
	geno_df = pd.concat([geno_df_del, geno_df_dup, geno_df_ins], ignore_index = True, sort = False)

	return var_df,geno_df


def calculate_frequency(var, geno, out="var_df"):

	with open(CLUSTERED_SAMPLES_PATH, "r") as path:
		tot_gen = len(path.readlines())
		path.close()

	V_het = geno[geno["GT"] == "0/1"].groupby("cluster_id")["cluster_id"].count()
	V_homo = geno[geno["GT"] == "1/1"].groupby("cluster_id")["cluster_id"].count()

	df = pd.DataFrame({"tot_gen": tot_gen, "V_het": V_het, "V_homo": V_homo}).fillna(0)
	df.index.name = 'cluster_id'
	df.reset_index(inplace=True)

	df["WGS_allele_freq"] = (df["V_het"] + (df["V_homo"] * 2)) / (tot_gen * 2)

	# add WGS_freq_factors column
	def add(a, b, c):
		return str(int(a)) + "]---[" + str(int(b)) + "]---[" + str(tot_gen - int(a) - int(b))

	df['WGS_freq_factors'] = df.apply(lambda row : add(row['V_homo'], row['V_het'], row['tot_gen']), axis = 1)

	# intersection
	merged = pd.merge(var, df[["cluster_id", "WGS_allele_freq", "WGS_freq_factors"]], on=["cluster_id"])
	merged.drop("WGS_allele_freq_x", axis=1, inplace=True)
	merged.columns = list(merged.columns[:-2]) + ["WGS_allele_freq"] + ["WGS_freq_factors"]


	global var_df
	var_df = merged

	var_df.to_csv("/hpcshare/genomics/plamagna/sv_backup/output/PPMI_ADNI/new_variants", index=False)
	
	if out == "var_df":
		return var_df
	elif out == "freq_col":
		var_df = var_df[["cluster_id", "WGS_allele_freq"]]
		return var_df



# initialize


# instantiate Connect_To_Database class
connection = Connect_To_Database(args)
connection.connect()

if args.action == "new":

	# merging variants and genotypes from different SV types
	merge_var_types()

	# compute frequency
	calculate_frequency(var=var_df, geno=geno_df)

	# connect && push clusters dataframe into a new database
	connection.create_tables()
	connection.pull("samples", df_name="samples_df")
	connection.push(table=f"variants_sv", df=var_df)

	time.sleep(10)

	connection.pull(f"variants_sv", df_name=f"variants_sv")

	connection.join_tables(samples=f"samples_df", variants=f"variants_sv", genotypes=geno_df)

	connection.push(table=f"genotypes_sv", df=geno_df)


elif args.action == "pull":

	connection.pull(args.table, df_name=f"{args.table}_df")


elif args.action == "upgrade":

	# download variants table & run add_to_cluster before

	_WORKING_DIR = args.output

	os.chdir(_WORKING_DIR)

	# load genotypes dataframe to upload
	geno_df = pd.read_csv(f"{_WORKING_DIR}genotypes_to_append", sep="\t")

	connection.pull("samples", df_name="samples_df")
	connection.pull(f"variants_sv_{args.type}", df_name=f"variants_sv_{args.type}_pulled")

	connection.join_tables(samples="samples_df", variants=f"variants_sv_{args.type}_pulled", genotypes=geno_df)

	# update new genotypes
	connection.push(table=f"genotypes_sv_{args.type}", df=geno_df)

	# pull variants table
	connection.pull(f"variants_sv_{args.type}", df_name="variants_sv_pulled")
	connection.pull(f"genotypes_sv_{args.type}", df_name="genotypes_sv_pulled")

	# compute frequency
	var_df = pd.read_csv(_WORKING_DIR+"variants_sv_pulled", sep=",")
	geno_df = pd.read_csv(_WORKING_DIR+"genotypes_sv_pulled", sep=",")
	calculate_frequency(var_df, geno_df, out="freq_col")

	# push updated frequencies
	connection.push(table=f"variants_sv_{args.type}", df=var_df, all_col=False)


elif args.action == "ask":

	_WORKING_DIR = "/hpcshare/genomics/plamagna/sv_backup/output/PPMI_ADNI/interrogated/"

	variants_to_interrogate = args.which.split(",")

	vartype = tuple([w.split("_")[3] for w in variants_to_interrogate])

	# creating clusters bed files
	for ty in vartype:

		connection.pull(f"variants_sv_{ty}", df_name=f"variants_sv_{ty}_pulled")
		connection.pull(f"genotypes_sv_{ty}", df_name=f"genotypes_sv_{ty}_pulled")

		bed_clusters = pd.read_csv(f"{_WORKING_DIR}variants_sv_{ty}_pulled", usecols=["chrom", "start", "end"])

		bed_clusters.to_csv(path_or_buf=f"{_WORKING_DIR}clusters_{ty}.bed", sep="\t", header=False, index=False)

	# creating variants bed files
	for vr in variants_to_interrogate:

		with open(f"{_WORKING_DIR}{vr}.bed", "w") as bed:

				chrom = vr.split("_")[0]
				start = vr.split("_")[1]
				end = vr.split("_")[2]
			
				bed.write(f"chr{chrom}\t{start}\t{end}\n")


elif args.action == "annotate":

	tsvs = []

	with open(args.tsv_path, "r") as tsv_path:
		for line in clst.readlines():
		    tsvs.append(line)
		tsv_path.close()

	for tsv in tsvs:

		connection.add_freq_to_tsv(tsv=tsv, clusters=args.clusters_path)




