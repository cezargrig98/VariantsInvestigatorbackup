#!/bin/env python

import sys
import argparse
import sqlite3
from sqlite3 import Error
import gzip


def log_to_file(message, file_path='file.log'):
    with open(file_path, 'a') as log_file:
        log_file.write(message + '\n')


## get and parse the arguments for the script
def get_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--vcf", help="vcf file to convert",
						action="store", dest = "vcf")
	parser.add_argument("-d", "--database", help="path to SQLite database - can be existing one",
						action="store", dest = "database")
	parser.add_argument("-v", "--variants_table", help="name of the table storing variant data - can be existing one",
						action="store", dest = "variants_table")
	parser.add_argument("-g", "--genotypes_table", help="name of the table storing genotype data - can be existing one",
						action="store", dest = "genotypes_table")
#	parser.add_argument("-l", "--log_file", help="path to log file", action="store", dest="log_file", default="error_log.txt")
	args = parser.parse_args()
	return(args)			

## parse info field into a dictionary
def parse_info_field(list_data):
	info_data = {}
	for data in list_data:
		if data.find("=") == -1: continue
		pair = data.split("=")
		info_data[pair[0]] = pair[1]
	return(info_data)

## parse sample data field into a dictionary
## dictated by the format field
def parse_sample_field(sample_data, format):
	geno_data = {}
	sample_data_sep = sample_data.split(":")
	for index in range(-1,len(format)-1):
		geno_data[format[index]] = sample_data_sep[index]
	return(geno_data) 

## parse each variant line
def parse_variant(string):

	fields = string.split("\t")

	## each line will return a dictionary structure 
	## this allows to access any object inside it with its key
	## and it can be also other lists or dictionaries as in the case of info or geno
	variant_data = {}

	variant_data['chr'] = fields[0]
	variant_data['start'] = fields[1]
	variant_data['end'] = fields[1]
	variant_data['id'] = fields[2]
	variant_data['ref'] = fields[3]
	variant_data['alt'] = fields[4].split(",")
	variant_data['qual'] = fields[5]
	variant_data['filter'] = fields[6]
    
	## use above function to parse info field data
	## returns a dictionary
	variant_data['info_data'] = parse_info_field(fields[7].split(";"))
	
	format = fields[8].split(":")
	sample_data_list = fields[9:]

    ## initialise sample values as a list
	## list contains as many dictionaries as many samples
	sample_values = []

	for sample_data in sample_data_list:
		sample_values.append(parse_sample_field(sample_data, format))
    
	variant_data['sample_values'] = sample_values

	
	## check if multiallelic and
	## if not, only uses the first value of the list
	type = ""
	alt = fields[4].split(",")
	if len(alt) > 1:
		type = "multiallelic"
		variant_data['width'] = 1
	elif len(variant_data['ref']) > len(alt[0]):
		type = "DEL"
		variant_data['width'] = len(variant_data['ref']) - len(alt[0])
	elif len(variant_data['ref']) < len(alt[0]):
		type = "INS"
		variant_data['width'] = len(alt[0]) - len(variant_data['ref'])
	elif len(variant_data['ref']) == len(alt[0]):
		type = "SNP"
		variant_data['width'] = 1

	## store this in dictionary
	variant_data['type'] = type
	
	return(variant_data)

## creates a connection with an SQLite database
def create_connection(db_file):
    conn = None
    try:
        con = sqlite3.connect(db_file)
        print(sqlite3.version)
        return(con)
    except Error as e:
        print(e)

## prepare SQL statement to create the necessary tables IF and only IF the tables do not exist
def create_tables(con, variants_table, genotypes_table):
    sql_statement_var = (
        "CREATE TABLE if not exists "
        + variants_table
        + "(id TEXT, chromosome TEXT, start NUM, "
        "end NUM, width NUM, ref TEXT, "
        "alt TEXT, qual NUM, filter TEXT, ac NUM, "
        "allele TEXT, consequence TEXT, impact TEXT, "
        "symbol TEXT, gene TEXT, feature TEXT, exon TEXT, intron TEXT, "
        "strand NUM, variantclass TEXT, sift TEXT, "
        "polyphen TEXT, gnomad_af NUM, "
        "clin_sig TEXT, lof TEXT, lof_filter TEXT, "
        "lof_flags TEXT, "
        "PRIMARY KEY(chromosome, start, ref, alt));"
    )

    sql_statement_gen = (
        "CREATE TABLE if not exists "
        + genotypes_table
        + "(samplename TEXT, chromosome TEXT, start NUM, end NUM, width NUM, "
        "ref TEXT, alt TEXT, "
        "gt TEXT, dp TEXT, gq TEXT, pl TEXT, ad TEXT,"
        "PRIMARY KEY(samplename, chromosome, start));"
    )

    cur = con.cursor()
    cur.execute(sql_statement_var)
    cur.execute(sql_statement_gen)
    con.commit()

def commit_data(con, variants_table, genotypes_table, variant_data, individuals):
	expected_info = ['AC', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature', 'EXON', 'INTRON', 
                     'STRAND', 'VARIANT_CLASS', 'SIFT', 'PolyPhen', 'gnomAD_AF', 
                     'CLIN_SIG', 'LoF', 'LoF_filter', 'LoF_flags']
	variant_info = variant_data['info_data']    
	for key in expected_info:
		if key not in variant_info:
			variant_info[key] = 'NA'  
	variant_info = variant_data['info_data'] # INSERT OR IGNORE INTO
	sql_variant = ("INSERT OR IGNORE INTO " + variants_table +
			        "(id, chromosome, start, "
					"end, width, ref, "
					"alt, qual, filter, AC, "
					"allele, consequence, impact, "
					"symbol, gene, feature, exon, intron, "
					"strand, variantclass, sift, "
					"polyphen, gnomad_af, "
					"clin_sig, lof, lof_filter, "
					"lof_flags) "
					"VALUES ("
					+ "\'" +  variant_data['id'] + "\'" +  ", "
					+ "\'" +  variant_data['chr'] + "\'" +  ", "
					+ "\'" +  variant_data['start'] + "\'" +  ", "
					+ "\'" +  variant_data['end'] + "\'" +  ", "
					+ "\'" +  str(variant_data['width']) + "\'" +  ", "
					+ "\'" +  variant_data['ref'] + "\'" +  ", "
					+ "\'" +  ",".join(variant_data['alt']) + "\'" +  ", "
					+ "\'" +  variant_data['qual'] + "\'" +  ", "
					+ "\'" +  variant_data['filter'] + "\'" +  ", "
					+ "\'" +  variant_info['AC'] + "\'" +  ", "
					+ "\'" +  variant_info['Allele'] + "\'" +  ", "
					+ "\'" +  variant_info['Consequence'] + "\'" +  ", "
					+ "\'" +  variant_info['IMPACT'] + "\'" +  ", "
					+ "\'" +  variant_info['SYMBOL'] + "\'" +  ", "
					+ "\'" +  variant_info['Gene'] + "\'" +  ", "
					+ "\'" +  variant_info['Feature'] + "\'" +  ", "
					+ "\'" +  variant_info['EXON'] + "\'" +  ", "
					+ "\'" +  variant_info['INTRON'] + "\'" +  ", "
					+ "\'" +  variant_info['STRAND'] + "\'" +  ", "
					+ "\'" +  variant_info['VARIANT_CLASS'] + "\'" +  ", "
					+ "\'" +  variant_info['SIFT'] + "\'" +  ", "
					+ "\'" +  variant_info['PolyPhen'] + "\'" +  ", "
					+ "\'" +  variant_info['gnomAD_AF'] + "\'" +  ", "
					+ "\'" +  variant_info['CLIN_SIG'] + "\'" +  ", "
					+ "\'" +  variant_info['LoF'] + "\'" +  ", "
					+ "\'" +  variant_info['LoF_filter'] + "\'" +  ", "
					+ "\'" +  variant_info['LoF_flags'] + "\'" +  ");"
				)
	cur = con.cursor()
	# print("STATEMENT is the following = " + sql_variant)
	cur.execute(sql_variant)
	con.commit()
	sample_values = variant_data['sample_values']
	for index in range(-1, len(sample_values)-1):
		geno_data = sample_values[index]
		gt_value = geno_data.get('GT', 'NA')
		dp_value = geno_data.get('DP', 'NA')
		gq_value = geno_data.get('GQ', 'NA')
		pl_value = geno_data.get('PL', 'NA')
		ad_value = geno_data.get('AD', 'NA')
		sql_genotype = ("INSERT OR IGNORE INTO " + genotypes_table +
			        "(samplename, chromosome, "
					"start, end, width, "
				    "ref, alt, "
					"gt, dp, gq, pl, ad) "
					"VALUES ("
					+ "\'" +  individuals[index] + "\'" +  ", "
					+ "\'" +  variant_data['chr'] + "\'" +  ", "
					+ "\'" +  variant_data['start'] + "\'" +  ", "
					+ "\'" +  variant_data['end'] + "\'" +  ", "
					+ "\'" +  str(variant_data['width']) + "\'" +  ", "
					+ "\'" +  variant_data['ref'] + "\'" +  ", "
					+ "\'" +  ",".join(variant_data['alt']) + "\'" +  ", "
#					+ "\'" +  geno_data['GT'] + "\'" +  ", "
#					+ "\'" +  geno_data['DP'] + "\'" +  ", "
#					+ "\'" +  geno_data['GQ'] + "\'" +  ", "
#					+ "\'" +  geno_data['PL'] + "\'" +  ");"
                        "'" + gt_value + "', " +
                        "'" + dp_value + "', " +
                        "'" + gq_value + "', " +
						"'" + pl_value + "', " +
                        "'" + ad_value + "');"
					)
		try:
			cur = con.cursor()
			# print("STATEMENT is the following = " + sql_genotype)
			#print variants causing the "sqlite3.IntegrityError: UNIQUE constraint failed: variants.chromosome, variants.start, variants.ref, variants.alt"		
			log_to_file(f"Attempting to insert: {sql_variant}")
			cur.execute(sql_variant)
			con.commit()
			log_to_file("Insertion successful or ignored if duplicate.")
		except sqlite3.IntegrityError as e:
			log_to_file(f"IntegrityError encountered: {e}")
		try:
			cur.execute(sql_genotype)
			con.commit()
		except sqlite3.IntegrityError as e_genotype:
			log_to_file(f"IntegrityError on genotype insert: {e_genotype}")

def main():
	args = get_arguments()
	if ".gz" in args.vcf:
		filein = gzip.open(args.vcf, 'rt')
	else:
		filein = open(args.vcf, 'r')
	con = create_connection(args.database)
	create_tables(con, args.variants_table, args.genotypes_table)
	for line in filein:
		if line.find("##") != -1: continue
		if "#CHROM" in line:
			fields = line.split("\t")
			individuals = fields[9:]
			continue
		## print("processing line")
		## print(line)
		## print("----------")
		variant_data = parse_variant(line)
		commit_data(con, args.variants_table, args.genotypes_table, variant_data, individuals)
	filein.close()


if __name__ == '__main__':
	main()


