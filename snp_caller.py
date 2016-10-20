################################################################
######## Julian Paris Morgan | October 16 2016
######## Homework 2 for CBB520
######## Uses BWA and Samtools to compare variants to the S288c	S.	cerevisiae	reference	sequence, available on ncbi
######## Finds private versus shared variants
######## Not meant for general pipeline usage - scripted specifically for this purpose, hence the use of shell=True and lack of input sanitization
################################################################

import subprocess
from fractions import Fraction
		
#ssh bitnami@colab-sbx-41.oit.duke.edu
#reference genome: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
#/opt/bitnami/apache2/htdocs

def generate_vcf_pipeline(ascension_number):
	#######################################
	#### Pipeline to align readsls to the reference genome and then call variants (snp's only, no indels)
	#### Input: ascension number of a strain
	#### Returns the path to the vcf file: '~/data/' + ascension_number + '.raw.vcf'
	#######################################

	reference_genome_fa = "~/data/reference_genome/GCF_000146045.2_R64_genomic.fna"
		
	#file download
	subprocess.call('fastq-dump ' + ascension_number + ' -O ~/data/temp/ --split-files', shell=True)
	file1_fastq = '~/data/temp/'+ascension_number+'_1.fastq'
	file2_fastq = '~/data/temp/'+ascension_number+'_2.fastq'
	
	print("\nfastq dump completed")
	print(ascension_number)	
	#bwa mem
	subprocess.call('bwa mem ' + reference_genome_fa + ' ' + file1_fastq + ' ' + file2_fastq + ' > ~/data/temp/alignment_' + ascension_number + '.sam', shell = True)
	alignment_sam = '~/data/temp/alignment_' + ascension_number + '.sam'
	print("\n bwa mem completed")

	#create bam files (binary sam file).
	subprocess.call('samtools view -S ' + alignment_sam + ' -b -o ~/data/temp/alignment_' + ascension_number + '.bam', shell = True) 
	alignment_bam = '~/data/temp/alignment_' + ascension_number + '.bam'
	print('\n samtools convert sam -> bam completed')
	
	#sort bam file
	subprocess.call('samtools sort ' + alignment_bam + ' ~/data/temp/alignment_' + ascension_number+'_sorted', shell=True)
	alignment_bam_sorted = '~/data/temp/alignment_' + ascension_number + '_sorted.bam'
	print('\n samtools sort completed')

	#create index of the bam file
	subprocess.call('samtools index ' + alignment_bam_sorted, shell=True)
	print('\n samtools index completed\n')
	
	# #call variants and write to a vcf file		
	subprocess.call('samtools mpileup -I -uf ' +  reference_genome_fa + ' ' + alignment_bam_sorted + ' | bcftools view -vcg - > ~/data/' + ascension_number + '.raw.vcf', shell = True)
	print('\n samtools mpileup completed')
	
	#quality control step - cutoff at 30	
	#delete all files except the .vcf file
	subprocess.call('rm ~/data/temp/*', shell = True)	
	print('deleted all extra files')

	return '~/data/' + ascension_number + '.raw.vcf'

def compress_index(ascension_number):
	#######################################
	#### Prepares a file to be used by vfc tools
	#### 	-Compresses file (.gz) and indexes it
	#### Input: ascension number of a strain
	#### Returns the path to the .gz file: '~/data/' + ascension_number + '.raw.vcf.gz'
	#######################################
	
	#compress by bgzip
	subprocess.call('bgzip ~/data/' + ascension_number + '.raw.vcf', shell=True)
	#index with tabix
	subprocess.call('tabix ~/data/' + ascension_number + '.raw.vcf.gz', shell=True)

	return '~/data/' + ascension_number + '.raw.vcf.gz'

def compare_snps(strain_paths_gz):
	#######################################
	#### Prepares a file to be used by vfc tools
	#### 	-Compresses file (.gz) and indexes it
	#### Input: A list of paths to .gz files that have been indexed
	#### Returns the path to the .gz file: '~/data/' + ascension_number + '.raw.vcf.gz'
	#######################################
	file_path = ' '.join(strain_paths_gz)
	subprocess.call('vcf-isec -f -n +2 ' + file_path+ ' | bgzip -c > ~/data/common_snp_file.vcf.gz', shell=True)
	return '~/data/common_snp_file.vcf.gz'

def generate_summary_stats(strains):
	#######################################
	#### Generate summary statistics on the number of private vs. common snp's for each strain. 
	#### Subdivided based on chromosome
	#### Input: None Output: None
	#######################################

	#return file: formatted html tables
	output_file = open('/home/bitnami/data/snp_table.html', 'w')

	#pre-processing step with the common snp file
	
	#unzip the common file
	subprocess.call('gunzip ~/data/common_snp_file.vcf.gz', shell = True)
	
	#generate a dictionary containing the common SNP's
	#format: common_snps = { chrom1= [snp_location1, snp_location2, snp_location2, etc. ], chrom2=[etc], etc }
	common_snp_file = open('/home/bitnami/data/common_snp_file.vcf', 'r')
	common_snps = {}
	for line in common_snp_file:
		#skip if starts with '#'
		if line[0] == '#':
			continue
		ln = line.strip().split('\t')
		cur_chrom = ln[0]
		if cur_chrom not in common_snps:
			common_snps[cur_chrom] = [ln[1]]
		else:
			temp = common_snps[cur_chrom]
			temp.append(ln[1])
			common_snps[cur_chrom] = temp
	
	for key, value in strains.iteritems():
		strain_name = key
		ascension_number = value

		#unzip the file
		subprocess.call('gunzip /home/bitnami/data/' + ascension_number + '.raw.vcf.gz', shell = True)
		
		#set up the html table
		table = '<br>' + strain_name +  ' : ' + ascension_number + "<table>  <tr> <th>Chromosome</th> <th>Private Number</th> <th>Private Density</th> <th>Shared Number</th> <th>Shared Density</th>  </tr>"

		#with open('~/data/' + ascension_number + '.raw.vcf', 'rb') as f:
		with open('/home/bitnami/data/'+ascension_number + '.raw.vcf', 'rb') as f:
			#check for correct open?
			chrom = 'NC_001133.9'
			private_snp_cnt = 0
			public_snp_cnt = 0
			chrom_length = 0
			for line in f:
				#skip if starts with '#'
				if line[0] == '#':
					continue
				ln = line.strip().split('\t')
				cur_chrom = ln[0]
				#check for current or new chromosome
				if cur_chrom != chrom:
					#store data
					chrom = int(chrom[5:9]) - 1132
					table += "<tr> <th>"+str(chrom)+"</th> <th>"+str(private_snp_cnt)+"</th><th>"+str(Fraction(private_snp_cnt,chrom_length))+"</th> <th>"+str(public_snp_cnt)+"</th> <th>"+str(Fraction(public_snp_cnt,chrom_length))+"</th>  </tr>"
					#reset snp_counts
					chrom = cur_chrom
					private_snp_cnt = 0
					public_snp_cnt = 0
					chrom_length = 0
				
				#check if common or private snp
				if ln[1] in common_snps[cur_chrom]:
					public_snp_cnt += 1
				else:
					private_snp_cnt += 1
				
				#increment chromosome length
				chrom_length += 1
			table += "<tr> <th> Mitochondria </th> <th>"+str(private_snp_cnt)+"</th><th>"+str(Fraction(private_snp_cnt,chrom_length))+"</th> <th>"+str(public_snp_cnt)+"</th> <th>"+str(Fraction(public_snp_cnt,chrom_length))+"</th>  </tr>"
	
		table += '</table>'
		output_file.write(table)
	
	common_snp_file.close()		
	return

if __name__ == "__main__":
	
	strains = {
		'yjm1190': 'SRR800799',
		'yjm1439': 'SRR800837',
		'yjm689': 'SRR800786',
		'yjm193': 'SRR800765',
		'yjm975': 'SRR800790',
		'yjm969': 'SRR800788',
		'yjm984': 'SRR800793',
		'yjm1199': 'SRR800800',	
		}

	#Preprocessing step: make index of reference genome
		subproceess.call('bwa index reference_genome.fasta', shell=True)
		subproceess.call('samtools faidx reference_genome.fasta', shell=True)

	#stores file paths of all vcf.gz files
	gz_strain_paths = []	
	for key, value in strains.iteritems():
		strain_name = key
		ascension_number = value
		
		#aligns strain, calls variants
		generate_vcf_pipeline(ascension_number)
		
		#prepare file for vcftools to use
		gz_strain_path = compress_index(ascension_number)

		gz_strain_paths.append(gz_strain_path)

	gz_strain_paths = ['~/data/SRR800793.raw.vcf.gz', '~/data/SRR800788.raw.vcf.gz', '~/data/SRR800800.raw.vcf.gz', '~/data/SRR800790.raw.vcf.gz', '~/data/SRR800837.raw.vcf.gz', '~/data/SRR800786.raw.vcf.gz', '~/data/SRR800799.raw.vcf.gz', '~/data/SRR800765.raw.vcf.gz']
		
	#find common snps
	common_snp_path = compare_snps(gz_strain_paths)

	#calculate summary statistics for all strains
	generate_summary_stats(strains)
	
