import os
import csv
import time
import pysam
import random
import argparse
from multiprocessing import Pool, Lock
   
class CIGAR_OPTIONS:
	ALG_MATCH = 0
	INSERTION = 1
	DELETION = 2
	REF_SKIP = 3
	SOFT_CLIP = 4
	HARD_CLIP = 5
	PADDING = 6
	SEQ_MATCH = 7
	SEQ_MISMATCH = 8

class FLAG_OPTIONS:
	MUTLIPLE_SEGMENTS_EXIST = 0x1
	ALL_PROPERLY_ALIGNED = 0x2
	SEG_UNMAPPED = 0x4
	NEXT_SEG_UNMAPPED = 0x8
	SEQ_REV_COMP = 0x10
	NEXT_SEQ_REV_COMP = 0x20
	FIRST_SEQ_IN_TEMP = 0x40
	LAST_SEQ_IN_TEMP = 0x80
	SECONDARY_ALG = 0x100
	POOR_QUALITY = 0x200
	DUPLICATE_READ = 0x400
	SUPPLEMENTARY_ALG = 0x800

def is_msi_insertion(pattern, read_str):
	if len(read_str) <= len(pattern):
		if read_str in pattern+pattern:
			return True
		else:
			return False

	else:
		read_val = read_str
		first_pattern = read_val.find(pattern)
		if first_pattern == -1:
			return False
		if first_pattern > 0:
			sub_read = read_val[:first_pattern]
			if pattern.endswith(sub_read):
				read_val = read_val[first_pattern:]
			else:
				return False

		while read_val.startswith(pattern):
			read_val = read_val[len(pattern):]

		if pattern.startswith(read_val):
			return True
		else:
			return False

def check_msi_insertion(pattern, read_str):
	conversion_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}

	new_p=""
	for letter in pattern:
		new_p += conversion_dict[letter]

	return (is_msi_insertion(pattern, read_str) or is_msi_insertion(new_p, read_str))

def process_single_read(locus_start, locus_end, ref_num_repeats, pattern, read):
	#The value being read by pysam is off by 1 from the actual reference start point
	read_start = read.reference_start + 1
	cigar_tuples = read.cigartuples

	current_ref_pos = read_start
	current_read_pos = 0

	current_num_repeats = ref_num_repeats
	fraction_flag = 0
	overflow_deletion_flag = 0
	non_msi_insertion_flag = 0

	for cigar_op in cigar_tuples:
		if cigar_op[0] in [CIGAR_OPTIONS.HARD_CLIP, CIGAR_OPTIONS.PADDING]:
			continue;

		#consumes read but not ref
		if cigar_op[0] == CIGAR_OPTIONS.SOFT_CLIP:
			current_read_pos += cigar_op[1]
			continue;

		#consumes both read and ref
		if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
			current_ref_pos += cigar_op[1]
			current_read_pos += cigar_op[1]
			if current_ref_pos > locus_end:
				return current_num_repeats, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag
			else:
				continue;

		#consumes ref but not read
		if cigar_op[0] ==  CIGAR_OPTIONS.REF_SKIP:
			current_ref_pos += cigar_op[1]
			if current_ref_pos > locus_end:
				return current_num_repeats, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag
			else:
				continue;

		#insertion consumes read but not ref
		if cigar_op[0] == CIGAR_OPTIONS.INSERTION:
			if current_ref_pos < locus_start:
				current_read_pos += cigar_op[1]
				continue;

			if current_ref_pos >= locus_start:
				if check_msi_insertion(pattern, read.query_sequence[current_read_pos:current_read_pos+cigar_op[1]]):
					if cigar_op[1]%len(pattern) != 0:
						fraction_flag = 1
					current_num_repeats += float(cigar_op[1]/len(pattern))
					current_read_pos += cigar_op[1]
				else:
					non_msi_insertion_flag = 1
					current_read_pos += cigar_op[1]
					continue;

			if current_ref_pos > locus_end:
				return current_num_repeats, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag

		if cigar_op[0] == CIGAR_OPTIONS.DELETION:
			if current_ref_pos < locus_start:
				current_ref_pos += cigar_op[1]
				continue;


			if current_ref_pos == locus_start:
				#special case where the deletion is bigger than the whole locus
				if cigar_op[1] >= locus_end - locus_start + 1:
					return 0, 0, 1, 0
				else:
					current_ref_pos += cigar_op[1]

					if cigar_op[1]%len(pattern) != 0:
						fraction_flag = 1
					current_num_repeats -= float(cigar_op[1]/len(pattern))

					continue;

			if current_ref_pos > locus_start:
				oper = cigar_op[1]
				if current_ref_pos + oper - 1 > locus_end:
					overflow_deletion_flag = 1
					oper = locus_end - current_ref_pos + 1
				if oper%len(pattern) != 0:
					fraction_flag = 1
				current_num_repeats -= float(oper/len(pattern))

				current_ref_pos += cigar_op[1]
				continue;

			if current_ref_pos > locus_end:
				return current_num_repeats, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag

	return current_num_repeats, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag




#decide whether the read should be processed or not.
#good reads are ones that contain the entire locus, belong to a primary alignment, and are of good quality.
def should_process_read(read, locus_start, locus_end, flanking,removed_bases,p_exclude):
	#the CIGAR information isn't available
	if not read.cigartuples:
		return False

	#read starts after the start of the locus
	if read.reference_start + 1 > locus_start - flanking:
		return False

	#the read's flags indicate that it's not a good candidate for processing
	if read.flag & FLAG_OPTIONS.SECONDARY_ALG or read.flag & FLAG_OPTIONS.POOR_QUALITY \
		or read.flag & FLAG_OPTIONS.DUPLICATE_READ or read.flag & FLAG_OPTIONS.SUPPLEMENTARY_ALG:
		return False

	#the read ends before the end of the locus
	current_pos = read.reference_start
	for cigar_op in read.cigartuples:
		if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.DELETION, CIGAR_OPTIONS.REF_SKIP, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
			current_pos += cigar_op[1]
		if cigar_op[0] in [CIGAR_OPTIONS.INSERTION ]:
			removed_bases -= cigar_op[1]
		if cigar_op[0] in [CIGAR_OPTIONS.DELETION]:
			removed_bases += cigar_op[1]

	if current_pos < locus_end + flanking + removed_bases:
		return False

	#Down sampling the covere
	if random.random() < p_exclude:
		return False

	return True


#process a batch of loci and return a string with all of their results.
def process_loci(loci_list, input_bam, flanking,removed_bases,p_exclude):
	bam_file = pysam.AlignmentFile(input_bam, "rb")
	final_str_list = []

	#perform processing on each locus separately
	for locus in loci_list:
		chrom = locus[0]
		start = int(locus[3])
		end = int(locus[4])
		pattern = locus[12]
		num_repeats = float(locus[6])

		reads_iter = bam_file.fetch(contig = chrom, start = start, end = end)

		res_dict = dict()
		fraction_flag = 0
		overflow_deletion_flag = 0
		non_msi_insertion_flag = 0
		for read in reads_iter:
			if should_process_read(read, start, end, flanking,removed_bases,p_exclude):
				res, read_fraction, read_deletion_overflow, read_non_msi_insertion = process_single_read(start, end, num_repeats, pattern, read)
				if res in res_dict:
					res_dict[res]+=1
				else:
					res_dict[res]=1
				fraction_flag = fraction_flag|read_fraction
				overflow_deletion_flag = overflow_deletion_flag|read_deletion_overflow
				non_msi_insertion_flag = non_msi_insertion_flag|read_non_msi_insertion

		#process results for all reads in a single locus
		res_string = ""
		for key in sorted(res_dict.keys()):
			res_string += "{0}_{1}, ".format(key, res_dict[key])
		#remove closing comma and space
		res_string=res_string[:-2]

		final_str_for_read = "{0}:{1}:{2}:{3}:{4}, {5}\t{6}\t{7}\t{8}".format(chrom, start, end, pattern, num_repeats, res_string, fraction_flag, overflow_deletion_flag, non_msi_insertion_flag)
		final_str_list.append(final_str_for_read)

	final_str = "\n".join(final_str_list)
	return final_str


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='create INDEL information from a sample on a specific set of loci.')
	parser.add_argument("-I", "--input_BAM", help="Input BAM file.")
	parser.add_argument("-O", "--output_file", help="Output file.", default="output.txt")
	parser.add_argument("-l", "--loci_list", help="List of loci to be processed and included in the output.")
	parser.add_argument("-b", "--batch_size", help="How many loci to process in a single subprocess.", default=1000)
	parser.add_argument("-f", "--flanking", help="Length of flanking on both sides of an accepted read.", default=10)
	parser.add_argument("-r", "--removed_bases", help="Number of bases to be removed at the end of the read.", default=0)
	parser.add_argument("-e", "--exclude", help="The probability that a read will be randomly removed.", default=0)
	args = parser.parse_args()

	loci_list = csv.reader(open(args.loci_list), dialect="excel-tab")

	with Pool() as p:
		#print start time of processing to stdout
		print(time.localtime())

		#list of all the result objects returned by the subprocesses
		results = []

		#build lists of loci of <batchSize> length, and send them to the executors
		lines_count = 0
		lines_list = []
		for line in loci_list:
			if lines_count == int(args.batch_size):
				lines_list.append(line)
				results.append(p.apply_async(process_loci, args=(lines_list, args.input_BAM, int(args.flanking),int(args.removed_bases),float(args.exclude) )))
				lines_count = 0
				lines_list = []
			else:
				lines_list.append(line)
				lines_count += 1
		#handle last batch of lines
		if lines_count > 0:
			results.append(p.apply_async(process_loci, args=(lines_list, args.input_BAM, int(args.flanking),int(args.removed_bases),float(args.exclude) )))

		p.close()
		p.join()

		with open(args.output_file, "w") as output_file:
			for res in results:
				output_file.write(res.get()+"\n")
	#print end time
	print(time.localtime())
