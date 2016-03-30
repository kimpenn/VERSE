/***************************************************************

   The Subread software package is free software package: 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
#ifndef __INPUT_FILES_H_
#define __INPUT_FILES_H_

#include "HelperFunctions.h"

#define FILE_TYPE_SAM     50
#define FILE_TYPE_BAM     500
#define FILE_TYPE_FAST_   100
#define FILE_TYPE_FASTQ   105
#define FILE_TYPE_FASTA   110
#define FILE_TYPE_GZIP_FAST_   1000
#define FILE_TYPE_GZIP_FASTQ   1105
#define FILE_TYPE_GZIP_FASTA   1110
#define FILE_TYPE_UNKNOWN 999
#define FILE_TYPE_EMPTY   999990
#define FILE_TYPE_NONEXIST 999999


#define SAM_SORT_BLOCKS 229
#define SAM_SORT_BLOCK_SIZE 512333303LLU
//#define SAM_SORT_BLOCK_SIZE 11123333LLU

#define SAM_FLAG_PAIRED_TASK	0x01
#define SAM_FLAG_FIRST_READ_IN_PAIR 0x40
#define SAM_FLAG_SECOND_READ_IN_PAIR 0x80
#define SAM_FLAG_MATE_UNMATCHED 0x08
#define SAM_FLAG_MATCHED_IN_PAIR 0x02
#define SAM_FLAG_REVERSE_STRAND_MATCHED 0x10
#define SAM_FLAG_MATE_REVERSE_STRAND_MATCHED 0x20
#define SAM_FLAG_SECONDARY_MAPPING 0x100
#define SAM_FLAG_DUPLICATE 0x400
#define SAM_FLAG_UNMAPPED 0x04

#define MAX_READ_LENGTH 1210
#define MAX_READ_NAME_LEN 100
#define MAX_CHROMOSOME_NAME_LEN 100
#define MAX_FILE_NAME_LENGTH 300

typedef struct
{
	unsigned long long int output_file_size;
	unsigned long long int current_chunk_size;
	unsigned int current_chunk;
	unsigned long long int written_reads;
	unsigned long long int unpaired_reads;
	FILE * current_block_fp_array [SAM_SORT_BLOCKS];
	FILE * all_chunks_header_fp;

	FILE * out_fp;
	char tmp_path[MAX_FILE_NAME_LENGTH];
} SAM_sort_writer;


// input file
int sort_SAM_create(SAM_sort_writer * writer, char * output_file, char * tmp_path);

int is_certainly_bam_file(char * fname, int * is_firstread_PE);

int sort_SAM_add_line(SAM_sort_writer * writer, char * SAM_line, int line_len);

unsigned long long int sort_SAM_hash(char * str);

char * fgets_noempty(char * buf, int maxlen, FILE * fp);

void sort_SAM_finalise(SAM_sort_writer * writer);

void sort_SAM_check_chunk(SAM_sort_writer * writer);

int warning_file_type(char * fname, int expected_type);
#endif