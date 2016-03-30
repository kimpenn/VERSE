/***************************************************************
 
    VERSE is developed based on the framework of featureCounts
    (SUBREAD).
    
    Current version supports the following modes of RNA-Seq
    quantification:
    1. featureCounts (Default)
    2. HTSeq Union (-z 1)
    3. HTSeq Intersection-strict (-z 2)
    4. HTSeq Intersection-nonempty (-z 3)
    5. VERSE Union-strict (-z 4)
    6. VERSE Cover-length (-z 5)
    
    VERSE also allows hierarchical-assign and independent-assign 
    for multiple feature types.
 
    For more information, please refer to the user manual.
 
    This work is supervised and generously supported by 
    Dr. Stephen Fisher and Professor Junhyong Kim.
 
    Qin Zhu
    Junhyong Kim Lab
    University of Pennsylvania
    2015
 
  ***************************************************************/

#define _GNU_SOURCE
#include <assert.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "sambam-file.h"
#include "HelperFunctions.h"
#include "input-files.h"

#define FEATURE_NAME_LENGTH  256 
#define CHROMOSOME_NAME_LENGTH 256 
#define MAX_LINE_LENGTH 3000
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

#define ALLOW_ALL_MULTI_MAPPING 1
#define ALLOW_PRIMARY_MAPPING 2

#define SECTION_LEVEL 0
#define READ_LEVEL 1

#define MAX_HIT_NUMBER 1000
#define MAX_CIGAR_SECTIONS 6
#define MAX_FEATURE_TYPE_NUM 12

#define BLOCK_CHUNK_NUM 128
#define MAX_UNPROCESSED_NUMBER 10000


typedef struct
{
	unsigned int feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;
	unsigned short chro_name_pos_delta;
	char is_negative_strand;
} fc_feature_info_t;

typedef struct
{
    char * read_name;
    int alignment_masks;
    char * CIGAR_str;
    char * read_chr;
    long read_pos;
    int allow_process;
} read_info_t;

typedef struct
{
	unsigned long long assigned_reads;
	unsigned long long unassigned_ambiguous;
	unsigned long long unassigned_multimapping;
	unsigned long long unassigned_nofeatures;
	unsigned long long unassigned_unmapped;
	unsigned long long unassigned_mappingquality;
	unsigned long long unassigned_fragmentlength;
	unsigned long long unassigned_chimericreads;
	unsigned long long unassigned_secondary;
	unsigned long long unassigned_nonjunction;
	unsigned long long unassigned_duplicate;
} fc_read_counters;

typedef struct
{
	unsigned short thread_id;
	char * line_buffer1;
	char * line_buffer2;
	unsigned long long int all_reads;
	unsigned int ** count_table;
	unsigned int chunk_read_ptr;
	pthread_t thread_object;

	char * input_buffer;
	unsigned int input_buffer_remainder;
	unsigned int input_buffer_write_ptr;	
	pthread_mutex_t input_buffer_lock;
    char step_back;

	short hits_total_length1[MAX_HIT_NUMBER];
	short hits_total_length2[MAX_HIT_NUMBER];
	long hits_indices1 [MAX_HIT_NUMBER];
	long hits_indices2 [MAX_HIT_NUMBER];
	long decision_table_ids [MAX_HIT_NUMBER];
	unsigned char decision_table_votes [MAX_HIT_NUMBER];
	long decision_table_exon_ids [MAX_HIT_NUMBER];
	long uniq_gene_exonid_table [MAX_HIT_NUMBER];
	long uniq_gene_table [MAX_HIT_NUMBER];
    long sec_gene_table [MAX_HIT_NUMBER];
    long sec_gene_exonid_table [MAX_HIT_NUMBER];
    long sec_table_ids [MAX_HIT_NUMBER];
    unsigned char sec_table_votes [MAX_HIT_NUMBER];
    long sec_table_exon_ids [MAX_HIT_NUMBER];

	char * chro_name_buff;
	z_stream * strm_buffer;

	fc_read_counters * read_counters;

	SamBam_Alignment aln_buffer;
    
    char ** unprocessed_read_ptrs;
    long unprocessed_read_cnts;
    unsigned long long int missing_mates;
} fc_thread_thread_context_t;

#define REVERSE_TABLE_BUCKET_LENGTH 131072 // 128kb
#define REDUCE_TO_5_PRIME_END 5
#define REDUCE_TO_3_PRIME_END 3

typedef struct
{
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

typedef struct
{
	int is_paired_end_data;
	int is_strand_checked;
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int is_multi_mapping_allowed;
	int is_input_file_resort_needed;
	int is_SAM_file;
	int is_read_details_out;
	int is_unpaired_warning_shown;
	int is_split_alignments_only;
	int is_duplicate_ignored;
	int reduce_5_3_ends_to_one;
	int isCVersion;

	int min_mapping_quality_score;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int feature_block_size;
	int read_length;
	int line_length;
	int * longest_chro_name;
	int five_end_extension;
	int three_end_extension;
    int read_M_len_required;
	int overlap_length_required;
    int nonoverlap_length_allowed;
    int min_dif_ambiguous;
    
    int run_mode;
    int is_independent_assign;
    int is_nonempty_modified;
    int multithread_unzipping;

	unsigned long long int all_reads;
    unsigned long long int missing_mates;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;        // this marks the completion of the bulk assignment
    int is_really_finished;     // if paired-end, this marks the completion of the assignment of the unprocessed reads
	unsigned int input_buffer_max_size;
	SamBam_Reference_Info * sambam_chro_table;

	char * debug_command;
	char ** unistr_buffer_space;
	unsigned int * unistr_buffer_size;
	unsigned int * unistr_buffer_used;

	HashTable ** gene_name_table;
	char input_file_name[300];
	char raw_input_file_name[300];
	char output_file_name[300];
	unsigned char *** gene_name_array;

	HashTable ** exontable_chro_table;
	int * exontable_nchrs;
	int * exontable_exons;
    int * original_exons;
	int ** exontable_geneid;
	char ** exontable_strand;
	long ** exontable_start;
	long ** exontable_stop;
    
	char feature_name_column[100];
    char ** feature_type_list;
    int feature_type_num;
    
	char gene_id_column[100];

	long ** exontable_block_end_index;
	long ** exontable_block_max_end;
	long ** exontable_block_min_start;

	FILE * detail_output_fp;
	double start_time;

	char * cmd_rebuilt;
    
    char ** unprocessed_reads;
    long long int unprocessed_cnts;
    
    long long int total_assigned;

	fc_read_counters * read_counters;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;

unsigned int unistr_cpy(fc_thread_global_context_t * global_context, char * str, int strl, int f_idx)
{
    unsigned int ret;
    if(global_context->unistr_buffer_used[f_idx] + strl >= global_context->unistr_buffer_size[f_idx]-1)
    {
        if( global_context->unistr_buffer_size[f_idx] < 3435973835u) // 4G / 5 * 4 - 5
        {
            global_context -> unistr_buffer_size[f_idx] = global_context->unistr_buffer_size[f_idx] /4 *5;
            global_context -> unistr_buffer_space[f_idx] = realloc(global_context -> unistr_buffer_space[f_idx], global_context->unistr_buffer_size[f_idx]);
        }
        else
        {
            SUBREADprintf("Error: exceed memory limit (4GB) for storing annotation data.\n");
            return 0xffffffffu;
        }
    }
    strcpy(global_context -> unistr_buffer_space[f_idx] + global_context->unistr_buffer_used[f_idx], str);
    ret = global_context->unistr_buffer_used[f_idx];
    global_context->unistr_buffer_used[f_idx] += strl +1;
    return ret;
}

void print_FC_configuration(fc_thread_global_context_t * global_context, char * annot, char * sam, char * out, int is_sam, int *n_input_files, int isReadDetailReport, char* feature_type)
{
	char * tmp_ptr1 = NULL , * next_fn, *sam_used = malloc(strlen(sam)+1), sam_ntxt[30],bam_ntxt[30], next_ntxt[50];
	int nfiles=1, nBAMfiles = 0, nNonExistFiles = 0;

	strcpy(sam_used, sam);

	VERSEputs("");
	print_verse_logo();
	VERSEputs("");
	print_in_box(80,1,1,"VERSE setting");
    switch(global_context->run_mode){
        case 0:print_in_box(80,0,0,"           Running mode : Default(featureCounts)"); break;
        case 1:print_in_box(80,0,0,"           Running mode : HTSeq-Union"); break;
        case 2:print_in_box(80,0,0,"           Running mode : HTSeq-Intersection_strict"); break;
        case 3:{
            if(!global_context->is_nonempty_modified)
                print_in_box(80,0,0,"           Running mode : HTSeq-Intersection_nonempty");
            else
                print_in_box(80,0,0,"           Running mode : HTSeq-Intersection_nonempty_modified");
            break;
        }
        case 4:print_in_box(80,0,0,"           Running mode : Union_strict"); break;
        case 5:print_in_box(80,0,0,"           Running mode : Cover_length"); break;
        default: print_in_box(80,0,0,"           Running mode : Default(featureCounts)");
    }
    
    print_in_box(80,0,0,"           Feature type : %s", feature_type);
    if(global_context -> feature_type_num > 1)
    {
        if(!global_context -> is_independent_assign)
            print_in_box(80,0,0,"Hierarchically assigning : %i feature types", global_context -> feature_type_num);
        else
            print_in_box(80,0,0,"Independently assigning : %i feature types", global_context -> feature_type_num);
    }
	nfiles = 0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		nfiles++;

		int file_probe = is_certainly_bam_file(next_fn, NULL);
		if(file_probe==-1) nNonExistFiles++;
		if(file_probe == 1) nBAMfiles++;		
	}

	sam_ntxt[0]=0;
	bam_ntxt[0]=0;
	next_ntxt[0]=0;

	if(nNonExistFiles)
		sprintf(next_ntxt, "%d unknown file%s", nNonExistFiles, nNonExistFiles>1?"s":"");
	if(nBAMfiles)
		sprintf(bam_ntxt, "%d BAM file%s  ", nBAMfiles, nBAMfiles>1?"s":"");
	if(nfiles-nNonExistFiles-nBAMfiles)
		sprintf(sam_ntxt, "%d SAM file%s  ", nfiles-nNonExistFiles-nBAMfiles , (nfiles-nNonExistFiles-nBAMfiles)>1?"s":"");

	strcpy(sam_used, sam);
	print_in_box(80,0,0,"             Input file : %s%s%s", sam_ntxt, bam_ntxt, next_ntxt);
	nfiles=0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		int is_first_read_PE = 0 , file_probe = is_certainly_bam_file(next_fn, &is_first_read_PE);

		char file_chr = 'S';
		if(file_probe == -1) file_chr = '?';
		else if(is_first_read_PE == 1) file_chr = 'P';
		//file_chr = 'o';

		print_in_box(94,0,0,"                          %c[32m%c%c[36m %s%c[0m",CHAR_ESC, file_chr,CHAR_ESC, next_fn,CHAR_ESC);
		nfiles++;
	}

	(*n_input_files) = nfiles;
	print_in_box(80,0,0,"            Output file : %s", out);
	print_in_box(80,0,0,"        Annotation file : %s", annot);
	if(isReadDetailReport)
		print_in_box(80,0,0,"     Assignment details : %s.detail.txt", out);

	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"                Threads : %d", global_context->thread_number);
    
    if(nBAMfiles && global_context->thread_number > 1)
        print_in_box(80,0,0,"Multithreaded BAM unzip : %s", global_context->multithread_unzipping?"yes":"no");
    
    print_in_box(80,0,0,"");
	print_in_box(80,0,0,"             Paired-end : %s", global_context->is_paired_end_data?"yes":"no");
	print_in_box(80,0,0,"        Strand specific : %s", global_context->is_strand_checked?(global_context->is_strand_checked==1?"yes":"inversed"):"no");
	char * multi_mapping_allow_mode = "not counted";
	if(global_context->is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING)
		multi_mapping_allow_mode = "primary only";
	else if(global_context->is_multi_mapping_allowed == ALLOW_ALL_MULTI_MAPPING)
		multi_mapping_allow_mode = "counted";
	print_in_box(80,0,0,"     Multimapping reads : %s", multi_mapping_allow_mode);
	if(global_context -> is_split_alignments_only)
		print_in_box(80,0,0,"       Split alignments : required");
    if(global_context -> read_M_len_required !=1)
    {
        print_in_box(80,0,0,"   Required Read Length : %d", global_context -> read_M_len_required);
    }
    if(global_context -> overlap_length_required !=1)
        print_in_box(80,0,0,"  Min Overlapping bases : %d", global_context -> overlap_length_required);
    if(global_context -> nonoverlap_length_allowed != 0x7fff)
    {
        print_in_box(80,0,0,"Max Nonoverlapping bases: %d", global_context -> nonoverlap_length_allowed);
        if(global_context -> run_mode == 2 || global_context -> run_mode == 4)
        {
            print_in_box(80,0,0,"      NOTE You are specifying maxReadNonoverlap in strict mode,");
            print_in_box(80,0,0,"      This option will have no effect.");
        }
    }
	if(global_context -> five_end_extension || global_context -> three_end_extension)
		print_in_box(80,0,0,"        Read extensions : %d on 5' and %d on 3' ends", global_context -> five_end_extension , global_context -> three_end_extension);
	if(global_context -> reduce_5_3_ends_to_one)
		print_in_box(80,0,0,"      Read reduction to : %d' end" , global_context -> reduce_5_3_ends_to_one == REDUCE_TO_5_PRIME_END ?5:3);
	if(global_context -> is_duplicate_ignored)
		print_in_box(80,0,0,"       Duplicated Reads : ignored");
    
    if(global_context -> min_dif_ambiguous)
    {
        print_in_box(80,0,0,"Min difference required : %d", global_context -> min_dif_ambiguous);
    }

	if(global_context->is_paired_end_data)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"         Chimeric reads : %s", global_context->is_chimertc_disallowed?"not counted":"counted");
		print_in_box(80,0,0,"       Both ends mapped : %s", global_context->is_both_end_required?"required":"not required");

		if(global_context->is_PE_distance_checked)
			print_in_box(80,0,0,"         TLEN range : %d - %d", global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
	}

	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"Based on the framework of featureCounts(SUBREAD)");
	VERSEputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,0,"");

	free(sam_used);
}

void print_FC_results(fc_thread_global_context_t * global_context)
{
	print_in_box(89,0,1,"%c[36mRead assignment finished.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"VERSE: a Versatile and Efficient Rna-SEq read assignment tool");
	VERSEputs("");
	return;
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}

int is_comment_line(const char * l, int file_type, unsigned int lineno)
{
	int tabs = 0, xk1 = 0;
	if(l[0]=='#') return 1;

	if(isalpha(l[0]) && file_type == FILE_TYPE_RSUBREAD)
	{
		char target_chr[16];
		memcpy(target_chr, l, 16);
		for(xk1=0; xk1<16; xk1++)
			target_chr[xk1] = tolower(target_chr[xk1]);

		if(memcmp(target_chr, "geneid\tchr\tstart",16)==0) return 1;
	}

	xk1=0;
	while(l[xk1]) tabs += (l[xk1++] == '\t');

	return tabs < ((file_type == FILE_TYPE_GTF)?8:4);
}

// This function loads annotations from the file.
// It returns the number of features loaded, or -1 if something is wrong.
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, fc_feature_info_t ** loaded_features, int f_idx)
{
    unsigned int features = 0, xk1 = 0, lineno=0;
    char * file_line = malloc(MAX_LINE_LENGTH+1);
    FILE * fp = fopen(annotation_file,"r");
    int is_GFF_warned = 0;
    if(!fp) return -1;
    
    HashTable * chro_name_table = HashTableCreate(1603);
    HashTableSetHashFunction(chro_name_table, fc_chro_hash);
    HashTableSetKeyComparisonFunction(chro_name_table, fc_strcmp_chro);
    global_context -> longest_chro_name[f_idx] = 0;
    
    // first scan: get the chromosome size, etc
    while(1)
    {
        char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
        char * token_temp, *chro_name;
        fc_chromosome_index_info * chro_stab;
        unsigned int feature_pos = 0;
        if(!fgets_ret) break;
        
        lineno++;
        if(is_comment_line(file_line, FILE_TYPE_GTF, lineno-1))continue;
        
        chro_name = strtok_r(file_line,"\t",&token_temp);
        strtok_r(NULL,"\t", &token_temp); // lib_name (not needed)
        char * feature_type = strtok_r(NULL,"\t", &token_temp);
        if(strcmp(feature_type, global_context -> feature_type_list[f_idx])==0)
        {
            strtok_r(NULL,"\t", &token_temp); // feature_start
            feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end
            features++;
        }
        else chro_name = NULL;
        
        if(chro_name)
        {
            if(strlen(chro_name)>=CHROMOSOME_NAME_LENGTH)
                chro_name[CHROMOSOME_NAME_LENGTH-1]=0;
            chro_stab = HashTableGet(chro_name_table, chro_name);
            
            if(chro_stab)
            {
                chro_stab -> chro_possible_length = max(chro_stab -> chro_possible_length , feature_pos+1);
            }else
            {
                char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
                term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
                chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
                chro_stab -> chro_number = chro_name_table->numOfElements;
                chro_stab -> chro_possible_length = feature_pos+1;
                chro_stab -> reverse_table_start_index = NULL;
                HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
            }
            
            chro_stab -> chro_features ++;
            //printf("chro_name = %s, features = %d, chro_features = %d, chro_possible_length = %d\n", chro_name, features,chro_stab -> chro_features, chro_stab ->chro_possible_length);
        }
        
    }
    
    fseek(fp,0,SEEK_SET);
    
    fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);
    
    lineno = 0;
    while(xk1 < features)
    {
        int is_gene_id_found = 0;
        fgets(file_line, MAX_LINE_LENGTH, fp);
        lineno++;
        char * token_temp;
        if(is_comment_line(file_line, FILE_TYPE_GTF, lineno-1))continue;
        
        
        char feature_name_tmp[FEATURE_NAME_LENGTH];
        sprintf(feature_name_tmp, "LINE_%07u", xk1 + 1);
        char * seq_name = strtok_r(file_line,"\t",&token_temp);
        strtok_r(NULL,"\t", &token_temp);       //source
        char * feature_type = strtok_r(NULL,"\t", &token_temp);     //feature_type
        if(strcmp(feature_type, global_context -> feature_type_list[f_idx])==0)
        {
            ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));       //start
            ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));         //end
            
            if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 ||  ret_features[xk1].start > 0x7fffffff ||  ret_features[xk1].end > 0x7fffffff || ret_features[xk1].start > ret_features[xk1].end)
                SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);
            
            
            strtok_r(NULL,"\t", &token_temp);   //score
            ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));      //strand
            ret_features[xk1].sorted_order = xk1;
            strtok_r(NULL,"\t",&token_temp);	// "frame"
            char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ...
            if(extra_attrs && (strlen(extra_attrs)>2))
            {
                int attr_val_len = GTF_extra_column_value(extra_attrs , global_context -> gene_id_column, feature_name_tmp, FEATURE_NAME_LENGTH);
                if(attr_val_len>0) is_gene_id_found=1;
                //		printf("V=%s\tR=%d\n", extra_attrs , attr_val_len);
            }
            
            if(is_gene_id_found)
            {
            }
            else
            {
                if(!is_GFF_warned)
                {
                    int ext_att_len = strlen(extra_attrs);
                    if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
                    SUBREADprintf("\nWarning: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nThe attributes included in your GTF annotation are '%s' \n\n",  global_context -> gene_id_column, extra_attrs);
                }
                is_GFF_warned++;
            }
            
            int feature_name_len = strlen(feature_name_tmp);
            if(feature_name_len > FEATURE_NAME_LENGTH) feature_name_tmp[FEATURE_NAME_LENGTH - 1] = 0;
            ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name_tmp, feature_name_len, f_idx);
            int chro_name_len = strlen(seq_name);
            if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
            unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len, f_idx);
            global_context -> longest_chro_name[f_idx] = max(chro_name_len, global_context -> longest_chro_name[f_idx]);
            
            ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;
            
            int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
            fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
            if(!chro_stab -> reverse_table_start_index)
            {
                chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
                memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
            }
            chro_stab -> reverse_table_start_index[bin_location]++;
            xk1++;
        }
    }
    fclose(fp);
    free(file_line);
    
    (*loaded_features) = ret_features;
    global_context -> exontable_nchrs[f_idx] = (int)chro_name_table-> numOfElements;
    
    global_context -> exontable_chro_table[f_idx] = chro_name_table;
    
    print_in_box(80,0,0,"   %ss : %d\n", global_context -> feature_type_list[f_idx], features);
    if(features < 1)
        SUBREADprintf("WARNING no features were loaded.\n");
    return features;
}

int find_or_insert_gene_name(fc_thread_global_context_t * global_context, unsigned char * feature_name, int f_idx)
{
	HashTable * genetable = global_context -> gene_name_table[f_idx];

	long long int gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else
	{
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		global_context -> gene_name_array[f_idx][gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

// This function puts into the reverse_table_start_index[#bin] the first block id in that bin.
void register_reverse_table(int block_no, long this_block_min_start, long this_block_max_end, fc_chromosome_index_info * chro_inf)
{

	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH; // start #bin of the block
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH; // end #bin of the block
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++) // from the start bin to the end bin
	{
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
		//chro_inf->reverse_table_end_index[x1] = max(chro_inf->reverse_table_end_index[x1], block_no+1);
	}

}

// These sort functions are for sorting unprocessed reads or features, using provided prameters (index). All these functions are passed to the merge_sort(...) function.
// Index = 0, sort by start; index = 2, sort by strand; index = 3, sort by gene entrez
void read_sort_merge(void * arrv, int start, int items, int items2, int index)
{
    
    void ** arr = (void **) arrv;
    char ** read_names = arr[0];
    char ** read_lines = arr[1];
    
    int total_items = items+items2;
    char ** tmp_names = malloc(sizeof(char *) * total_items);
    char ** tmp_lines = malloc(sizeof(char *) * total_items);
    
    int read_1_ptr = start;
    int read_2_ptr = start+items;
    int write_ptr;
    
    for(write_ptr=0; write_ptr<total_items; write_ptr++)
    {
        if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && strcmp(read_names[read_1_ptr], read_names[read_2_ptr]) >= 0))
        {
            tmp_names[write_ptr] = read_names[read_2_ptr];
            tmp_lines[write_ptr] = read_lines[read_2_ptr];
            read_2_ptr++;
        }
        else
        {
            tmp_names[write_ptr] = read_names[read_1_ptr];
            tmp_lines[write_ptr] = read_lines[read_1_ptr];
            read_1_ptr++;
        }
    }
    
    memcpy(read_names + start, tmp_names, sizeof(char *) * total_items);
    memcpy(read_lines + start, tmp_lines, sizeof(char *) * total_items);
    
    free(tmp_names);
    free(tmp_lines);
}

int read_sort_compare(void * arrv, int l, int r, int index)
{
    void ** arr = (void **) arrv;
    char ** read_names = arr[0];
    char * ll = read_names[l];
    char * rl = read_names[r];
    return strcmp(ll, rl);
}

void read_sort_exchange(void * arrv, int l, int r)
{
    void ** arr = (void **) arrv;
    char * tmp;
    char ** read_names = arr[0];
    char ** read_lines = arr[1];
    
    tmp = read_names[r];
    read_names[r]=read_names[l];
    read_names[l]=tmp;
    
    tmp = read_lines[r];
    read_lines[r]=read_lines[l];
    read_lines[l]=tmp;
}

void feature_sort_merge(void * arrv, int start, int items, int items2, int index)
{
    
    void ** arr = (void **) arrv;
    
    long * ret_start = (long *) arr[0];
    long * ret_end = (long *) arr[1];
    unsigned char * ret_strand = (unsigned char *) arr[2];
    int * ret_entyrez = (int *) arr[3];
    fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];
    
    int total_items = items+items2;
    long * tmp_start = malloc(sizeof(long) * total_items);
    long * tmp_end = malloc(sizeof(long) * total_items);
    unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
    int * tmp_entyrez = malloc(sizeof(int) * total_items);
    fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);
    
    int read_1_ptr = start;
    int read_2_ptr = start+items;
    int write_ptr;
    
    for(write_ptr=0; write_ptr<total_items; write_ptr++)
    {
        if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ((index == 0 && ret_start[read_1_ptr] >= ret_start[read_2_ptr])||(index == 2 && ret_strand[read_1_ptr] >= ret_strand[read_2_ptr]) || (index == 3 && ret_entyrez[read_1_ptr] >= ret_entyrez[read_2_ptr]))))
        {
            tmp_start[write_ptr] = ret_start[read_2_ptr];
            tmp_end[write_ptr] = ret_end[read_2_ptr];
            tmp_strand[write_ptr] = ret_strand[read_2_ptr];
            tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
            tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
            read_2_ptr++;
        }
        else
        {
            tmp_start[write_ptr] = ret_start[read_1_ptr];
            tmp_end[write_ptr] = ret_end[read_1_ptr];
            tmp_strand[write_ptr] = ret_strand[read_1_ptr];
            tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
            tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
            read_1_ptr++;
        }
    }
    
    memcpy(ret_start+ start, tmp_start, sizeof(long) * total_items);
    memcpy(ret_end+ start, tmp_end, sizeof(long) * total_items);
    memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
    memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
    memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);
    
    free(tmp_start);
    free(tmp_end);
    free(tmp_strand);
    free(tmp_entyrez);
    free(tmp_info_ptr);
}

int feature_sort_compare(void * arrv, int l, int r, int index)
{
	void ** arr = (void **) arrv;
    long * ret_start = (long *)arr[0];
    unsigned char * ret_strand = (unsigned char *) arr[2];
    int * ret_entrez = (int *)arr[3];
    long ll = 0, rl = 0;
    switch(index)
    {
        case 0:
            ll = ret_start[l];
            rl = ret_start[r];
            break;
        case 2:
            ll = (long)(ret_strand[l]);
            rl = (long)(ret_strand[r]);
            break;
        case 3:
            ll = (long)(ret_entrez[l]);
            rl = (long)(ret_entrez[r]);
            break;
        default:
            printf("Invalid comparison index, something's wrong!\n");
            break;
    }
    
	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long tmp;
	fc_feature_info_t * tmpptr;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}

// This function merges overlapping exons that belong to the same gene.
void verse_feature_merge(fc_thread_global_context_t * global_context, void * arrv, void * arrv2, int * new, fc_chromosome_index_info * this_chro_inf)
{
    void ** old_arr = (void **) arrv;
    long * ret_start = (long *) old_arr[0];
    long * ret_end = (long *) old_arr[1];
    unsigned char * ret_strand = (unsigned char *) old_arr[2];
    int * ret_entrez = (int *) old_arr[3];
    fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) old_arr[4];
    
    void ** chro_arr = (void **) arrv2;
    long * chro_start = (long *) chro_arr[0];
    long * chro_end = (long *) chro_arr[1];
    unsigned char * chro_strand = (unsigned char *) chro_arr[2];
    int * chro_entrez = (int *) chro_arr[3];
    fc_feature_info_t ** chro_feature_ptr = (fc_feature_info_t **) chro_arr[4];
    
    int xk1, cnt = 0;
    
    for(xk1 = 0; xk1 < this_chro_inf -> chro_features; xk1++)
    {
        if (cnt != 0 && ret_entrez[xk1] == chro_entrez[cnt-1] && (ret_start[xk1] <= chro_end[cnt-1] + 1) && (!global_context -> is_strand_checked || (global_context -> is_strand_checked && ret_strand[xk1] == chro_strand[cnt-1])))     // overlap!
        {
            if (ret_end[xk1] > chro_end[cnt-1])         // new end!
            {
                chro_end[cnt-1] = ret_end[xk1];
                chro_feature_ptr[cnt-1] -> end = ret_end[xk1];
            }
            else                                        // skip
                continue;
        }
        else
        {
            chro_start[cnt] = ret_start[xk1];
            chro_end[cnt] = ret_end[xk1];
            chro_strand[cnt] = ret_strand[xk1];
            chro_entrez[cnt] = ret_entrez[xk1];
            chro_feature_ptr[cnt++] = old_info_ptr[xk1];
            (*new)++;
        }
        //printf("chro = %i\toriginal_feature_cnts=%i\txk1=%i\tcnt=%i\tnew=%i\tret_start=%li\tchro_start=%li\tret_end=%li\tchro_end=%li\tret_strand=%i\tchro_strand=%i\tret_entrez=%i\tchro_entrez=%i\n", this_chro_inf -> chro_number, this_chro_inf -> chro_features, xk1, cnt, (*new) - 1, ret_start[xk1], chro_start[cnt-1], ret_end[xk1], chro_end[cnt-1], ret_strand[xk1], chro_strand[cnt-1], ret_entrez[xk1], chro_entrez[cnt-1]);
    }
    this_chro_inf -> chro_features = cnt;
}

// Process annotation information, including mergeing exons and sorting by start position, etc.
int sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, long ** block_end_index, long ** block_min_start_pos, long ** block_max_end_pos, fc_feature_info_t *** merged_features_ptr, int f_idx)
{
    unsigned int chro_pnt;
    unsigned int xk1,xk2;
    int * ret_entrez = malloc(sizeof(int) * features);
    long * ret_start = malloc(sizeof(long) * features);
    long * ret_end = malloc(sizeof(long) * features);
    int current_block_buffer_size = 2000;
    
    long * ret_block_end_index = malloc(sizeof(long) * current_block_buffer_size);
    long * ret_block_min_start = malloc(sizeof(long) * current_block_buffer_size);
    long * ret_block_max_end = malloc(sizeof(long) * current_block_buffer_size);
    unsigned char * ret_strand = malloc(features);
    fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);

    unsigned int * chro_feature_ptr = calloc(sizeof(int) * global_context -> exontable_nchrs[f_idx],1);
    fc_chromosome_index_info ** tmp_chro_info_ptrs = malloc(global_context -> exontable_nchrs[f_idx] * sizeof(fc_chromosome_index_info *));
    
    global_context -> gene_name_array[f_idx] = malloc(sizeof(char *) * features);	// there should be much less identical names.
    global_context -> gene_name_table[f_idx] = HashTableCreate(5000);
    HashTableSetHashFunction(global_context -> gene_name_table[f_idx], HashTableStringHashFunction);
    HashTableSetKeyComparisonFunction(global_context -> gene_name_table[f_idx], fc_strcmp);
    
    long * new_start = malloc(sizeof(long) * features);
    long * new_end = malloc(sizeof(long) * features);
    int * new_entrez = malloc(sizeof(int) * features);
    unsigned char * new_strand = malloc(features);
    fc_feature_info_t ** new_info_ptr = malloc(sizeof(void *) * features);
    int new_features = 0;
    
    int * new_chro_start = malloc(sizeof(int) * global_context -> exontable_nchrs[f_idx]);
    
    // init start positions of each chromosome block.
    if(1)
    {
        KeyValuePair * cursor;
        int bucket;
        unsigned int sum_ptr = 0;
        for(bucket=0; bucket < global_context -> exontable_chro_table[f_idx] -> numOfBuckets; bucket++)
        {
            cursor = global_context -> exontable_chro_table[f_idx] -> bucketArray[bucket];
            while (1)
            {
                if (!cursor) break;
                fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
                cursor = cursor->next;
                //tmp_chro_inf -> reverse_table_end_index = calloc(sizeof(int), tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
                chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
                tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
            }
        }
        
        for(xk1 = 0; xk1 < global_context -> exontable_nchrs[f_idx]; xk1++)
        {
            unsigned int tmpv = chro_feature_ptr[xk1];
            chro_feature_ptr[xk1] = sum_ptr;
            tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
            sum_ptr += tmpv;
        }
        
    }
    int current_block_id = 0, sort_i = 0;
    
    (*sorted_entrezid) = new_entrez;
    (*sorted_start) = new_start;
    (*sorted_end) = new_end;
    (*sorted_strand) = new_strand;
    
    for(chro_pnt=0; chro_pnt < features; chro_pnt++) // for each feature, assign them to genes, chro_feature_ptr[chr] = #feature, ret_entrez[#feature] = #gene
    {
        char * this_chro_name = global_context -> unistr_buffer_space[f_idx] + loaded_features[chro_pnt].feature_name_pos + loaded_features[chro_pnt].chro_name_pos_delta;
        fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table[f_idx] , this_chro_name);
        assert(this_chro_info);
        unsigned int this_chro_number = this_chro_info -> chro_number;
        unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];
        //arrays storing info of each feature
        ret_entrez[this_chro_table_ptr] = find_or_insert_gene_name(global_context, (unsigned char *)(global_context -> unistr_buffer_space[f_idx] + loaded_features[chro_pnt].feature_name_pos), f_idx);
        ret_start[this_chro_table_ptr] = loaded_features[chro_pnt].start;
        ret_end[this_chro_table_ptr] = loaded_features[chro_pnt].end;
        ret_strand[this_chro_table_ptr] = loaded_features[chro_pnt].is_negative_strand;
        old_info_ptr[this_chro_table_ptr] = &loaded_features[chro_pnt];
        
        chro_feature_ptr[this_chro_number]++;
    }
    
    
    for(xk1 = 0; xk1 < global_context -> exontable_nchrs[f_idx]; xk1++)
    {
        fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
        int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
        short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);
        for(xk2=0; xk2<bins_in_chr; xk2++)
        {
            features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
        }
        
        memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) * bins_in_chr);
        
        tmp_chro_inf -> chro_block_table_start = current_block_id;
        unsigned int this_block_items = 0;
        long this_block_min_start = 0x7fffffff, this_block_max_end = 0;
        unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;
        int feature_i = 0, mark = 0;
        
        void * in_array[5];
        
        in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start;
        in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start;
        in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start;
        in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start;
        in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start;
        
        
        if (global_context -> is_strand_checked)
        {
            // Sort by strand
            merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, 2, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
            
            // For plus and minus strand, sort each by gene entrez
            while (feature_i < tmp_chro_inf -> chro_features - 1)
            {
                feature_i++;
                if (ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i] != ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i - 1] || feature_i == tmp_chro_inf -> chro_features - 1)
                {
                    void * gene_array[5];
                    gene_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start + mark;
                    gene_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start + mark;
                    gene_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start + mark;
                    gene_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start + mark;
                    gene_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start + mark;
                    if(feature_i == tmp_chro_inf -> chro_features - 1 && ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i] == ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i - 1])  // if last feature in the chromosome, and same strand as previous one, sort it together, otherwise leave it
                        merge_sort(gene_array, feature_i - mark + 1, 3, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
                    else
                        merge_sort(gene_array, feature_i - mark, 3, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
                    mark = feature_i;
                }
            }
            
        }
        else
        {
            merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, 3, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
        }
        
        
        // For each gene, sort by start position
        feature_i = 0;
        mark = 0;
        while (feature_i < tmp_chro_inf -> chro_features - 1)
        {
            feature_i++;
            if (ret_entrez[tmp_chro_inf -> chro_feature_table_start + feature_i] != ret_entrez[tmp_chro_inf -> chro_feature_table_start + feature_i - 1] || feature_i == tmp_chro_inf -> chro_features - 1 || (global_context-> is_strand_checked && ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i] != ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i - 1]))
            {
                void * gene_array[5];
                
                gene_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start + mark;
                gene_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start + mark;
                gene_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start + mark;
                gene_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start + mark;
                gene_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start + mark;
                if(feature_i == tmp_chro_inf -> chro_features - 1 && ((ret_entrez[tmp_chro_inf -> chro_feature_table_start + feature_i] == ret_entrez[tmp_chro_inf -> chro_feature_table_start + feature_i - 1])&& (!global_context-> is_strand_checked||(global_context-> is_strand_checked && ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i] == ret_strand[tmp_chro_inf -> chro_feature_table_start + feature_i - 1])))) // if last feature in the chromosome, and same gene and strand as previous one, sort it together, otherwise leave it
                    merge_sort(gene_array, feature_i - mark + 1, 0, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
                else
                    merge_sort(gene_array, feature_i - mark, 0, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
                mark = feature_i;
            }
        }
        
        new_chro_start[xk1] = new_features;
        
        void * new_array[5];
        new_array[0] = new_start + new_features;
        new_array[1] = new_end + new_features;
        new_array[2] = new_strand + new_features;
        new_array[3] = new_entrez + new_features;
        new_array[4] = new_info_ptr + new_features;
        
        // Merge module
        verse_feature_merge(global_context, in_array, new_array, &new_features, tmp_chro_inf);
        
        // Sort by start position after merge
        merge_sort(new_array, new_features - new_chro_start[xk1], 0, feature_sort_compare, feature_sort_exchange, feature_sort_merge);
        
        // block bin module
        for(sort_i = new_chro_start[xk1]; sort_i < new_features ; sort_i++)
        {
            // NOW THE FEATURES (ret_start, ret_end, ret_strand, ret_entrez, old_info_ptr) ARE ALL SORTED!
            // Check print here:
            // SUBREADprintf("CHRO=%d\tsort_i = %i\tentrez = %i\tgene_name = %s\tstrand = %i\tstart[sort_i]=%lu\tend[sort_i]=%lu\n", tmp_chro_inf->chro_number, sort_i, new_entrez[sort_i], global_context->gene_name_array[f_idx][new_entrez[sort_i]], new_strand[sort_i], new_start[sort_i], new_end[sort_i]);
            
            new_info_ptr[sort_i]->sorted_order = sort_i;
            
            int feature_bin_location = new_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH; //bin_local of the feature
            int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;
            
            if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location))
            {
                
                if(current_block_id >= current_block_buffer_size - 1)
                {
                    current_block_buffer_size *= 1.3;
                    ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
                    ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
                    ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
                }
                
                
                ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID //block_end_index[block#] = last feature# in the block
                ret_block_min_start[current_block_id] = this_block_min_start;
                ret_block_max_end[current_block_id] = this_block_max_end;
                register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
                //printf("B=%d; ST=%ld, END=%ld, ITM=%d\n", current_block_id, this_block_min_start, this_block_max_end, this_block_items);
                current_block_id++;
                this_block_max_end = 0;
                this_block_items = 0;
                this_block_min_start = 0x7fffffff;
            }
            
            this_block_max_end = max(this_block_max_end, new_end[sort_i]); //
            this_block_min_start = min(this_block_min_start, new_start[sort_i]);
            this_block_items ++;
        }
        
        if(this_block_items)
        {
            if(current_block_id >= current_block_buffer_size)
            {
                current_block_buffer_size *= 1.3;
                ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
                ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
                ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
            }
            
            ret_block_end_index[current_block_id] = new_features;
            ret_block_min_start[current_block_id] = this_block_min_start;
            ret_block_max_end[current_block_id] = this_block_max_end;
            register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
            current_block_id++;
        }
        
        tmp_chro_inf -> chro_block_table_end = current_block_id; 
        free(features_per_block_bins);
    }
    
    for(xk1 = 0; xk1 < global_context -> exontable_nchrs[f_idx]; xk1++)
    {
        fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
        tmp_chro_inf -> chro_feature_table_start = new_chro_start[xk1];
    }
    
    (*block_end_index) = ret_block_end_index;
    (*block_min_start_pos) = ret_block_min_start;
    (*block_max_end_pos) = ret_block_max_end;
    
    (*merged_features_ptr) = new_info_ptr;
    
    free(new_chro_start);
    free(ret_start);
    free(ret_end);
    free(ret_strand);
    free(ret_entrez);
    free(old_info_ptr);
    free(tmp_chro_info_ptrs);
    free(chro_feature_ptr);
    return new_features;
}

int strcmp_slash(char * s1, char * s2)
{
	char nch;
	while(0!=(nch = *(s1++))){
		if(nch == '/') break;
		if(nch != (*s2)) return 1;
		s2++;
	}
	return nch != *s2;
}

char * vs_tab_strtok(char * start)
{
    char * end = strchr(start,'\t');
    char * token = malloc(end - start + 1);
    sprintf(token, "%.*s", (int)(end - start), start);
    return token;
}

// This function returns the number of decision table items after processing the input section/read. It is called in section result combination and read result combination.
int build_decision_table(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, unsigned int counts, long * indices, int decision_table_items, int is_read_level, int f_idx)
{
    long * unique_gene_table = is_read_level? thread_context-> uniq_gene_table : thread_context -> sec_gene_table;
    long * unique_gene_exonid_table = is_read_level? thread_context -> uniq_gene_exonid_table : thread_context -> sec_gene_exonid_table;
    long * decision_table_ids = is_read_level? thread_context -> decision_table_ids : thread_context -> sec_table_ids;
    unsigned char * decision_table_votes = is_read_level? thread_context -> decision_table_votes : thread_context ->sec_table_votes;
    long * decision_table_exon_ids = is_read_level? thread_context -> decision_table_exon_ids : thread_context -> sec_table_exon_ids;
    int unique_genes = 0, xk1, xk2;
    
    // First loop, get the unique gene list of the current section/read
    for(xk1=0;xk1<counts;xk1++)
    {
        int gene_id = global_context -> exontable_geneid[f_idx][indices[xk1]];
        int is_unique = 1;
        for(xk2=0; xk2<unique_genes; xk2++)
        {
            if(gene_id == unique_gene_table[xk2])
            {
                is_unique = 0;
                break;
            }
        }
        if(is_unique){
            unique_gene_exonid_table[unique_genes] = indices[xk1];
            unique_gene_table[unique_genes++] = gene_id;
            if(unique_genes >= MAX_HIT_NUMBER) break;
        }
    }
    
    // Second loop, 1 vote for each of the unique genes. The decision table keeps track of the votes and gene_ids for all processed sections/reads.
    for(xk1=0;xk1<unique_genes; xk1++)
    {
        long gene_id = unique_gene_table[xk1];
        int is_fresh = 1;
        if(decision_table_items >= MAX_HIT_NUMBER) break;
        for(xk2=0; xk2<decision_table_items; xk2++)
        {
            if(gene_id == decision_table_ids[xk2])
            {
                decision_table_votes[xk2]++;
                is_fresh = 0;
                break;
            }
            
        }
        if(is_fresh)
        {
            decision_table_votes[decision_table_items] = 1;
            decision_table_exon_ids[decision_table_items] = unique_gene_exonid_table[xk1];
            decision_table_ids[decision_table_items++] = gene_id;
        }
    }
    return decision_table_items;
}

// This function will throw away features that don't reach the end of the last exon region. It will return the number of genes left. These genes are potential final features. If there is a next exon region, the features in that region will be checked to see if any of them is the exon of the potential features.
int region_clean_up(fc_thread_global_context_t * global_context, long region_end, long sec_end, int counts, long * indices, int f_idx)
{
    int k = 0, i;
    long end = 0;
    long tmp[counts];
    for (i = 0; i < counts; i++)
    {
        end = min(global_context -> exontable_stop[f_idx][indices[i]], sec_end);
        if (end == region_end)
            tmp[k++] = indices[i];
    }
    for (i = 0; i < k; i++)
    {
        indices[i] = tmp[i];
    }
    return k;
}

#define NO_CRITERIA 0
#define OVERLAP_THRESHOLD 1
#define NONOVERLAP_THRESHOLD 2
// This function calculate the total number of overlapped bases of each gene hit by the read. The hits_indices table and hits_total_length table will be updated. Exons will be merged so that there's going to be one record for one gene. Additional criteria can be applied to remove records that do not satisfy.
int length_sum_select(fc_thread_global_context_t * global_context, long * indices, short * hits_length, int nhits, short read_length, int criteria, int f_idx)
{
    // calculating the summed lengths of overlapping exons
    int x1;
    for(x1 = 0; x1 < nhits; x1++)
    {
        long exon_no = indices[x1];
        if(exon_no >= 0x7fffffff) continue;	//already removed
        int gene_no = global_context -> exontable_geneid[f_idx][exon_no];
        int x2;
        for(x2 = x1 + 1; x2 < nhits; x2++)
        {
            if(indices[x2] != 0x7fffffff && global_context -> exontable_geneid[f_idx][indices[x2]] == gene_no)
            {
                hits_length[x1] += hits_length[x2];
                indices[x2] = 0x7fffffff;
            }
        }
        if(criteria)
        {
            if(criteria == OVERLAP_THRESHOLD)
            {
                if(hits_length[x1] < global_context -> overlap_length_required)
                    indices[x1]=0x7fffffff;
            }
            if(criteria == NONOVERLAP_THRESHOLD)
            {
                if((read_length - hits_length[x1]) > global_context -> nonoverlap_length_allowed)
                    indices[x1]=0x7fffffff;
            }
        }
    }
    
    // remove the exons in indices & length tables when it is marked as removed (0x7fffffff)
    int new_hits = 0;
    for(x1 = 0; x1 < nhits; x1++)
    {
        if(indices[x1] >= 0x7fffffff) continue;
        if(new_hits != x1)
        {
            indices[new_hits] = indices[x1];
            hits_length[new_hits] = hits_length[x1];
        }
        new_hits++;
    }
    return new_hits;
}


// These exit functions are called by the read_assign(...) function. Any reads being processed by the read_assign(...) will have one of these fates.

void no_feature_exit(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * read_name, int f_idx)
{
    thread_context->read_counters[f_idx].unassigned_nofeatures ++;
    if(global_context -> detail_output_fp && f_idx == global_context -> feature_type_num - 1)
        fprintf(global_context -> detail_output_fp,"%s\tUnassigned_NoFeatures\t%s\t*\n", read_name, global_context->feature_type_list[f_idx]);
}

void ambiguous_exit(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * read_name, int f_idx)
{
    if(global_context -> detail_output_fp)
        fprintf(global_context -> detail_output_fp,"%s\tUnassigned_Ambiguit\t%s\t*\n", read_name, global_context->feature_type_list[f_idx]);
    thread_context->read_counters[f_idx].unassigned_ambiguous ++;
}

void assignment_exit(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * read_name, long hit_exon_id, int f_idx)
{
    thread_context->count_table[f_idx][hit_exon_id]++;
    if(global_context -> detail_output_fp)
    {
        int final_gene_number = global_context -> exontable_geneid[f_idx][hit_exon_id];
        unsigned char * final_feture_name = global_context -> gene_name_array[f_idx][final_gene_number];
        fprintf(global_context -> detail_output_fp,"%s\tAssigned\t%s\t%s\n", read_name, global_context->feature_type_list[f_idx], final_feture_name);
    }
    thread_context->read_counters[f_idx].assigned_reads ++;
}

// This function probe the read names of a buffer pair. It returns 0 if they are mates, 1 if they are not, -1 if no reads.
int probe_pair_name(fc_thread_thread_context_t * thread_context)
{
    char * line1 = thread_context -> line_buffer1;
    char * line2 = thread_context -> line_buffer2;
    
    char * read_name1 = vs_tab_strtok(line1);
    char * read_name2 = vs_tab_strtok(line2);
    int return_value;
    
    if(!read_name1 || !read_name2)
        return_value = -1;
    
    if(strcmp(read_name1,read_name2))
    {
        thread_context -> step_back = 1;
        return_value = 1;
    }
    else
        return_value = 0;
    free(read_name1);
    free(read_name2);
    return return_value;
}

// 1 sucessfully assigned to feature and pass the read, -1 not assigned and pass the read，0 not assigned but need further process
int read_assign(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, int f_idx, read_info_t read1, read_info_t read2)
{
    int nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id;
    long search_item_id, read_pos;
    long * hits_indices1 = thread_context -> hits_indices1, * hits_indices2 = thread_context -> hits_indices2;
    short * hits_total_length1 = thread_context -> hits_total_length1 ,  * hits_total_length2 = thread_context -> hits_total_length2;
    int is_second_read, allow_process = 1;
    unsigned int search_start = 0, search_end;
    short read_total_length1 = 0, read_total_length2 = 0, read_overlap_length1 = 0, read_overlap_length2 = 0;
    char * read_name = read1.read_name, * CIGAR_str, * read_chr;
    int treat_as_pair = global_context -> is_paired_end_data?1:0;
    if(thread_context -> step_back) treat_as_pair = 0;
    
    for(is_second_read = 0 ; is_second_read < treat_as_pair + 1; is_second_read++)
    {
        if(!is_second_read)
        {
            alignment_masks = read1.alignment_masks;
            CIGAR_str = read1.CIGAR_str;
            read_chr = read1.read_chr;
            read_pos = read1.read_pos;
            allow_process = read1.allow_process;
        }
        else
        {
            alignment_masks = read2.alignment_masks;
            CIGAR_str = read2.CIGAR_str;
            read_chr = read2.read_chr;
            read_pos = read2.read_pos;
            allow_process = read2.allow_process;
        }
        if(!CIGAR_str || (SAM_FLAG_UNMAPPED & alignment_masks) || !allow_process)
            continue;
        
		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0;
		int is_second_read_in_pair = alignment_masks & SAM_FLAG_SECOND_READ_IN_PAIR;
		int is_fragment_negative_strand = is_second_read_in_pair?(!is_this_negative_strand):is_this_negative_strand;

		fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table[f_idx], read_chr);

		if(this_chro_info == NULL)
		{
			if(this_chro_info == NULL && memcmp(read_chr, "chr", 3)==0)
			{
				this_chro_info = HashTableGet(global_context -> exontable_chro_table[f_idx], read_chr+3);
			}
			//printf("NL=%s, CI=%p\n", (read_chr), this_chro_info);
			if(this_chro_info == NULL && strlen(read_chr)<=2)
			{
				strcpy(thread_context -> chro_name_buff, "chr");
				strcpy(thread_context -> chro_name_buff+3, read_chr);
				this_chro_info = HashTableGet(global_context -> exontable_chro_table[f_idx], thread_context -> chro_name_buff);
			}
		}

		if(this_chro_info)
		{
			unsigned int nhits = 0;
            short read_total_length = 0;
            short read_overlap_length = 0;

			int cigar_section_id, cigar_sections;
			unsigned int Starting_Points[MAX_CIGAR_SECTIONS];
			unsigned short Section_Lengths[MAX_CIGAR_SECTIONS];
            unsigned short sec_overlap_length[MAX_CIGAR_SECTIONS] = {0};
			long * hits_indices = (is_second_read?hits_indices2:hits_indices1);
			short * hits_total_length = (is_second_read?hits_total_length2:hits_total_length1);
            int is_cigar_N = 0;
            
            // If the CIGAR string contains D or N, the read is divided into different sections.
			cigar_sections = RSubread_parse_CIGAR_string(CIGAR_str, Starting_Points, Section_Lengths, &is_cigar_N);
            long counts_indices[cigar_sections][MAX_HIT_NUMBER];
            unsigned int ncounts[cigar_sections]; // this records how many counts for each section
            memset(ncounts, 0, cigar_sections*sizeof(int));

            
            if(global_context -> reduce_5_3_ends_to_one)
            {
                if((REDUCE_TO_5_PRIME_END == global_context -> reduce_5_3_ends_to_one) + is_this_negative_strand == 1) // reduce to 5' end (small coordinate if positive strand / large coordinate if negative strand)
                {
                    Section_Lengths[0]=1;
                }
                else
                {
                    Starting_Points[0] = Starting_Points[cigar_sections-1] + Section_Lengths[cigar_sections-1] - 1;
                    Section_Lengths[0]=1;
                }
                
                cigar_sections = 1;
            }
            
            // Extending the reads to the 3' and 5' ends. (from the read point of view)
            if(global_context -> five_end_extension)
            {
                if(is_this_negative_strand){
                    Section_Lengths [cigar_sections - 1] += global_context -> five_end_extension;
                }else{
                    //SUBREADprintf("5-end extension: %d [%d]\n", Starting_Points[0], Section_Lengths[0]);
                    if( read_pos > global_context -> five_end_extension)
                    {
                        Section_Lengths [0] += global_context -> five_end_extension;
                        read_pos  -= global_context -> five_end_extension;
                    }
                    else
                    {
                        Section_Lengths [0] += read_pos-1;
                        read_pos = 1;
                    }
                }
            }
            
            if(global_context -> three_end_extension)
            {
                
                if(is_this_negative_strand)
                {
                    if( read_pos > global_context -> three_end_extension)
                    {
                        Section_Lengths [0] += global_context -> three_end_extension;
                        read_pos -= global_context -> three_end_extension;
                    }
                    else
                    {
                        Section_Lengths [0] += read_pos - 1;
                        read_pos = 1;
                    }
                }
                else	Section_Lengths [cigar_sections - 1] += global_context -> three_end_extension;
            }
            
            //#warning "=================== COMMENT THESE 2 LINES ================================"
            //for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
            //	SUBREADprintf("ACCC: %llu , sec[%d] %u ~ %d ; secs=%d\n", read_pos, cigar_section_id, Starting_Points[cigar_section_id], Section_Lengths[cigar_section_id], cigar_sections);
            if( global_context -> debug_command [0] == 'D')
                SUBREADprintf("\n\nRead name = %s ; Second read = %d\n", read_name, is_second_read);
            
            for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
            {
                read_total_length += Section_Lengths[cigar_section_id];
                
                long section_begin_pos = read_pos + Starting_Points[cigar_section_id];
                long section_end_pos = Section_Lengths[cigar_section_id] + section_begin_pos - 1;
                
                long region[2] = {0};   // This keeps the start and end position of the last exon region (mode 3 only)
                long sec_mark[2] = {0}; // This marks the start and end position of the exon region that needs to be covered (mode 3 only)
                long sec_mark_id = 0; // This stores the sec_marked exon id
                int exon_region_counts = 0;    // This tracks the number of exon regions in this sections
                
                if(global_context -> debug_command [0] == 'D')
                    SUBREADprintf("Section [%d]: %d ~ %d\n", cigar_section_id , Starting_Points[cigar_section_id], Starting_Points[cigar_section_id]+Section_Lengths[cigar_section_id]);
                
                int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH; // start bin of the section
                int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH; // end bin of the section
                
                start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
                end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);
                
                while(start_reverse_table_index<=end_reverse_table_index)
                {
                    search_start = this_chro_info -> reverse_table_start_index [start_reverse_table_index]; //search_start:  the start block in the bin
                    if(search_start<0xffffff00)break;
                    start_reverse_table_index++;
                }
                if(search_start>0xffffff00) continue;
                
                search_end = this_chro_info -> chro_block_table_end;
                
                for(search_block_id=search_start;search_block_id<search_end;search_block_id++)
                {
                    if (global_context -> exontable_block_min_start[f_idx][search_block_id] > section_end_pos) break;
                    if (global_context -> exontable_block_max_end[f_idx][search_block_id] < section_begin_pos) continue;
                    // the section(read) must be within the block.
                    
                    int search_item_start = 0, search_item_end = global_context -> exontable_block_end_index[f_idx][search_block_id];
                    if(search_block_id>0)search_item_start = global_context -> exontable_block_end_index[f_idx][search_block_id-1];
                    
                    for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++)
                    {
                        if (global_context -> exontable_stop[f_idx][search_item_id] >= section_begin_pos)
                        {
                            if (global_context -> exontable_start[f_idx][search_item_id] > section_end_pos) break;
                            // if the start position of feature i is larger than end pos of section, then no overlap, otherwise there is an overlap
                            // the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
                            
                            int is_strand_ok = 1;
                            if(global_context->is_strand_checked){
                                if(global_context->is_strand_checked == 1)
                                    is_strand_ok = (is_fragment_negative_strand == global_context -> exontable_strand[f_idx][search_item_id]);
                                else
                                    is_strand_ok = (is_fragment_negative_strand != global_context -> exontable_strand[f_idx][search_item_id]);
                                //SUBREADprintf("%d = %d == %d\n", is_strand_ok, is_fragment_negative_strand, global_context -> exontable_strand[search_item_id]);
                            }
                            
                            if(is_strand_ok){
                                if(nhits<=MAX_HIT_NUMBER - 1)
                                {
                                    if (global_context -> run_mode < 2 || global_context -> run_mode == 5) // featurecounts, union, and cover_length
                                    {
                                        hits_indices[nhits] = search_item_id;
                                        
                                        if(global_context -> overlap_length_required != 1 || global_context -> nonoverlap_length_allowed != 0x7fff || global_context -> run_mode == 5)
                                        {
                                            int section_overlapped = min(global_context -> exontable_stop[f_idx][search_item_id], section_end_pos)
                                            - max(global_context -> exontable_start[f_idx][search_item_id], section_begin_pos) + 1;
                                            hits_total_length[nhits] = (short)section_overlapped;
                                        }
                                        nhits++;
                                    }
                                    
                                    if (global_context -> run_mode == 2 || global_context -> run_mode == 4)     // If is strict mode, This throw away any read with incomplete feature cover
                                    {
                                        long r[2] = {0};
                                        r[0] = max(global_context -> exontable_start[f_idx][search_item_id] , section_begin_pos);
                                        r[1] = min(global_context -> exontable_stop[f_idx][search_item_id] , section_end_pos);
                                        if (r[0] > section_begin_pos)   // The feature do not start at the section start
                                        {
                                            if (ncounts[cigar_section_id] == 0)
                                            {
                                                no_feature_exit(global_context, thread_context, read_name, f_idx);
                                                return 0;
                                            }
                                            else
                                                goto section_loop_end;  // at least one cover for the read, skip the section, go on to check next one;
                                        }
                                        else if (r[1] < section_end_pos)    // r[0] == section_begin_pos, The feature do not reach the section end
                                        {
                                            continue;
                                        }
                                        else        // a covering feature!
                                        {
                                            counts_indices[cigar_section_id][ncounts[cigar_section_id]++] = search_item_id;
                                        }
                                    }
                                    
                                    else if (global_context -> run_mode == 3)   // Intersection-nonempty mode
                                    {
                                        long r[2] = {0};
                                        int i = 0;
                                        r[0] = max(global_context -> exontable_start[f_idx][search_item_id] , section_begin_pos);
                                        r[1] = min(global_context -> exontable_stop[f_idx][search_item_id] , section_end_pos);
                                        
                                        if (region[1] == 0) // first hit of the section
                                        {
                                            region[0] = r[0];
                                            region[1] = r[1];
                                            sec_overlap_length[cigar_section_id] = r[1] - r[0] + 1;
                                            counts_indices[cigar_section_id][ncounts[cigar_section_id]++] = search_item_id;
                                        }
                                        else if (!sec_mark[1])       // no section needs to be covered
                                        {
                                            if (r[0] == region[0])
                                            {
                                                if (r[1] == region[1])
                                                {
                                                    if (!exon_region_counts)   // in the first exon region
                                                        counts_indices[cigar_section_id][ncounts[cigar_section_id]++] = search_item_id;
                                                    else        // in exon region other than the first
                                                    {
                                                        for(i = 0; i < ncounts[cigar_section_id]; i++)
                                                        {
                                                            if (global_context -> exontable_geneid[f_idx][search_item_id] == global_context -> exontable_geneid[f_idx][counts_indices[cigar_section_id][i]])
                                                                counts_indices[cigar_section_id][i] = search_item_id;
                                                        }
                                                    }
                                                }
                                                else if (r[1] > region[1])
                                                {
                                                    sec_overlap_length[cigar_section_id] += r[1] - region[1];
                                                    if (!exon_region_counts)   // in the first exon region
                                                    {
                                                        counts_indices[cigar_section_id][0] = search_item_id;
                                                        ncounts[cigar_section_id] = 1;
                                                        region[1] = r[1];
                                                    }
                                                    else        // in exon region other than the first
                                                    {
                                                        int exon_flag = 0;
                                                        for(i = 0; i < ncounts[cigar_section_id]; i++)
                                                        {
                                                            if (global_context -> exontable_geneid[f_idx][search_item_id] == global_context -> exontable_geneid[f_idx][counts_indices[cigar_section_id][i]])
                                                            {
                                                                exon_flag = 1;
                                                                counts_indices[cigar_section_id][i] = search_item_id;
                                                                region[1] = r[1];
                                                            }
                                                        }
                                                        if (!exon_flag)
                                                        {
                                                            sec_mark[0] = region[0];
                                                            sec_mark[1] = r[1];
                                                            sec_mark_id = search_item_id;
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            else if (r[1] > region[1])      // && (r[0] > region[0])
                                            {
                                                if (r[0] <= region[1])
                                                {
                                                    if(!global_context->is_nonempty_modified)
                                                    {
                                                        no_feature_exit(global_context, thread_context, read_name, f_idx);
                                                        return 0;
                                                    }
                                                    else
                                                    {
                                                        ambiguous_exit(global_context, thread_context, read_name, f_idx);
                                                        return -1;
                                                    }
                                                }
                                                else // r[0] > region[1]
                                                {
                                                    // clean up previous exon region
                                                    ncounts[cigar_section_id] = region_clean_up(global_context, region[1], section_end_pos, ncounts[cigar_section_id], counts_indices[cigar_section_id], f_idx);
                                                    
                                                    // set up new exon region
                                                    sec_overlap_length[cigar_section_id] += r[1] - r[0] + 1;
                                                    exon_region_counts += 1;
                                                    int exon_flag = 0;
                                                    for (i = 0; i < ncounts[cigar_section_id]; i++)
                                                    {
                                                        if (global_context -> exontable_geneid[f_idx][search_item_id] == global_context -> exontable_geneid[f_idx][counts_indices[cigar_section_id][i]])
                                                        {
                                                            exon_flag = 1;
                                                            counts_indices[cigar_section_id][i] = search_item_id;
                                                            region[1] = r[1];
                                                            region[0] = r[0];
                                                        }
                                                    }
                                                    if(!exon_flag)
                                                    {
                                                        sec_mark[0] = r[0];
                                                        sec_mark[1] = r[1];
                                                        sec_mark_id = search_item_id;
                                                    }
                                                }
                                            }
                                            else
                                                continue;
                                        }
                                        else                    // if sec_mark
                                        {
                                            if (r[0] > sec_mark[0])     // the "foreign" exon will not be covered
                                            {
                                                if(!global_context->is_nonempty_modified)
                                                {
                                                    no_feature_exit(global_context, thread_context, read_name, f_idx);
                                                    return 0;
                                                }
                                                else
                                                {
                                                    ambiguous_exit(global_context, thread_context, read_name, f_idx);
                                                    return -1;
                                                }
                                            }
                                            else if (r[1] >= sec_mark[1])   // && r[0] == sec_mark[0]
                                            {
                                                int exon_flag = 0;
                                                sec_overlap_length[cigar_section_id] += r[1] - sec_mark[1];
                                                for(i = 0; i < ncounts[cigar_section_id]; i++)
                                                {
                                                    if (global_context -> exontable_geneid[f_idx][search_item_id] == global_context -> exontable_geneid[f_idx][counts_indices[cigar_section_id][i]])
                                                    {
                                                        exon_flag = 1;
                                                        counts_indices[cigar_section_id][i] = search_item_id;
                                                        region[0] = r[0];
                                                        region[1] = r[1];
                                                        sec_mark[0] = 0;
                                                        sec_mark[1] = 0;
                                                    }
                                                }
                                                if (!exon_flag)
                                                {
                                                    sec_mark[1] = r[1];
                                                    sec_mark_id = search_item_end;
                                                }
                                            }
                                        }
                                    }
                                    
                                }
                                else break;
                            }
                        }
                    }
                    
                }
                if (global_context -> run_mode == 3)
                {
                    if (sec_mark[1])
                    {
                        if(!global_context->is_nonempty_modified)
                        {
                            no_feature_exit(global_context, thread_context, read_name, f_idx);
                            return 0;
                        }
                        else
                        {
                            ambiguous_exit(global_context, thread_context, read_name, f_idx);
                            return -1;
                        }
                    }
                    else if (exon_region_counts)    // The last exon region needs to be cleaned up
                        ncounts[cigar_section_id] = region_clean_up(global_context, region[1], section_end_pos, ncounts[cigar_section_id], counts_indices[cigar_section_id], f_idx);
                }
            section_loop_end:;
            } // section loop end;
            
            if (global_context -> run_mode == 3)    // Calculate the overlap length
            {
                for(cigar_section_id = 0; cigar_section_id < cigar_sections; cigar_section_id++)
                {
                    read_overlap_length += sec_overlap_length[cigar_section_id];
                }
            }
            
            // merge counts_indices (for each section) to hits_indices (for each read)
            if (global_context -> run_mode == 2 || global_context -> run_mode == 3 || global_context -> run_mode == 4)   
            {
                nhits = 0;
                int sec_table_items = 0, xk1;
                unsigned char * sec_table_votes = thread_context ->sec_table_votes;
                long * sec_table_exon_ids = thread_context -> sec_table_exon_ids;
                int nsections = cigar_sections;
                int sum_counts = 0, sec, i = 0;
                
                if (cigar_sections == 1)        // only one section
                {
                    nhits = ncounts[0];
                    for (i = 0; i < nhits; i++)
                    {
                        hits_indices[i] = counts_indices[0][i];
                    }
                }
                
                else        // combining section results to read result
                {
                    for (sec = 0; sec < cigar_sections; sec++)
                    {
                        if (ncounts[sec]<1)
                        {
                            if (global_context -> run_mode != 3)
                            {
                                no_feature_exit(global_context, thread_context, read_name, f_idx);
                                return 0;
                            }
                            else        // run_mode == 3
                            {
                                nsections--;
                                continue;
                            }
                        }
                        sec_table_items = build_decision_table(global_context, thread_context, ncounts[sec], counts_indices[sec], sec_table_items, SECTION_LEVEL, f_idx);
                        sum_counts += ncounts[sec];
                    }
                    
                    if (global_context -> run_mode == 4)
                    {
                        if (sec_table_items < 1)
                        {
                            no_feature_exit(global_context, thread_context, read_name, f_idx);
                            return 0;
                        }
                        else if(sec_table_items > 1)
                        {
                            ambiguous_exit(global_context, thread_context, read_name, f_idx);
                            return -1;
                        }
                        else                        // one gene one read
                        {
                            nhits = 1;
                            hits_indices[0] = sec_table_exon_ids[0];
                        }
                    }
                    else if (sec_table_items > 0)       // run mode 2, 3 do the intersection
                    {
                        int max_votes = 0;
                        int top_voters = 0;
                        long top_voter_id = 0;
                        
                        for(xk1 = 0; xk1 < sec_table_items; xk1++)
                        {
                            if(sec_table_votes[xk1] > max_votes)
                            {
                                max_votes = sec_table_votes[xk1];
                                top_voters = 1;
                                top_voter_id = sec_table_exon_ids[xk1];
                            }
                            else
                                if(sec_table_votes[xk1] == max_votes) top_voters++;
                        }
                        
                        if (max_votes < nsections) // No gene overlap all voted sections, empty set
                        {
                            if(!global_context->is_nonempty_modified)
                            {
                                no_feature_exit(global_context, thread_context, read_name, f_idx);
                                return 0;
                            }
                            else
                            {
                                ambiguous_exit(global_context, thread_context, read_name, f_idx);
                                return -1;
                            }
                        }
                        else
                        {
                            nhits = top_voters; // Set nhits to the number of hit genes for the read
                            for (xk1 = 0; xk1 < sec_table_items; xk1++)
                            {
                                if (sec_table_votes[xk1] == max_votes)
                                {
                                    hits_indices[i++] = sec_table_exon_ids[xk1];
                                }
                            }
                        }
                    }
                }
                
                if (nhits == 0 && (global_context->run_mode == 2 || global_context->run_mode == 4 || (global_context->run_mode == 3 && sum_counts > 1)))
                {
                    if(!global_context->is_nonempty_modified)
                    {
                        no_feature_exit(global_context, thread_context, read_name, f_idx);
                        return 0;
                    }
                    else
                    {
                        ambiguous_exit(global_context, thread_context, read_name, f_idx);
                        return -1;
                    }
                }
                
            }
            
            if(is_second_read)
            {
                nhits2 = nhits;
                read_total_length2 = read_total_length;
                read_overlap_length2 = read_overlap_length;
            }
            else
            {
                nhits1 = nhits;
                read_total_length1 = read_total_length;
                read_overlap_length1 = read_overlap_length;
            }
        }
	}       // read level
    
    // [Hidden function]: This is used to select the reads with total M section length larger than a threshold (say 80 bp). Currently it cannot be specified by user.
    if(global_context -> read_M_len_required != 1)
    {
        if(read_total_length1 < global_context -> read_M_len_required || (treat_as_pair && read_total_length2 < global_context -> read_M_len_required))
        {
            no_feature_exit(global_context, thread_context, read_name, f_idx);
            return 0;
        }
    }
    
	if(global_context -> overlap_length_required != 1)
	{
        if(global_context -> run_mode < 2 || global_context -> run_mode == 5)      // featurecounts, union, and cover_length
        {
            // two ends in a read pair is considered independently; the overlapping bases are not added up.
            int ends;
            for(ends = 0; ends < treat_as_pair + 1; ends++)
            {
                long * hits_indices = ends?hits_indices2:hits_indices1;
                short * hits_total_length = ends?hits_total_length2:hits_total_length1;
                int nhits = ends?nhits2:nhits1;
                short read_total_length = ends?read_total_length2:read_total_length1;
                int new_hits = 0;
                
                new_hits = length_sum_select(global_context, hits_indices, hits_total_length, nhits, read_total_length, OVERLAP_THRESHOLD, f_idx);
                if(ends) nhits2 = new_hits;
                else nhits1 = new_hits;
            }
        }
        
        // In intersection mode, all assigned features must overlap the same bases, so it wouldn't be necessary to calculate the overlap length for each hit
        if(global_context -> run_mode == 2 || global_context -> run_mode == 4)      // strict mode require the assigned feature to overlap every single base of the read
        {
            if(read_total_length1 < global_context -> overlap_length_required || (treat_as_pair && read_total_length2 < global_context -> overlap_length_required))
            {
                no_feature_exit(global_context, thread_context, read_name, f_idx);
                return 0;
            }
        }
        
        if(global_context -> run_mode == 3)     // intersection_nonempty mode, length stored in read_overlap_length;
        {
            if(read_overlap_length1 < global_context -> overlap_length_required || (treat_as_pair && read_overlap_length2 < global_context -> overlap_length_required))
            {
                no_feature_exit(global_context, thread_context, read_name, f_idx);
                return 0;
            }
        }
	}
    
    if(global_context -> nonoverlap_length_allowed != 0x7fff)
    {
        if(global_context -> run_mode < 2 || global_context -> run_mode == 5)
        {
            // two ends in a read pair is considered independently; the overlapping bases are not added up.
            int ends;
            for(ends = 0; ends < treat_as_pair + 1; ends++)
            {
                long * hits_indices = ends?hits_indices2:hits_indices1;
                short * hits_total_length = ends?hits_total_length2:hits_total_length1;
                int nhits = ends?nhits2:nhits1;
                short read_total_length = ends?read_total_length2:read_total_length1;
                int new_hits = 0;

                new_hits = length_sum_select(global_context, hits_indices, hits_total_length, nhits, read_total_length, NONOVERLAP_THRESHOLD, f_idx);
                if(ends) nhits2 = new_hits;
                else nhits1 = new_hits;
            }
        }
        
        if(global_context -> run_mode == 3)
        {
            short non_overlap_length1 = read_total_length1 - read_overlap_length1;
            short non_overlap_length2 = read_total_length2 - read_overlap_length2;
            if ((non_overlap_length1 > global_context -> nonoverlap_length_allowed || (treat_as_pair && non_overlap_length2 > global_context -> nonoverlap_length_allowed)) && (read_total_length1 != 0 && (treat_as_pair && read_total_length2 != 0)))
            {
                no_feature_exit(global_context, thread_context, read_name, f_idx);
                return 0;
            }
        }
    }
    
    // assign reads
    
    if(nhits2 + nhits1 == 0)
    {
        no_feature_exit(global_context, thread_context, read_name, f_idx);
        return 0;
    }

	else if(nhits2 + nhits1 == 1)
	{
		long hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];
        assignment_exit(global_context, thread_context, read_name, hit_exon_id, f_idx);
        return 1;
	}
    
	else if(nhits2 == 1 && nhits1 == 1 && (global_context -> exontable_geneid[f_idx][hits_indices2[0]] == global_context ->exontable_geneid[f_idx][hits_indices1[0]]))
	{
		long hit_exon_id = hits_indices1[0];
        assignment_exit(global_context, thread_context, read_name, hit_exon_id, f_idx);
        return 1;
	}
    
	else
	{
		if(nhits2+nhits1>=MAX_HIT_NUMBER)       // hits exceeds MAX_HIT_NUMBER, warn
		{
			SUBREADprintf("A %s overlapped with %d features.\n", global_context -> is_paired_end_data?"read pair":"read", nhits2+nhits1);
			SUBREADprintf("Please increase MAX_HIT_NUMBER in the program.\n");
		}

		unsigned char * decision_table_votes = thread_context ->decision_table_votes;
		long * decision_table_exon_ids = thread_context -> decision_table_exon_ids;
		int decision_table_items = 0;

        for(is_second_read = 0; is_second_read < treat_as_pair + 1; is_second_read++)
        {
            long * hits_indices = is_second_read?hits_indices2:hits_indices1;
            int nhits = is_second_read?nhits2:nhits1;
            if (nhits<1) continue;
            decision_table_items = build_decision_table(global_context, thread_context, nhits, hits_indices, decision_table_items, READ_LEVEL, f_idx);
        }

         // Union mode, if nhits1 and nhits 2 together hit more than 2 features, 1:1, 2:0, 0:2, 1:2 ..., then after previous step, if decision table items is 1, then only one gene is hit (multiple exons), if more than 1, read pair is set ambiguous.
        if (global_context -> run_mode == 1 || global_context -> run_mode == 4)
        {
            if(decision_table_items == 1) {
                long hit_exon_id = decision_table_exon_ids[0];
                assignment_exit(global_context, thread_context, read_name, hit_exon_id, f_idx);
                return 1;
            }
            else if (decision_table_items > 1) {
                ambiguous_exit(global_context, thread_context, read_name, f_idx);
                return -1;
            }
        }
        
        // Cover_length mode, if hit 1 gene, assign; if more than 1, calculate length coverage
        else if (global_context -> run_mode == 5)
        {
            if(decision_table_items == 1) {
                long hit_exon_id = decision_table_exon_ids[0];
                assignment_exit(global_context, thread_context, read_name, hit_exon_id, f_idx);
                return 1;
            }
            else if (decision_table_items > 1)
            {
                // append hits_indices2 to hits_indices1, hits_total_length2 to hits_total_length1
                int i = 0;
                for(i = 0; i < nhits2; i++)
                {
                    hits_indices1[nhits1] = hits_indices2[i];
                    hits_total_length1[nhits1++] = hits_total_length2[i];
                }
                
                // calculate length
                nhits1 = length_sum_select(global_context, hits_indices1, hits_total_length1, nhits1, read_total_length1 + read_total_length2, NO_CRITERIA, f_idx);
                
                // get the longest and second longest
                short max1_len = hits_total_length1[0], max2_len = 0;
                long max1_idx = hits_indices1[0], max2_idx = 0;
                for(i = 1; i < nhits1; i++)
                {
                    if(hits_total_length1[i] <= max1_len && hits_total_length1[i] > max2_len)
                    {
                        max2_len = hits_total_length1[i];
                        max2_idx = hits_indices1[i];
                    }
                    if(hits_total_length1[i] > max1_len)
                    {
                        max2_len = max1_len;
                        max2_idx = max1_idx;
                        max1_len = hits_total_length1[i];
                        max1_idx = hits_indices1[i];
                    }
                }
                
                // assign
                if (max1_len - max2_len > global_context -> min_dif_ambiguous)  // by default > 0
                {
                    long hit_exon_id = max1_idx;
                    assignment_exit(global_context, thread_context, read_name, hit_exon_id, f_idx);
                    return 1;
                }
                else
                {
                    ambiguous_exit(global_context, thread_context, read_name, f_idx);
                    return -1;
                }
                
            }
        }
        
        // featureCounts/intersection mode, find out the result with max vote, otherwise set ambiguous
		else if((global_context -> run_mode == 0 || global_context -> run_mode == 2 || global_context -> run_mode == 3) && decision_table_items > 0)
		{
			int max_votes = 0, xk1;
			int top_voters = 0;
			long top_voter_id = 0;

			for(xk1 = 0; xk1 < decision_table_items; xk1++)
			{
				if(decision_table_votes[xk1] > max_votes)
				{
					max_votes = decision_table_votes[xk1];
					top_voters = 1;
					top_voter_id = decision_table_exon_ids[xk1];
				}
				else if(decision_table_votes[xk1] == max_votes) top_voters++;
			}
            
            // One gene was hit by 2 reads
			if(top_voters == 1)
			{
                assignment_exit(global_context, thread_context, read_name, top_voter_id, f_idx);
                return 1;
			}
            
            else if(top_voters > 1)
            {
                // If featureCounts mode or single end read, this indicate that the read hit multiple genes.
                if(global_context -> run_mode == 0 || !treat_as_pair)
                {
                    ambiguous_exit(global_context, thread_context, read_name, f_idx);
                    return -1;
                }
                if(global_context -> run_mode == 2)
                {
                    // No gene was hit by both reads
                    if(max_votes == 1)
                    {
                        no_feature_exit(global_context, thread_context, read_name, f_idx);
                        return 0;
                    }
                    // More than one gene were hit by both reads
                    else
                    {
                        ambiguous_exit(global_context, thread_context, read_name, f_idx);
                        return -1;
                    }
                }
                if (global_context -> run_mode == 3)
                {
                    // More than one gene were hit by both reads, or (n:0), (0:n)
                    if (max_votes > 1 || (nhits1 == 0 || nhits2 == 0))
                    {
                        ambiguous_exit(global_context, thread_context, read_name, f_idx);
                        return -1;
                    }
                    else
                    {
                        if(!global_context->is_nonempty_modified)
                        {
                            no_feature_exit(global_context, thread_context, read_name, f_idx);
                            return 0;
                        }
                        else
                        {
                            ambiguous_exit(global_context, thread_context, read_name, f_idx);
                            return -1;
                        }
                    }
                }
            }
		}
        
		else{
            no_feature_exit(global_context, thread_context, read_name, f_idx);
            return 0;
		}
	}
    return 0;
}

// This function is the main logic step for read assignment scheme (hierachical or independent). It also preprocess the read buffer and prevent unsplit, unmapped, multimapper, low_qual, duplicate reads from entering the read_assign(...) step, if specifed (or by default).
void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{
    thread_context->all_reads++;
    int f_idx;
    read_info_t read1, read2;
    char * tmp_tok_ptr= NULL;
    long fragment_length = 0;
    int is_second_read;
    int first_read_quality_score = 0;
    int skipped_read = 0;
    
    int treat_as_pair = global_context -> is_paired_end_data?1:0;
    if(thread_context -> step_back) treat_as_pair = 0;
    
    read1.read_name = NULL; read2.read_name = NULL; read1.alignment_masks = 0; read2.alignment_masks = 0; read1.CIGAR_str = NULL; read2.CIGAR_str = NULL; read1.read_pos = 0; read2.read_pos = 0; read1.read_chr = NULL; read2.read_chr = NULL; read1.allow_process = 1; read2.allow_process = 1;
    
    for(is_second_read = 0 ; is_second_read < treat_as_pair + 1; is_second_read++)
    {
        
        char * line = is_second_read? thread_context -> line_buffer2:thread_context -> line_buffer1;
        
        //printf("LINE_BUF=%s\n",line);
        
        read2.read_name = strtok_r(line,"\t", &tmp_tok_ptr);
        if(!read2.read_name)return;
        if(!is_second_read)
            read1.read_name = read2.read_name;
        
        char * mask_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
        if((!mask_str) || !isdigit(mask_str[0])) return;
        read2.alignment_masks = atoi(mask_str);
        
        if(!is_second_read)
            read1.alignment_masks = read2.alignment_masks;
        
        if(is_second_read == 0)
        {
            // skip the read if unmapped
            if(((!treat_as_pair) &&  (read1.alignment_masks & SAM_FLAG_UNMAPPED) ) ||
               ((read1.alignment_masks & SAM_FLAG_UNMAPPED)   &&  (read1.alignment_masks & SAM_FLAG_MATE_UNMATCHED) && treat_as_pair) ||
               (((read1.alignment_masks & SAM_FLAG_UNMAPPED) || (read1.alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && treat_as_pair && global_context -> is_both_end_required))
            {
                thread_context->read_counters[0].unassigned_unmapped ++;
                if(global_context -> detail_output_fp)
                    fprintf(global_context -> detail_output_fp,"%s\tUnassigned_Unmapped\t*\t*\t*\n", read1.read_name);
                return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
            }
        }
        
        read2.read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);
        if(!read2.read_chr) return;
        char * read_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
        if(!read_pos_str) return;
        
        read2.read_pos = atoi(read_pos_str);
        if(read2.read_pos < 1 && read_pos_str[0]!='0') return;
        if(!is_second_read)
        {
            read1.read_chr = read2.read_chr;
            read1.read_pos = read2.read_pos;
        }
        char * mapping_qual_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
        
        read2.CIGAR_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
        if(!is_second_read)
            read1.CIGAR_str = read2.CIGAR_str;
        if(!read2.CIGAR_str)
            continue;
        
        if(global_context->is_split_alignments_only)
        {
            unsigned int Starting_Points[MAX_CIGAR_SECTIONS];
            unsigned short Section_Lengths[MAX_CIGAR_SECTIONS];
            int is_cigar_N = 0; // This indicates whether this CIGAR is N divided or not (there are also D divided sections).
            RSubread_parse_CIGAR_string(read2.CIGAR_str, Starting_Points, Section_Lengths, &is_cigar_N);
            if (!is_cigar_N)    // This is not a split read
                skipped_read ++;
            // If both ends in read pair is skipped (not split), or the read missing a mate is not a split read, or the read containing an unmapped mate is not a split read, then throw it to nonjunction.
            if((is_second_read && skipped_read == 2) || (skipped_read && (thread_context -> step_back || (read2.alignment_masks & 0x8))))
            {
                thread_context->read_counters[0].unassigned_nonjunction ++;
                if(global_context -> detail_output_fp)
                    fprintf(global_context -> detail_output_fp,"%s\tUnassigned_Nonjunction\t*\t*\t*\n", read1.read_name);
                return;
            }
        }
        
        if(global_context -> min_mapping_quality_score>0)
        {
            int mapping_qual =atoi(mapping_qual_str);
            if(( mapping_qual < global_context -> min_mapping_quality_score  && ! treat_as_pair)||( is_second_read  && max( first_read_quality_score, mapping_qual ) < global_context -> min_mapping_quality_score))
            {
                thread_context->read_counters[0].unassigned_mappingquality ++;
                if(global_context -> detail_output_fp)
                {
                    fprintf(global_context -> detail_output_fp,"%s\tUnassigned_MappingQuality\t*\t*\tMapping_Quality=%d,%d\n", read1.read_name, first_read_quality_score, mapping_qual);
                }
                return;
            }
            if(is_second_read==0 && treat_as_pair)
            {
                first_read_quality_score = mapping_qual;
            }
        }
        
        long mate_pos = 0;
        char * mate_chr = NULL;
        
        if(is_second_read)
        {
            mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);// mate_chr
            if(mate_chr[0]=='=') mate_chr = read2.read_chr;
            char * mate_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);	// mate_pos
            mate_pos = atol(mate_pos_str);
        }
        
        if(is_second_read == 0 && treat_as_pair &&
           (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed))
        {
            int is_half_mapped = (read1.alignment_masks & SAM_FLAG_UNMAPPED) || (read1.alignment_masks & SAM_FLAG_MATE_UNMATCHED);
            
            if(!is_half_mapped)
            {
                char * mate_chrx = strtok_r(NULL,"\t", &tmp_tok_ptr);
                if(!mate_chrx) return;
                strtok_r(NULL,"\t", &tmp_tok_ptr);
                if(!tmp_tok_ptr) return;
                char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
                if(!tmp_tok_ptr) return;
                
                fragment_length = abs(atoi(frag_len_str));
                
                int is_first_read_negative_strand = (read1.alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0;
                int is_second_read_negative_strand = (read1.alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0;
                
                if(mate_chrx[0]=='=' && is_first_read_negative_strand!=is_second_read_negative_strand)
                {
                    if(global_context -> is_PE_distance_checked && ((fragment_length > global_context -> max_paired_end_distance) || (fragment_length < global_context -> min_paired_end_distance)))
                    {
                        thread_context->read_counters[0].unassigned_fragmentlength ++;
                        
                        if(global_context -> detail_output_fp)
                            fprintf(global_context -> detail_output_fp,"%s\tUnassigned_TLEN\t*\t*\tLength=%ld\n", read1.read_name, fragment_length);
                        return;
                    }
                }
                else
                {
                    if(global_context -> is_chimertc_disallowed)
                    {
                        thread_context->read_counters[0].unassigned_chimericreads ++;
                        
                        if(global_context -> detail_output_fp)
                            fprintf(global_context -> detail_output_fp,"%s\tUnassigned_Chimera\t*\t*\t*\n", read1.read_name);
                        return;
                    }
                }
            }
        }
        
        if(!tmp_tok_ptr) return;
        
        
        // This filter has to be put here because the 0x400 FLAG is not about mapping but about sequencing.
        // A unmapped read with 0x400 FLAG should be able to kill the mapped mate which may have no 0x400 FLAG.
        if(global_context -> is_duplicate_ignored)
        {
            if(read2.alignment_masks & SAM_FLAG_DUPLICATE)
            {
                thread_context->read_counters[0].unassigned_duplicate ++;
                if(global_context -> detail_output_fp)
                    fprintf(global_context -> detail_output_fp,"%s\tUnassigned_Duplicate\t*\t*\t*\n", read1.read_name);
                return;
            }
        }
        
        if(SAM_FLAG_UNMAPPED & read2.alignment_masks) continue;
        
        char * NH_pos = strstr(tmp_tok_ptr,"\tNH:i:");
        if(NH_pos)
        {
            if(NH_pos[6]>'1' || isdigit(NH_pos[7]))
            {
                if(global_context -> is_multi_mapping_allowed == 0)
                {
                    // now it is a NH>1 read
                    if(!is_second_read)
                        read1.allow_process = 0;
                    else
                        read2.allow_process = 0;
                    if(is_second_read && !read1.allow_process)      // both reads are multimappers!
                    {
                        thread_context->read_counters[0].unassigned_multimapping ++;
                        if(global_context -> detail_output_fp)
                            fprintf(global_context -> detail_output_fp,"%s\tUnassigned_MultiMapping\t*\t*\t*\n", read1.read_name);
                        return;
                    }
                    continue;
                }
            }
        }
    }

    // logic step to determine the read is assigned hierachically or independently.
    for(f_idx = 0; f_idx < global_context->feature_type_num; f_idx++)
    {
        if(read_assign(global_context, thread_context, f_idx, read1, read2))
            if(!global_context->is_independent_assign)
                break;
    }
}

// Thread worker which handles buffer from BAM.
void * feature_count_worker(void * vargs)
{
	void ** args = (void **) vargs;

	fc_thread_global_context_t * global_context = args[0];
	fc_thread_thread_context_t * thread_context = args[1];

	free(vargs);

	if(global_context -> is_SAM_file)
	{
		while (1)
		{
			while(1)
			{
				int is_retrieved = 0;
				pthread_mutex_lock(&thread_context->input_buffer_lock);
				if(thread_context->input_buffer_remainder)
				{
					int is_second_read;
					unsigned int buffer_read_bytes ;
					unsigned int buffer_read_ptr;
					if(thread_context->input_buffer_remainder <= thread_context->input_buffer_write_ptr)
						buffer_read_ptr = thread_context->input_buffer_write_ptr - thread_context->input_buffer_remainder; 
					else
						buffer_read_ptr = thread_context->input_buffer_write_ptr + global_context->input_buffer_max_size - thread_context->input_buffer_remainder;

					for(is_second_read = 0; is_second_read < (global_context->is_paired_end_data? 2:1); is_second_read++)
					{
						char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;
						
						for(buffer_read_bytes=0; ; buffer_read_bytes++)
						{
							char nch =  thread_context->input_buffer[buffer_read_ptr ++];
							curr_line_buff[buffer_read_bytes] = nch;
							if(buffer_read_ptr >= global_context->input_buffer_max_size)
								buffer_read_ptr = 0; 
							if(nch=='\n' || buffer_read_bytes>2998){
								curr_line_buff[buffer_read_bytes+1]=0;
								break;
							}
						}
						thread_context->input_buffer_remainder -= buffer_read_bytes + 1;
					}
					is_retrieved = 1;
				}

				pthread_mutex_unlock(&thread_context->input_buffer_lock);
				if(global_context->is_all_finished && !is_retrieved) return NULL;

				if(is_retrieved) break;
				else
					usleep(tick_time);
			}
			process_line_buffer(global_context, thread_context);
		}
	}
	else
	{	// if is BAM: decompress the chunk and process reads.
    
        /******** Retrieve compressed chunk *********/
        
		char * CDATA[BLOCK_CHUNK_NUM]; // compressed chunk
        int xk1 = 0;
        if(global_context -> multithread_unzipping) {
            for (xk1 = 0; xk1 < BLOCK_CHUNK_NUM; xk1++) {
                CDATA[xk1] = malloc(65537) ;
            }
        }
        
        char * PDATA = malloc(65537 * BLOCK_CHUNK_NUM + 10000); // decompressed chunk
        long current_unprocessed_buf = MAX_UNPROCESSED_NUMBER;
		
		//thread_context -> current_read_length1 = global_context -> read_length;
		//thread_context -> current_read_length2 = global_context -> read_length;
		while(1)
		{
            int PDATA_len = 0;
            // Multithreaded decompressing
            if(global_context -> multithread_unzipping) {
                int CDATA_size[BLOCK_CHUNK_NUM] = {0};
                int chunk_num = 0;
                while(1)
                {
                    int is_retrieved = 0;
                    int check_size = 0;
                    
                    //retrieve the next chunk.
                    
                    pthread_mutex_lock(&thread_context->input_buffer_lock);
                    
                    //SUBREADprintf("IT:re=%i\n", thread_context->input_buffer_remainder);
                    
                    if(thread_context->input_buffer_remainder)
                    {
                        assert(thread_context->input_buffer_remainder>4);
                        unsigned int tail_bytes = global_context->input_buffer_max_size - thread_context -> chunk_read_ptr ;
                        if(tail_bytes<4)
                        {
                            thread_context -> chunk_read_ptr = 0;
                            thread_context -> input_buffer_remainder -= tail_bytes;
                            
                        }
                        else
                        {
                            memcpy(&check_size, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
                            if(check_size==0)
                            {
                                thread_context -> chunk_read_ptr = 0;
                                thread_context -> input_buffer_remainder -= tail_bytes;
                            }
                        }
                        
                        // Read in chunk block, if remainder <=0, done
                        int xk2 = 0;
                        
                        for(xk2 = 0; xk2 < BLOCK_CHUNK_NUM ; xk2++) {
                            if(thread_context -> input_buffer_remainder <= 0) {
                                break;
                            }
                            memcpy(&CDATA_size[xk2], thread_context->input_buffer + thread_context -> chunk_read_ptr, 4);
                            //SUBREADprintf("%i\t", CDATA_size[xk2]);
                            thread_context -> chunk_read_ptr+=4;
                            thread_context -> input_buffer_remainder -= 4;
                            if(CDATA_size[xk2]<0 || CDATA_size[xk2] > 65536)
                            {
                                SUBREADprintf("size = %i\n", CDATA_size[xk2]);
                                SUBREADprintf("THREAD ABNORMALLY QUIT\n");
                                return NULL;
                            }
                            memcpy(CDATA[xk2], thread_context -> input_buffer + thread_context -> chunk_read_ptr , CDATA_size[xk2]);
                            thread_context -> chunk_read_ptr += CDATA_size[xk2];
                            thread_context -> input_buffer_remainder -= CDATA_size[xk2];
                            
                        }
                        
                        chunk_num = xk2;
                        //SUBREADprintf("\ncn = %i\n", chunk_num);
                        
                        if( chunk_num > 0 )
                            is_retrieved = 1;
                    }
                    
                    pthread_mutex_unlock(&thread_context->input_buffer_lock);
                    if(global_context->is_all_finished && !is_retrieved){
                        int xk3 = 0;
                        for (xk3 = 0; xk3 < BLOCK_CHUNK_NUM; xk3++) {
                            free(CDATA[xk3]);
                        }
                        free(PDATA);
                        return NULL;
                    }
                    
                    if(is_retrieved) break;
                    else
                        usleep(tick_time);
                    
                }
                
                /********** Decompress chunk **********/
                
                int chunk_size = 0;
                int i = 0;
                for(i = 0; i < chunk_num; i++) {
                    chunk_size = SamBam_unzip(PDATA + PDATA_len , CDATA[i] , CDATA_size[i]);
                    //SUBREADprintf("%i\t", chunk_size);
                    if(chunk_size<0) {
                        SUBREADprintf("ERROR: BAM GZIP FORMAT ERROR.\n");
                        int j = 0;
                        for (j = 0; j < BLOCK_CHUNK_NUM; j++) {
                            free(CDATA[j]);
                        }
                        free(PDATA);
                        return NULL;
                    }
                    PDATA_len += chunk_size;
                }
            }
            // Mainthread decompressing mode (Default)
            else {
                while(1)
                {
                    int is_retrieved = 0;
                    
                    pthread_mutex_lock(&thread_context->input_buffer_lock);
                    if(thread_context->input_buffer_remainder)
                    {
                        assert(thread_context->input_buffer_remainder>4);
                        unsigned int tail_bytes = global_context->input_buffer_max_size - thread_context -> chunk_read_ptr ;
                        if(tail_bytes<4)
                        {
                            thread_context -> chunk_read_ptr = 0;
                            thread_context -> input_buffer_remainder -= tail_bytes;
                            memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
                        }
                        else
                        {
                            memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
                            if(PDATA_len==0)
                            {
                                thread_context -> chunk_read_ptr = 0;
                                thread_context -> input_buffer_remainder -= tail_bytes;
                                memcpy(&PDATA_len, thread_context->input_buffer , 4);
                            }
                        }
                        thread_context -> chunk_read_ptr+=4;
                        thread_context -> input_buffer_remainder -= 4;
                        
                        if(PDATA_len<0)
                        {
                            SUBREADprintf("THREAD ABNORMALLY QUIT\n");
                            return NULL;
                        }
                        
                        memcpy(PDATA, thread_context -> input_buffer + thread_context -> chunk_read_ptr , PDATA_len);
                        thread_context -> chunk_read_ptr += PDATA_len;
                        thread_context -> input_buffer_remainder -= PDATA_len;
                        
                        if( PDATA_len > 0 )
                            is_retrieved = 1;
                    }
                    
                    pthread_mutex_unlock(&thread_context->input_buffer_lock);
                    if(global_context->is_all_finished && !is_retrieved){
                        free(PDATA);
                        return NULL;
                    }
                    
                    if(is_retrieved) break;
                    else
                        usleep(tick_time);
                }
            }
            
            //SUBREADprintf("\nCHUNK_NUM = %i\tPDATA_LEN = %i\n", chunk_num, PDATA_len);
        
            
            /********** Get reads and call processor **********/
            
            SamBam_Alignment * aln = &thread_context->aln_buffer;
			
			// convert binary reads into sam lines and process;
			int PDATA_ptr = 0;
            char tmp_line[global_context -> line_length + 2];
            int first_read_mark = 1;
        
			while(PDATA_ptr < PDATA_len)
			{
				int is_second_read, binary_read_len = 0;
				for(is_second_read = 0; is_second_read <= global_context -> is_paired_end_data; is_second_read++)
				{
					int local_PDATA_ptr = PDATA_ptr;
					char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;

					memcpy(&binary_read_len, PDATA + PDATA_ptr, 4);
                    //SUBREADprintf("%i\n", binary_read_len);
                    if(binary_read_len > 10000)
                    {
                        SUBREADprintf("   A format error was detected in this BAM file.\n");
                        if (global_context->multithread_unzipping) {
                            SUBREADprintf("   Please rerun without multithreaded decompression.\n");
                        } else
                            SUBREADprintf("   Please check the file format using samtools.\n");
                        return NULL;
                    }
                    
					int ret = PBam_chunk_gets(PDATA, &local_PDATA_ptr, PDATA_len, global_context -> sambam_chro_table, curr_line_buff, 2999, aln,0);
					if(ret<0)
						SUBREADprintf("READ DECODING ERROR!\n");
                    
                    if(global_context -> is_paired_end_data && PDATA_ptr == 0) // paired end mode, store first read of the chunk in case not in pair
                    {
                        sprintf(tmp_line, "%s", thread_context -> line_buffer1);
                    }

					PDATA_ptr += 4 + binary_read_len;
                    
                    if (global_context -> is_paired_end_data && is_second_read == 0 && PDATA_ptr >= PDATA_len) // last read unpaired
                    {
                        first_read_mark = 0; // Even if it is first read, because it is at the end of the chunk, it must be treated as last read unpaired!
                        thread_context -> step_back = 1;    // this is just a mark, step back is not performed
                        break;
                    }
				}
                
                if(!global_context -> is_paired_end_data)       // single-end, directly process
                    process_line_buffer(global_context, thread_context);
                else if((!first_read_mark && is_second_read) || (first_read_mark && !probe_pair_name(thread_context)))  // paired-end, only process the middle ones.
                {
                    probe_pair_name(thread_context);
                    process_line_buffer(global_context, thread_context);
                }
                
                

                // If unpaired and at boundary
                if(thread_context -> step_back)
                {
                    if(is_second_read == 0)         // last read unpaired, put in unprocessed pool
                    {
                        char * last_line = malloc(global_context -> line_length + 2);
                        sprintf(last_line, "%s", thread_context -> line_buffer1);
                        thread_context -> unprocessed_read_ptrs[thread_context -> unprocessed_read_cnts ++] = last_line;
                        if(thread_context -> unprocessed_read_cnts >= current_unprocessed_buf) {
                            //SUBREADprintf("buf = %li\n", current_unprocessed_buf);
                            current_unprocessed_buf = 1.3 * current_unprocessed_buf;
                            thread_context -> unprocessed_read_ptrs = realloc(thread_context -> unprocessed_read_ptrs, current_unprocessed_buf * sizeof(char *));
                        }
                        //printf("last in chunk, buff1: %s\n", thread_context -> line_buffer1);
                    }
                    else
                    {
                        PDATA_ptr -= (4+binary_read_len);
                        thread_context -> step_back = 0;
                        if(first_read_mark)         // first read unpaired
                        {
                            char * first_line = malloc(global_context -> line_length + 2);
                            sprintf(first_line, "%s", tmp_line);
                            thread_context -> unprocessed_read_ptrs[thread_context -> unprocessed_read_cnts ++] = first_line;
                            if(thread_context -> unprocessed_read_cnts >= current_unprocessed_buf) {
                                //SUBREADprintf("buf = %li\n", current_unprocessed_buf);
                                current_unprocessed_buf = 1.3 * current_unprocessed_buf;
                                thread_context -> unprocessed_read_ptrs = realloc(thread_context -> unprocessed_read_ptrs, current_unprocessed_buf* sizeof(char *));
                            }
                            //printf("first in chunk, buff1: %s\n", tmp_line);
                        }
                        else
                            thread_context -> missing_mates++;
                    }
                }
                first_read_mark = 0;
			}
            thread_context -> step_back = 0;    // reset step_back after the chunk is done
		}
	}
}

// In the bulk process, the first and last unpaired read in BAM chunks won't be processed.
// They are pooled and sorted here after all threads finishes.
void combine_sort_unprocessed(fc_thread_global_context_t * global_context)
{
    int xk1, xk2, up = 0;
    long long int total_unprocessed = 0;
    char * line = NULL;
    
    // get the total number of unassigned reads
    for(xk1=0; xk1<global_context-> thread_number; xk1++)
        total_unprocessed += global_context -> thread_contexts[xk1].unprocessed_read_cnts;
    global_context -> unprocessed_cnts = total_unprocessed;
    
    // combine the thread arrays to one large array and get names out
    char ** all_unprocessed_ptrs = malloc(sizeof(char*) * (total_unprocessed + 1));
    char * read_names[total_unprocessed + 1];
    for(xk1=0; xk1<global_context-> thread_number; xk1++) {
        for(xk2 = 0; xk2 < global_context -> thread_contexts[xk1].unprocessed_read_cnts; xk2++)
        {
            line = global_context -> thread_contexts[xk1].unprocessed_read_ptrs[xk2];
            all_unprocessed_ptrs[up] = line;
            read_names[up++] = vs_tab_strtok(line);
        }
        free(global_context -> thread_contexts[xk1].unprocessed_read_ptrs);
    }
    
    // sort the reads by name
    void * in_array[2];
    in_array[0] = read_names;
    in_array[1] = all_unprocessed_ptrs;
    merge_sort(in_array, total_unprocessed, 0, read_sort_compare, read_sort_exchange, read_sort_merge);
    global_context -> unprocessed_reads = all_unprocessed_ptrs;
    for(xk1=0; xk1 < total_unprocessed; xk1++)
    {
        free(read_names[xk1]);
    }
    //SUBREADprintf("tn = %lli\n", total_unprocessed);
}

// Merge read assignment results from all threads.
// Note if it is BAM file this function will be called twice. The second time is to add results from the unprocessed_read handling thread to the bulk results.
void fc_thread_merge_results(fc_thread_global_context_t * global_context, unsigned int ** nreads)
{
	int xk1, xk2, f_idx;

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
        for(f_idx = 0; f_idx < global_context -> feature_type_num; f_idx++)
        {
            for(xk2=0; xk2<global_context -> exontable_exons[f_idx]; xk2++)
            {
                nreads[f_idx][xk2]+=global_context -> thread_contexts[xk1].count_table[f_idx][xk2];
            }
            if(f_idx == 0)
            {
                global_context->all_reads += global_context -> thread_contexts[xk1].all_reads;
                global_context->missing_mates += global_context -> thread_contexts[xk1].missing_mates;
            }
            
            global_context -> read_counters[f_idx].unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_ambiguous;
            global_context -> read_counters[f_idx].unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_nofeatures;
            global_context -> read_counters[f_idx].unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_unmapped;
            global_context -> read_counters[f_idx].unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_mappingquality;
            global_context -> read_counters[f_idx].unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_fragmentlength;
            global_context -> read_counters[f_idx].unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_chimericreads;
            global_context -> read_counters[f_idx].unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_multimapping;
            global_context -> read_counters[f_idx].unassigned_nonjunction += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_nonjunction;
            global_context -> read_counters[f_idx].unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_duplicate;
            global_context -> read_counters[f_idx].assigned_reads += global_context -> thread_contexts[xk1].read_counters[f_idx].assigned_reads;
        }
	}

    if(!global_context->is_paired_end_data || global_context->is_really_finished)
    {
        char pct_str[10];
        if(global_context->all_reads>0)
        {
            for(f_idx = 0; f_idx < global_context -> feature_type_num; f_idx++)
            {
                global_context -> total_assigned += global_context -> read_counters[f_idx].assigned_reads;
            }
            sprintf(pct_str,"(%.1f%%%%)", global_context -> total_assigned * 100./global_context->all_reads);
        }
        else
            pct_str[0]=0;
        
        if (global_context -> is_paired_end_data && global_context->missing_mates)
            print_in_box(80,0,0,"   Read pairs with missing mates : %llu", global_context->missing_mates);
        print_in_box(80,0,0,"   Total %s : %llu", global_context -> is_paired_end_data?"read pairs":"reads", global_context->all_reads);
        print_in_box(pct_str[0]?81:80,0,0,"   Successfully assigned %s : %llu %s", global_context -> is_paired_end_data?"read pairs":"reads", global_context -> total_assigned, pct_str);
        print_in_box(80,0,0,"   Running time : %.2f minutes", (miltime() - global_context -> start_time)/60);
        print_in_box(80,0,0,"");
    }
}

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_strand_checked, char * output_fname, int is_detail_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column, int min_map_qual_score, int is_multi_mapping_allowed, int is_SAM, char * cmd_rebuilt, int is_input_file_resort_needed, int feature_block_size, int isCVersion, int fiveEndExtension,  int threeEndExtension, int minReadLength, int minReadOverlap, int maxReadNonoverlap, int is_split_alignments_only, int reduce_5_3_ends_to_one, char * debug_command, int is_duplicate_ignored, int run_mode, int minDifAmbiguous, int featureTypeNum, char ** featureTypeList, int isindependentAssign, int nonempty_modified, int multithreadunzip)
{

	global_context -> input_buffer_max_size = buffer_size;
	global_context -> all_reads = 0;
    global_context -> missing_mates = 0;
    global_context -> unprocessed_cnts = 0;
    global_context -> total_assigned = 0;
	global_context -> detail_output_fp = NULL;

	global_context -> isCVersion = isCVersion;
	global_context -> is_read_details_out = is_detail_out;
	global_context -> is_paired_end_data = is_PE_data;
	global_context -> is_strand_checked = is_strand_checked;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	global_context -> is_multi_mapping_allowed = is_multi_mapping_allowed;
	global_context -> is_split_alignments_only = is_split_alignments_only;
	global_context -> is_duplicate_ignored = is_duplicate_ignored;
	global_context -> reduce_5_3_ends_to_one = reduce_5_3_ends_to_one;
	global_context -> is_SAM_file = is_SAM;
    global_context -> multithread_unzipping = multithreadunzip;

	global_context -> thread_number = threads;
	global_context -> min_mapping_quality_score = min_map_qual_score;

	global_context -> cmd_rebuilt = cmd_rebuilt;
	global_context -> is_input_file_resort_needed = is_input_file_resort_needed;
	global_context -> feature_block_size = feature_block_size;
	global_context -> five_end_extension = fiveEndExtension;
	global_context -> three_end_extension = threeEndExtension;
    global_context -> read_M_len_required = minReadLength;
	global_context -> overlap_length_required = minReadOverlap;
    global_context -> nonoverlap_length_allowed = maxReadNonoverlap;
	global_context -> debug_command = debug_command;
    global_context -> min_dif_ambiguous = minDifAmbiguous;

	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);
	strcpy(global_context -> output_file_name, output_fname);

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;
    
    global_context -> run_mode = run_mode;
    global_context -> is_independent_assign = isindependentAssign;
    global_context -> is_nonempty_modified = nonempty_modified;
    
    global_context -> feature_type_num = featureTypeNum;
    global_context -> feature_type_list = featureTypeList;
    
    global_context -> longest_chro_name = malloc(sizeof(int) * featureTypeNum);
    global_context -> exontable_nchrs = malloc(sizeof(int) * featureTypeNum);
    global_context -> gene_name_array = malloc(sizeof(unsigned char **) * featureTypeNum);
    global_context -> gene_name_table = malloc(sizeof(HashTable *) * featureTypeNum);
    global_context -> exontable_chro_table = malloc(sizeof(HashTable *) * featureTypeNum);
    global_context -> original_exons = malloc(sizeof(int) * featureTypeNum);
    global_context -> exontable_exons = malloc(sizeof(int) * featureTypeNum);
    global_context -> read_counters = malloc(sizeof(fc_read_counters) * featureTypeNum);
    global_context -> unistr_buffer_space = malloc(sizeof(char *) * featureTypeNum);
    global_context -> unistr_buffer_used = malloc(sizeof(int) * featureTypeNum);
    global_context -> unistr_buffer_size = malloc(sizeof(int) * featureTypeNum);
    
    int f_idx;
    for(f_idx = 0; f_idx < featureTypeNum; f_idx++)
    {
        global_context -> unistr_buffer_size[f_idx] = 1024*1024*2;
        global_context -> unistr_buffer_space[f_idx] = malloc(global_context -> unistr_buffer_size[f_idx]);
        global_context -> unistr_buffer_used[f_idx] = 0;
        global_context -> read_counters[f_idx].unassigned_ambiguous=0;
        global_context -> read_counters[f_idx].unassigned_nofeatures=0;
        global_context -> read_counters[f_idx].unassigned_unmapped=0;
        global_context -> read_counters[f_idx].unassigned_mappingquality=0;
        global_context -> read_counters[f_idx].unassigned_fragmentlength=0;
        global_context -> read_counters[f_idx].unassigned_chimericreads=0;
        global_context -> read_counters[f_idx].unassigned_multimapping=0;
        global_context -> read_counters[f_idx].unassigned_secondary=0;
        global_context -> read_counters[f_idx].unassigned_nonjunction=0;
        global_context -> read_counters[f_idx].unassigned_duplicate=0;
        global_context -> read_counters[f_idx].assigned_reads=0;
    }
}

int fc_thread_start_threads(fc_thread_global_context_t * global_context, int * et_exons, int ** et_geneid, long ** et_start, long ** et_stop, unsigned char ** et_strand, long ** et_bk_end_index, long ** et_bk_min_start, long ** et_bk_max_end, int read_length)
{
	int xk1, f_idx;
	global_context -> read_length = read_length;
	global_context -> exontable_geneid = et_geneid;
	global_context -> exontable_start = et_start;
	global_context -> exontable_stop = et_stop;
	global_context -> exontable_strand = (char **)et_strand;
	global_context -> exontable_block_end_index = et_bk_end_index;
	global_context -> exontable_block_max_end = et_bk_max_end;
	global_context -> exontable_block_min_start = et_bk_min_start;
    if(!global_context -> unprocessed_cnts)     // first bulk process
    {
        global_context -> is_all_finished = 0;
        global_context -> is_unpaired_warning_shown = 0;
    }
    else        // bulk finished, deal with unprocessed
        global_context -> thread_number = 1;
    global_context -> is_really_finished = 0;
	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * global_context -> thread_number);
	for(xk1=0; xk1<global_context -> thread_number; xk1++)
	{
	//	printf("CHRR_MALLOC\n");
		pthread_mutex_init(&global_context->thread_contexts[xk1].input_buffer_lock, NULL);
        global_context -> thread_contexts[xk1].step_back = 0;
		global_context -> thread_contexts[xk1].input_buffer_remainder = 0;
		global_context -> thread_contexts[xk1].input_buffer_write_ptr = 0;
		global_context -> thread_contexts[xk1].input_buffer = malloc(global_context -> input_buffer_max_size);
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].chunk_read_ptr = 0;
        global_context -> thread_contexts[xk1].missing_mates = 0;
		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].line_buffer1 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].line_buffer2 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].chro_name_buff = malloc(CHROMOSOME_NAME_LENGTH);
		global_context -> thread_contexts[xk1].strm_buffer = malloc(sizeof(z_stream));
        global_context -> thread_contexts[xk1].unprocessed_read_cnts = 0;
        global_context -> thread_contexts[xk1].count_table = malloc(sizeof(unsigned int *) * global_context -> feature_type_num);
        global_context -> thread_contexts[xk1].read_counters = malloc(sizeof(fc_read_counters) * global_context -> feature_type_num);
        
        if(!global_context -> unprocessed_cnts)
            global_context -> thread_contexts[xk1].unprocessed_read_ptrs = malloc(sizeof(char * ) * MAX_UNPROCESSED_NUMBER);

        for(f_idx = 0; f_idx < global_context -> feature_type_num; f_idx++)
        {
            global_context -> thread_contexts[xk1].count_table[f_idx] = calloc(sizeof(unsigned int), et_exons[f_idx]);
            if(!global_context -> thread_contexts[xk1].count_table[f_idx]) return 1;
            
            global_context -> thread_contexts[xk1].read_counters[f_idx].assigned_reads = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_ambiguous = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_nofeatures = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_unmapped = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_mappingquality = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_fragmentlength = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_chimericreads = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_multimapping = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_secondary = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_nonjunction = 0;
            global_context -> thread_contexts[xk1].read_counters[f_idx].unassigned_duplicate = 0;
        }

		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];

		if((global_context -> thread_number>1 || ! global_context -> is_SAM_file) && !global_context->is_all_finished)
			pthread_create(&global_context -> thread_contexts[xk1].thread_object, NULL, feature_count_worker, thread_args);
	}
	return 0;
}

void fc_thread_destroy_thread_context(fc_thread_global_context_t * global_context)
{
	int xk1, f_idx;
	if(global_context -> is_read_details_out && global_context -> is_really_finished)
	{
		fclose(global_context -> detail_output_fp);
		global_context -> detail_output_fp = NULL;
	}

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		//printf("CHRR_FREE\n");
        for(f_idx = 0; f_idx < global_context -> feature_type_num; f_idx++)
        {
            free(global_context -> thread_contexts[xk1].count_table[f_idx]);
        }
        free(global_context -> thread_contexts[xk1].count_table);
        free(global_context -> thread_contexts[xk1].read_counters);
		free(global_context -> thread_contexts[xk1].line_buffer1);
		free(global_context -> thread_contexts[xk1].line_buffer2);	
		free(global_context -> thread_contexts[xk1].input_buffer);
		free(global_context -> thread_contexts[xk1].chro_name_buff);
		free(global_context -> thread_contexts[xk1].strm_buffer);
		pthread_mutex_destroy(&global_context -> thread_contexts[xk1].input_buffer_lock);
	}
	free(global_context -> thread_contexts);
}

void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
		pthread_join(global_context -> thread_contexts[xk1].thread_object, NULL);
}

// This function is called only when the input BAM/SAM is not name sorted (the user must specify -S to sort).
int resort_input_file(fc_thread_global_context_t * global_context)
{
	char * temp_file_name = malloc(300), * fline = malloc(3000);
	SamBam_FILE * sambam_reader;

    print_in_box(80,0,0,"   Resorting the input file ...");
	sprintf(temp_file_name, "./temp-core-%06u-%08X.sam", getpid(), rand());
	sambam_reader = SamBam_fopen(global_context-> input_file_name, global_context-> is_SAM_file?SAMBAM_FILE_SAM:SAMBAM_FILE_BAM);

	if(!sambam_reader){
		SUBREADprintf("Unable to open %s.\n", global_context-> input_file_name);
		return -1;
	}

	SAM_sort_writer writer;
	int ret = sort_SAM_create(&writer, temp_file_name, ".");
	if(ret)
	{
		SUBREADprintf("Unable to sort input file because temporary file '%s' cannot be created.\n", temp_file_name);
		return -1;
	}
	int is_read_len_warned = 0;

	while(1)
	{
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) break;
		int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
		if(ret<0) 
		{
			if(!is_read_len_warned)
				SUBREADprintf("WARNING: reads with very long names were found.\n");
			is_read_len_warned = 1;
		//	break;
		}
	//printf("N1=%llu\n",  writer.unpaired_reads);
	}

	sort_SAM_finalise(&writer);
	print_in_box(80,0,0,"   Input was converted to a format accepted by VERSE.");

	SamBam_fclose(sambam_reader);
	strcpy(global_context-> input_file_name, temp_file_name);
	global_context->is_SAM_file = 1;
	free(temp_file_name);
	free(fline);
	return 0;
}

// Write final counts to outname.feature_type.txt. Output gene length (total exon/intron/... length of each gene) if -l is specified.
void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, long * et_start, long * et_stop, unsigned char * et_strand, const char * out_file, int features, unsigned int * column_numbers, char * file_used, int f_length, int f_idx)
{
	int xk1;
	int genes = global_context -> gene_name_table[f_idx] -> numOfElements;
	unsigned int * gene_columns = calloc(sizeof(unsigned int) , genes);
    unsigned long * gene_lengths = calloc(sizeof(unsigned long) , genes);
    
	FILE * fp_out = fopen(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
		return;
	}
	char * tmp_ptr = NULL, * next_fn;
    fprintf(fp_out,"gene\tcount");
    if(f_length)
        fprintf(fp_out,"\tlength");
    fprintf(fp_out, "\n");
    
	next_fn = strtok_r(file_used, ";", &tmp_ptr);
    
    for(xk1 = 0; xk1 < features; xk1++)
    {
        int gene_id = et_geneid[xk1];
        gene_columns[gene_id] += column_numbers[xk1];
        if(f_length)
            gene_lengths[gene_id] += (et_stop[xk1] - et_start[xk1] + 1);
    }
    
    for(xk1 = 0 ; xk1 < genes; xk1++)
    {
        unsigned char * gene_symbol = global_context -> gene_name_array[f_idx][xk1];
        fprintf(fp_out, "%s\t%u", gene_symbol, gene_columns[xk1]);
        if(f_length)
            fprintf(fp_out, "\t%ld", gene_lengths[xk1]);
        fprintf(fp_out,"\n");
    }
    
    free(gene_columns);
    free(gene_lengths);
    fclose(fp_out);
}

// Write stats to outname.summary.txt.
void fc_write_final_counts(fc_thread_global_context_t * global_context, const char * out_file, char * file_used, char * annot)
{
	char fname[350];

    if(global_context->feature_type_num > 1)        // one single summary file if Hierarchical assign
        sprintf(fname, "%s.summary.txt", out_file);
    else
        sprintf(fname, "%s.%s.summary.txt", out_file, global_context->feature_type_list[0]);        // this is to prevent overwriting if run separately for each feature type
	FILE * fp_out = fopen(fname,"w");

	if(!fp_out){
		SUBREADprintf("Unable to create summary file %s.\n", fname);
		return;
	}
    
    int f_idx;
    fprintf(fp_out,"AnnotationFile\t%s\n", annot);
    for(f_idx = 0; f_idx < global_context->feature_type_num; f_idx++)
    {
        fputc(toupper(global_context->feature_type_list[f_idx][0]), fp_out);
        fprintf(fp_out,"%s:\n", global_context->feature_type_list[f_idx]+1);
        fprintf(fp_out,"\tUnmerged\t%u\n", global_context -> original_exons[f_idx]);
        fprintf(fp_out,"\tMerged\t%u\n", global_context -> exontable_exons[f_idx]);
        fprintf(fp_out,"\tGenes\t%ld\n", global_context -> gene_name_table[f_idx] -> numOfElements);
        fprintf(fp_out,"\tChromosomes\t%d\n", global_context -> exontable_nchrs[f_idx]);
    }
    
	fprintf(fp_out,"\nInputFile\t%s\n", file_used);
    fprintf(fp_out,"Total%s\t%llu\n", global_context -> is_paired_end_data?"ReadPairs":"Reads", global_context->all_reads);
    fprintf(fp_out,"MissingMates\t%llu\n", global_context->missing_mates);
    fprintf(fp_out,"TotalAssigned\t%llu\n", global_context->total_assigned);
    fprintf(fp_out,"TotalAssignedFraction\t%.1f%%\n\n", (global_context->total_assigned)*100./(global_context->all_reads));
    
    for(f_idx = 0; f_idx < global_context->feature_type_num; f_idx++)
    {
        fprintf(fp_out,"Assigned%c%s\t%llu\n", toupper(global_context->feature_type_list[f_idx][0]), global_context->feature_type_list[f_idx]+1, global_context -> read_counters[f_idx].assigned_reads);
        fprintf(fp_out,"Assigned%c%sFraction\t%.1f%%\n", toupper(global_context->feature_type_list[f_idx][0]), global_context->feature_type_list[f_idx]+1, (global_context->read_counters[f_idx].assigned_reads)*100./(global_context->all_reads));
        fprintf(fp_out,"Ambiguous%c%s\t%llu\n", toupper(global_context->feature_type_list[f_idx][0]), global_context->feature_type_list[f_idx]+1, global_context -> read_counters[f_idx].unassigned_ambiguous);
        if(global_context -> is_independent_assign)
            fprintf(fp_out,"NoFeature%c%s\t%llu\n", toupper(global_context->feature_type_list[f_idx][0]), global_context->feature_type_list[f_idx]+1, global_context -> read_counters[f_idx].unassigned_nofeatures);
        fprintf(fp_out,"\n");
    }
    
    if(!global_context -> is_independent_assign)
        fprintf(fp_out,"NoFeature\t%llu\n", global_context -> read_counters[global_context->feature_type_num - 1].unassigned_nofeatures);
    
    fprintf(fp_out,"MultiMapping\t%llu\n", global_context -> read_counters[0].unassigned_multimapping);
    fprintf(fp_out,"Unmapped\t%llu\n", global_context -> read_counters[0].unassigned_unmapped);
    fprintf(fp_out,"LowMAPQ\t%llu\n", global_context -> read_counters[0].unassigned_mappingquality);
    fprintf(fp_out,"TLENTooLong/TooShort\t%llu\n", global_context -> read_counters[0].unassigned_fragmentlength);
    fprintf(fp_out,"Chimera\t%llu\n", global_context -> read_counters[0].unassigned_chimericreads);
    fprintf(fp_out,"Nonjunction\t%llu\n", global_context -> read_counters[0].unassigned_nonjunction);
    fprintf(fp_out,"Duplicate\t%llu\n", global_context -> read_counters[0].unassigned_duplicate);
    
	fclose(fp_out);
}

static struct option long_options[] =
{
    {"singleEnd",no_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension3", required_argument, 0, 0},
	{"read2pos", required_argument, 0, 0},
	{"minReadOverlap", required_argument, 0, 0},
    {"maxReadNonoverlap", required_argument, 0, 0},
	{"countSplitAlignmentsOnly", no_argument, 0, 0},
	{"debugCommand", required_argument, 0, 0},
	{"ignoreDup", no_argument, 0, 0},
    {"minDifAmbiguous", required_argument, 0, 0},
    {"assignIndependently",no_argument, 0, 0},
    {"nonemptyModified",no_argument, 0, 0},
    {"multithreadDecompress",no_argument, 0, 0},
	{0, 0, 0, 0}
};

void print_usage()
{
	VERSEprintf("\nVersion %s\n\n", VERSE_VERSION);

	VERSEputs("\nUsage: verse [options] -a <annotation_file> -o <output_file> input_file\n");
	VERSEputs("    Required parameters:\n");
	VERSEputs("    -a <input>\tGive the name of the annotation file. The program currently only");
	VERSEputs("              \tsupports GTF format.");
	VERSEputs("    "); 
	VERSEputs("    -o <input>\tGive the general name of the output file, e.g., 'Sample_A'.");
	VERSEputs("              \tThe summary file will be named 'Sample_A.summary.txt.'");
    VERSEputs("              \tThe file containing gene counts will be named 'Sample_A.exon.txt',");
    VERSEputs("              \t'Sample_A.intron.txt', etc.");
	VERSEputs("    "); 
	VERSEputs("    input_file\tGive the name of input read file that contains the read mapping");
	VERSEputs("              \tresults. Format of input file is automatically determined (SAM/BAM)");
	VERSEputs("              \tIf it is paired-end data, the file MUST be name-sorted, otherwise");
	VERSEputs("              \tthe user MUST specify -S, or use samtools to sort it by name.");
	VERSEputs("    "); 
	VERSEputs("    Optional parameters:"); 
	VERSEputs("    ");
    VERSEputs("    -z <int>  \tThe Running Mode: 0 by default (featureCounts), 1 (HTSeq-Union),");
    VERSEputs("              \t 2 (HTSeq-Intersection_strict), 3 (HTSeq-Intersection_nonempty),");
    VERSEputs("              \t 4 (Union_strict), 5 (Cover_length).");
    VERSEputs("              \tPlease refer to the manual or use '-Z' to check the details of");
    VERSEputs("              \teach mode.");
    VERSEputs("    ");
	VERSEputs("    -t <input>\tSpecify the feature type. Only rows which have the matched"); 
	VERSEputs("              \tfeature type in the provided GTF annotation file will be included");
	VERSEputs("              \tfor read counting. 'exon' by default.");
    VERSEputs("              \tProviding a list of feature types (e.g., -t 'exon;intron;mito')");
    VERSEputs("              \twill automatically enter hierarchical_assign mode. If the user");
    VERSEputs("              \twants to assign independently for each feature type, please");
    VERSEputs("              \tspecify '--assignIndependently'. Use -Z to check the details.");
	VERSEputs("    "); 
	VERSEputs("    -g <input>\tSpecify the gene_identifier attribute, which is normally 'gene_id'");
	VERSEputs("              \tor 'gene_name'. 'gene_id' by default.");
	VERSEputs("    ");
    VERSEputs("    -S        \tIf the input file is paired-end data but not sorted by read name,");
    VERSEputs("              \tthis option MUST be specified.");
    VERSEputs("    ");
	VERSEputs("    -s <int>  \tIndicate if strand-specific read counting should be performed.");
	VERSEputs("              \tIt has three possible values:  0 (unstranded), 1 (stranded) and");
	VERSEputs("              \t2 (reversely stranded). 0 by default.");
	VERSEputs("    ");
	VERSEputs("    -Q <int>  \tThe minimum mapping quality score a read must satisfy in order");
	VERSEputs("              \tto be counted. For paired-end reads, at least one end should");
	VERSEputs("              \tsatisfy this criteria. 0 by default."); 
	VERSEputs("    ");
    VERSEputs("    -l        \tOutput the gene length. This length is feature_type specific,");
    VERSEputs("              \twhich means the length in an exon_count file will be the total");
    VERSEputs("              \texon length of each gene, the length in an intron_count file will");
    VERSEputs("              \tbe the total intron length, etc.");
    VERSEputs("    ");
	VERSEputs("    -R        \tOutput read assignment details for each read/read pairs. The");
	VERSEputs("              \tdetails will be saved to a tab-delimited file that contains five");
	VERSEputs("              \tcolumns including read name, status(assigned or the reason if not");
	VERSEputs("              \tassigned), feature type and assigned gene name.");
	VERSEputs("    ");
    VERSEputs("    -T <int>  \tNumber of threads used for read assignment. 1 by default.");
    VERSEputs("              \tNote that when running, VERSE will initiate one main thread and");
    VERSEputs("              \tat least one helper thread for read assignment. This option");
    VERSEputs("              \tspecifies the number of helper threads.");
    VERSEputs("    ");
    VERSEputs("    --singleEnd                 If specified, VERSE will assume the input is");
    VERSEputs("              \tsingle-end data, and assign every reads individually. If this is");
    VERSEputs("              \tnot specified(default), the input will be treated as paired-end");
    VERSEputs("              \tdata. The 'missing mate' count will show how many reads are not");
    VERSEputs("              \tin pairs. VERSE can correctly assign single-end data in paired-");
    VERSEputs("              \tend mode, but it will take longer and the reads will all be");
    VERSEputs("              \tcounted as missing mates. So it is recommended to specify this");
    VERSEputs("              \tif the user knows it is single-end.");
    VERSEputs("    ");
    VERSEputs("    --assignIndependently       If specified, VERSE will assign reads to listed");
    VERSEputs("              \tfeature types independently. This has the same effect as running");
    VERSEputs("              \tVERSE separately for each feature type, except that it only");
    VERSEputs("              \trequires one run, thus is more efficient.");
    VERSEputs("    ");
	VERSEputs("    --readExtension5 <int>      Reads are extended upstream by <int> bases from");
	VERSEputs("              \ttheir 5' end."); 
	VERSEputs("    "); 
	VERSEputs("    --readExtension3 <int>      Reads are extended upstream by <int> bases from");
	VERSEputs("              \ttheir 3' end."); 
	VERSEputs("    ");
	VERSEputs("    --minReadOverlap <int>      Specify the minimum number of overlapped bases");
	VERSEputs("              \trequired to assign a read to a feature. 1 by default. ");
	VERSEputs("    ");
    VERSEputs("    --maxReadNonoverlap <int>   Specify the maximum number of non-overlapped bases");
    VERSEputs("              \ta read can have. A read is considered no_feature if its number of");
    VERSEputs("              \tnon-overlapped bases exceeds this threshold. This option is not");
    VERSEputs("              \tuseful in strict mode since it requires the assigned feature to");
    VERSEputs("              \toverlap every base of the read.");
    VERSEputs("    ");
	VERSEputs("    --countSplitAlignmentsOnly  If specified, only split alignments (CIGAR");
	VERSEputs("              \tstrings containing letter `N') will be counted. All the other");
	VERSEputs("              \talignments will be ignored. An example of split alignments is");
	VERSEputs("              \tthe exon-spanning reads in RNA-seq data.");
	VERSEputs("    "); 
	VERSEputs("    --read2pos <5:3>            The read is reduced to its 5' most base or 3'");
	VERSEputs("              \tmost base. Read summarization is then performed based on the");
	VERSEputs("              \tsingle base which the read is reduced to."); 
	VERSEputs("    "); 
	VERSEputs("    --ignoreDup                 If specified, reads that were marked as");
	VERSEputs("              \tduplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM");
	VERSEputs("              \tfile is used for identifying duplicate reads. In paired end");
	VERSEputs("              \tdata, the entire read pair will be ignored if at least one end");
	VERSEputs("              \tis found to be a duplicate read.");
	VERSEputs("    ");
    VERSEputs("    --minDifAmbiguous <int>     This option can only be used in VERSE-cover_length");
    VERSEputs("              \tmode. When the read or the read pair hits more than one genes,");
    VERSEputs("              \tthe number of bases overlapped by each gene will be calculated.");
    VERSEputs("              \tIf the covered_length difference between the most covered gene");
    VERSEputs("              \tand the second most covered gene is smaller than this specified");
    VERSEputs("              \tvalue, the read will be set ambiguous. 1 base difference by default.");
    VERSEputs("    ");
    VERSEputs("    --nonemptyModified          This option can only be used in intersection_");
    VERSEputs("              \tnonempty mode. In cases where HTSeq would assign multi-hit reads");
    VERSEputs("              \tto no_feature, VERSE will assign those to ambiguous.");
    VERSEputs("    ");
    VERSEputs("    --multithreadDecompress     BAM files generated with SAMTools or STAR after");
    VERSEputs("              \tyear 2013 has a slight format change which allows multithreaded");
    VERSEputs("              \tdecompression. BAM processing will be faster if this option is");
    VERSEputs("              \tspecifed and multiple core is used.");
    VERSEputs("    ");
	VERSEputs("    Optional paired-end parameters:");
	VERSEputs("    -P        \tIf specified, template length (TLEN in SAM specification) will be");
	VERSEputs("              \tchecked when assigning read pairs (templates) to genes. This option");
	VERSEputs("              \tis only applicable in paired-end mode. The distance thresholds");
	VERSEputs("              \tshould be specified using -d and -D options.");
	VERSEputs("    "); 
	VERSEputs("    -d <int>  \tMinimum template(read pair) length, 50 by default.");
	VERSEputs("    "); 
	VERSEputs("    -D <int>  \tMaximum template(read pair) length, 600 by default.");
	VERSEputs("    "); 
	VERSEputs("    -B        \tIf specified, only read pairs that have both ends successfully");
	VERSEputs("              \taligned will be considered for summarization. This option is only");
	VERSEputs("              \tapplicable for paired-end reads.");
	VERSEputs("    "); 
	VERSEputs("    -C        \tIf specified, the chimeric read pairs (those that have their two");
	VERSEputs("              \tends aligned to different chromosomes) will NOT be included for");
	VERSEputs("              \tsummarization. This option is only applicable for paired-end data.");
	VERSEputs("    ");
    VERSEputs("    ");
    VERSEputs("    Additional Information:");
    VERSEputs("    ");
	VERSEputs("    -v        \tOutput version of the program.");
	VERSEputs("    ");
    VERSEputs("    -Z        \tShow details about the running mode or scheme.");
    VERSEputs("    ");
}

void print_run_mode()
{
    VERSEprintf("\nVersion %s\n\n", VERSE_VERSION);
    VERSEputs("    Running Mode Information:");
    VERSEputs("    ");
    VERSEputs("    The main differences between modes are the stringency of assignment and the method");
    VERSEputs("    of assigning reads that overlap multiple features.");
    VERSEputs("    ");
    VERSEputs("    ");
    VERSEputs("    featureCounts Mode:");
    VERSEputs("    ");
    VERSEputs("    -z 0   featureCounts(default): Quantify by overlapping and voting.");
    VERSEputs("               If the read pair overlaps multiple genes, it will assign the read pair");
    VERSEputs("               to the gene that is overlapped by both reads.");
    VERSEputs("    ");
    VERSEputs("               Please refer to the SUBREAD users guide for more information. ");
    VERSEputs("               This mode is developed by Drs. Yang Liao, Gordon K Smyth and Wei Shi.");
    VERSEputs("               http://bioinf.wehi.edu.au/featureCounts.");
    VERSEputs("    ");
    VERSEputs("    ");
    VERSEputs("    HTseq Modes:");
    VERSEputs("    ");
    VERSEputs("    -z 1   HTSeq-Union: Quantify by overlapping.");
    VERSEputs("               If the read (or read pair) overlaps multiple genes, it will be set");
    VERSEputs("               ambiguous. Only reads that overlap one gene will be assigned.");
    VERSEputs("    ");
    VERSEputs("    -z 2   HTSeq-Intersection_strict: Quantify by overlapping and intersection.");
    VERSEputs("               This mode requires the assigned gene to cover every base of the read.");
    VERSEputs("               If more than one such genes exist,  the read is set ambiguous.");
    VERSEputs("    ");
    VERSEputs("    -z 3   HTSeq-Intersection_nonempty: Quantify by overlapping and intersection.");
    VERSEputs("               This mode does NOT require the assigned gene to cover every base of");
    VERSEputs("               the read, but the gene must cover all sections of the read that overlap");
    VERSEputs("               genes. If more than one such genes exist, the read is set ambiguous.");
    VERSEputs("    ");
    VERSEputs("               Please refer to the HTSeq documentation for more information.");
    VERSEputs("               HTSeq is developed by Dr. Simon Anders at EMBL Heidelberg.");
    VERSEputs("               http://www-huber.embl.de/HTSeq/doc/count.html#count.");
    VERSEputs("               The actual implementation of the HTseq scheme is different in VERSE.");
    VERSEputs("    ");
    VERSEputs("    ");
    VERSEputs("    VERSE Modes:");
    VERSEputs("    ");
    VERSEputs("    -z 4   Union_strict: A combination of HTSeq-Union and HTSeq-Intersection_strict.");
    VERSEputs("               This mode requires every base of the read to overlap one and only one gene.");
    VERSEputs("               This mode is the most conservative.");
    VERSEputs("    ");
    VERSEputs("    -z 5   Cover_length: Quantify by overlapping length comparison.");
    VERSEputs("               After getting a list of overlapping genes, VERSE will calculate the");
    VERSEputs("               overlapping length of each gene and assign the read to the most covered");
    VERSEputs("               gene. One can use --minDifAmbiguous to set the minimum allowed coverage");
    VERSEputs("               difference between the most covered gene and the gene with the second");
    VERSEputs("               highest coverage. If the difference is not large enough the read will");
    VERSEputs("               be set ambiguous.");
    VERSEputs("    ");
    VERSEputs("    ");
    VERSEputs("    Running Scheme Information:");
    VERSEputs("    ");
    VERSEputs("    Hierarchical_assign:");
    VERSEputs("               VERSE offers a hierarchical_assign scheme, which ensures that each read");
    VERSEputs("               is assigned only once to the most interested feature type. When user");
    VERSEputs("               inputs a list of feature types, e.g., 'exon;intron;mito', VERSE will");
    VERSEputs("               assume the first feature type (exon) in the list is of the highest priority,");
    VERSEputs("               and will try assigning reads to it first, using the specified mode (e.g.,");
    VERSEputs("               union mode). The reads that cannot be assigned to exon will enter the");
    VERSEputs("               next round of assignment to the next feature type, i.e., intron.");
    VERSEputs("               This procedure will repeat until every read is successfully assigned,");
    VERSEputs("               or fails to be assigned to any of the feature types.");
    VERSEputs("    ");
    VERSEputs("    Independent_assign:");
    VERSEputs("               If '--assignIndependently' is specified, instead of running hierarchical");
    VERSEputs("               _assign, VERSE will assign reads to feature types independently.");
    VERSEputs("               This has the same effect as running VERSE separately for each feature type, ");
    VERSEputs("               except that it only requires one run, thus is more efficient.");
    VERSEputs("    ");
}

void write_compressed_block(fc_thread_global_context_t * global_context, FILE * fp_in, char ** chunk_in_buff){
    while (1){
        unsigned int real_len = 0;
        int current_thread_id = 0;
        int cdata_size[BLOCK_CHUNK_NUM] = {0};
        int total_size = 0;
        int chunk_num = 0, i = 0;
        
        for (i = 0; i < BLOCK_CHUNK_NUM; i++) {
            cdata_size[i] = PBam_get_next_zchunk(fp_in, chunk_in_buff[i], 65537, & real_len);
            //SUBREADprintf("cs=%i; ", cdata_size[i]);
            if(cdata_size[i] < 0) break;
            total_size += cdata_size[i];
        }
        
        chunk_num = i;
        //SUBREADprintf("\nchunk_num = %i\ttotal_size = %i\n", chunk_num, total_size);
        
        if(chunk_num) {
            // write this chunk to thread buffer
            while(1){
                int is_finished = 0;
                
                fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;
                
                pthread_mutex_lock(&thread_context->input_buffer_lock);
                
                // the number of bytes can be utilised given the two_chunk_len.
                int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder;
                int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr;
                
                //SUBREADprintf("T=%i\te=%i\tt=%i\n", current_thread_id, empty_bytes, tail_bytes);
                
                if(thread_context->input_buffer_remainder > global_context->input_buffer_max_size)
                {
                    SUBREADprintf("RMD=%d\n", thread_context->input_buffer_remainder );
                    assert(0);
                }
                
                if(tail_bytes < total_size + chunk_num * 4)
                    empty_bytes -= tail_bytes;
                
                // copy the new buffer to thread buffer.
                // format: read_number=n, read_chunk1, read_chunk2, ..., read_chunk_n
                if(empty_bytes >= total_size + chunk_num * 4)
                {
                    
                    if(tail_bytes < total_size + chunk_num * 4)
                    {
                        if(tail_bytes>=4)
                            memset(thread_context->input_buffer + thread_context->input_buffer_write_ptr, 0, 4);
                        thread_context->input_buffer_write_ptr = 0;
                        thread_context->input_buffer_remainder += tail_bytes;
                    }
                    int i = 0;
                    for(i = 0; i < chunk_num; i++) {
                        memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, & cdata_size[i], 4);
                        memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr + 4, chunk_in_buff[i] , cdata_size[i]);
                        thread_context->input_buffer_write_ptr += 4 + cdata_size[i];
                        thread_context->input_buffer_remainder += 4 + cdata_size[i];
                    }
                    
                    is_finished = 1;
                }
                
                pthread_mutex_unlock(&thread_context->input_buffer_lock);
                current_thread_id++;
                if(current_thread_id >= global_context->thread_number) current_thread_id = 0;
                
                if(is_finished) break;
                else usleep(tick_time);
            }
        } else {
            break; // Done or nothing to do
        }
    }
    return;
}

void write_decompressed_chunk(fc_thread_global_context_t * global_context, FILE * fp_in, char ** chunk_in_buff, int header_remain, char * binary_in_buff){
    
    int binary_remainder = header_remain;
    int binary_read_ptr = 0;
    
    while(1) {
        unsigned int real_len = 0;
        // most of the data must have been given out before this step.
        int current_thread_id = 0;
        int cdata_size[BLOCK_CHUNK_NUM] = {0};
        int i = 0;
        int no_of_reads = 0;
        
        // Get BLOCK_CHUNK_NUM blocks and decompress when there's enough space in buffer
        if(binary_remainder < 10000) {
            for (i = 0; i < BLOCK_CHUNK_NUM; i++) {
                cdata_size[i] = PBam_get_next_zchunk(fp_in, chunk_in_buff[i], 65537, & real_len);
                //SUBREADprintf("cs=%i; ", cdata_size[i]);
                if(cdata_size[i] < 0) break;
                int x1;
                
                // Shift remainder data to the front of the buffer
                if(binary_read_ptr>0)
                {
                    for(x1=0; x1< binary_remainder; x1++)
                        binary_in_buff[x1] = binary_in_buff [x1 + binary_read_ptr];
                    binary_read_ptr = 0;
                }
                
                if(cdata_size[i]>0)
                {
                    //printf("rmd = %i\n ", binary_remainder);
                    int new_binary_bytes = SamBam_unzip(binary_in_buff + binary_remainder , chunk_in_buff[i] , cdata_size[i]);
                    if(new_binary_bytes>=0)
                        binary_remainder += new_binary_bytes;
                    else	SUBREADprintf("ERROR: BAM GZIP FORMAT ERROR.\n");
                }
            }
        }
        
        // Get the whole buffer size.
        
        while(binary_remainder>4)
        {
            unsigned int binary_read_len = 0;
            memcpy(& binary_read_len , binary_in_buff + binary_read_ptr , 4);
            //printf("RLEN=%d; PTR=%d; RMD=%d\n", binary_read_len , binary_read_ptr, binary_remainder);
            if(binary_read_len > 10000)
            {
                binary_remainder = -1;
                SUBREADprintf("   A format error was detected in this BAM file.\n");
                SUBREADprintf("   The remaining part in the file is skipped.\n");
                SUBREADprintf("   Please check the file format using samtools.\n");
                break;
            }
            
            if(binary_read_len + 4<= binary_remainder)
            {
                no_of_reads ++;
                binary_read_ptr  += 4 + binary_read_len;
                binary_remainder -= 4 + binary_read_len;
            }
            else break;
        }
        
        if(binary_remainder <0)break;
        
        if(no_of_reads>0)
        {
            while(1)
            {
                int is_finished = 0;
                
                fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;
                
                pthread_mutex_lock(&thread_context->input_buffer_lock);
                
                // the number of bytes can be utilised given the two_chunk_len.
                int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder;
                int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr;
                
                if(thread_context->input_buffer_remainder > global_context->input_buffer_max_size)
                {
                    SUBREADprintf("RMD=%d\n", thread_context->input_buffer_remainder );
                    assert(0);
                }
                
                if(tail_bytes < binary_read_ptr + 4)
                    empty_bytes -= tail_bytes;
                
                // copy the new buffer to thread buffer.
                // format: buffer_size, read_buffer
                if(empty_bytes >= binary_read_ptr + 4)
                {
                    
                    if(tail_bytes < binary_read_ptr + 4)
                    {
                        if(tail_bytes>=4)
                            memset(thread_context->input_buffer + thread_context->input_buffer_write_ptr, 0, 4);
                        thread_context->input_buffer_write_ptr = 0;
                        thread_context->input_buffer_remainder += tail_bytes;
                    }
                    
                    memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, & binary_read_ptr, 4);
                    memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr + 4, binary_in_buff , binary_read_ptr);
                    thread_context->input_buffer_write_ptr += 4 + binary_read_ptr;
                    thread_context->input_buffer_remainder += 4 + binary_read_ptr;
                    is_finished = 1;
                }
                
                pthread_mutex_unlock(&thread_context->input_buffer_lock);
                current_thread_id++;
                if(current_thread_id >= global_context->thread_number) current_thread_id = 0;
                
                if(is_finished)
                    break;
                else usleep(tick_time);
                
            }
            
        }
        else
            break;
    }
    return;
}


int readSummary_single_file(fc_thread_global_context_t * global_context, unsigned int ** column_numbers, int * nexons,  int ** geneid, long ** start, long ** stop, unsigned char ** sorted_strand, long ** block_end_index, long ** block_min_start , long ** block_max_end);

int readSummary(int argc,char *argv[]){

	/*
	   This function counts the number of reads falling into each exon region.
	   The order of exons in the output is the same as that of exons included in the annotation.
	   The annotation, if provided as a file, should be sorted by chromosome name.
    */

	int isStrandChecked, isCVersion, isChimericDisallowed, isPEDistChecked, minMappingQualityScore=0, isInputFileResortNeeded, feature_block_size = 20, reduce_5_3_ends_to_one;

	char *nameFeatureTypeColumn, *nameGeneIDColumn,*debug_command;
    char * nameFeatureTypeList[MAX_FEATURE_TYPE_NUM + 1];
    int featureTypeNum = 0;
    int is_independent_assign, multithreadunzip;
    
    int run_mode; // The running mode of VERSE
    int f_length; // output length info?

	long curchr, curpos;
	char * curchr_name;
	
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadDetailReport, isBothEndRequired, isMultiMappingAllowed, fiveEndExtension, threeEndExtension, minReadLength, minReadOverlap, maxReadNonoverlap, isSplitAlignmentOnly, is_duplicate_ignored, minDifAmbiguous, nonempty_modified;

	int isSAM, isGTF, n_input_files=0;
	char * cmd_rebuilt = NULL;

	isCVersion = ((argv[0][0]=='C')?1:0);
    run_mode = atoi(argv[4]);
	isPE = atoi(argv[5]);
	minPEDistance = atoi(argv[6]);
	maxPEDistance = atoi(argv[7]);

	isSAM = atoi(argv[8]);
	unsigned short thread_number;
	if(argc > 10)
		thread_number = atoi(argv[10]);
	else	thread_number = 4;
	if(argc > 11)
		isGTF = atoi(argv[11]);
	else	isGTF = 0;
	if(argc > 12)
		isStrandChecked = atoi(argv[12]);
	else	isStrandChecked = 0;
	if(argc > 13)
		isReadDetailReport = atoi(argv[13]);
	else	isReadDetailReport = 0;
	if(argc > 14)
		isBothEndRequired = atoi(argv[14]);
	else	isBothEndRequired = 0;
	if(argc > 15)
		isChimericDisallowed = atoi(argv[15]);
	else	isChimericDisallowed = 0;
	if(argc > 16)
		isPEDistChecked = atoi(argv[16]);
	else	isPEDistChecked = 0;
	if(argc > 17)
		nameFeatureTypeColumn = argv[17];
	else	nameFeatureTypeColumn = "exon";
	if(argc > 18)
		nameGeneIDColumn = argv[18];
	else	nameGeneIDColumn = "gene_id";
	if(argc > 19)
		minMappingQualityScore = atoi(argv[19]);
	else	minMappingQualityScore = 0;
	if(argc > 20)
		isMultiMappingAllowed = atoi(argv[20]);
	else	isMultiMappingAllowed = 1;
    if(argc > 21)
        is_independent_assign = atoi(argv[21]);
    else	is_independent_assign = 0;
	if(argc > 22)
	{
		cmd_rebuilt = argv[22];
		if(cmd_rebuilt == NULL || cmd_rebuilt[0]==' '||cmd_rebuilt[0]==0)
			cmd_rebuilt=NULL;
	}
	else	cmd_rebuilt = NULL;
	if(argc>23)
		isInputFileResortNeeded = atoi(argv[23]);
	else	isInputFileResortNeeded = 0;
	if(thread_number<1) thread_number=1;

	if(argc>25)
		fiveEndExtension = atoi(argv[25]);
	else 	fiveEndExtension = 0;

	if(argc>26)
		threeEndExtension = atoi(argv[26]);
	else	threeEndExtension = 0;

    if(argc>27)
        minReadLength = atoi(argv[27]);
    else	minReadLength = 1;
    
	if(argc>28)
		minReadOverlap = atoi(argv[28]);
	else	minReadOverlap = 1;
	
    if(argc>29)
        maxReadNonoverlap = atoi(argv[29]);
    else    maxReadNonoverlap = 0x7fff;
        
	if(argc>30)
		isSplitAlignmentOnly = atoi(argv[30]);
	else	isSplitAlignmentOnly = 0;

	if(argc>31)
		reduce_5_3_ends_to_one = atoi(argv[31]);	// 0 : no reduce; 1: reduce to 5' end; 2: reduce to 3' end.
	else	reduce_5_3_ends_to_one = 0;

	if(argc>32 && strlen(argv[32])>0 && argv[32][0]!=' ')
		debug_command = argv[32];
	else
		debug_command = " ";

	if(argc>33)
		is_duplicate_ignored = atoi(argv[33]);
	else
		is_duplicate_ignored = 0;
    
    if(argc>34)
        minDifAmbiguous = atoi(argv[34]);
    else	minDifAmbiguous = 0;

    if(argc>35)
        f_length = atoi(argv[35]);
    else    f_length = 0;
    
    if(argc>36)
        nonempty_modified = atoi(argv[36]);
    else    nonempty_modified = 0;

    if(argc>37)
        multithreadunzip = atoi(argv[37]);
    else    multithreadunzip = 0;
    

	unsigned int buffer_size = 1024*1024*64;

    int i = 0;
    char * start_ptr = nameFeatureTypeColumn;
    while(1)
    {
        char * end_ptr = strchr(start_ptr,';');
        if(end_ptr == NULL) end_ptr = nameFeatureTypeColumn + strlen(nameFeatureTypeColumn);
        nameFeatureTypeList[i] = malloc(end_ptr - start_ptr + 1);
        sprintf(nameFeatureTypeList[i], "%.*s", (int)(end_ptr - start_ptr), start_ptr);
        i++;
        if(i > MAX_FEATURE_TYPE_NUM)
        {
            SUBREADprintf("Too many feature types specified. You specified %i, current maximum is %i.\n", i, MAX_FEATURE_TYPE_NUM);
            SUBREADprintf("Try increasing the MAX_FEATURE_TYPE_NUM in the program.\n");
            return -1;
        }
        start_ptr = end_ptr + 1;
        if(start_ptr - nameFeatureTypeColumn >= strlen(nameFeatureTypeColumn)) break;
    }
    featureTypeNum = i;

	fc_thread_global_context_t global_context;
	fc_thread_init_global_context(& global_context, buffer_size, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance, isStrandChecked, (char *)argv[3] , isReadDetailReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn, minMappingQualityScore,isMultiMappingAllowed, isSAM, cmd_rebuilt, isInputFileResortNeeded, feature_block_size, isCVersion, fiveEndExtension, threeEndExtension, minReadLength, minReadOverlap, maxReadNonoverlap, isSplitAlignmentOnly, reduce_5_3_ends_to_one, debug_command, is_duplicate_ignored, run_mode, minDifAmbiguous, featureTypeNum, nameFeatureTypeList, is_independent_assign, nonempty_modified, multithreadunzip);
	print_FC_configuration(&global_context, argv[1], argv[2], argv[3], global_context.is_SAM_file, & n_input_files, isReadDetailReport, nameFeatureTypeColumn);

	// Loading the annotations.
	// Nothing is done if the annotation does not exist.

    int nfeatures[featureTypeNum];
    int nexons[featureTypeNum];
    fc_feature_info_t * loaded_features[featureTypeNum];
    fc_feature_info_t ** merged_features[featureTypeNum];
    
    long * start[featureTypeNum], * stop[featureTypeNum];
    int * geneid[featureTypeNum];
    unsigned char * sorted_strand[featureTypeNum];
    long * block_min_start[featureTypeNum], * block_max_end[featureTypeNum], * block_end_index[featureTypeNum];
    
	print_in_box(84,0,0,"Load annotation file %s %c[0m...", argv[1], CHAR_ESC);
    
    int f_idx;
    for(f_idx = 0; f_idx < featureTypeNum; f_idx++)
    {
        nfeatures[f_idx] = load_feature_info(&global_context,argv[1], &loaded_features[f_idx], f_idx);
        if(nfeatures[f_idx]<1){
            SUBREADprintf("Failed to open the annotation file %s, maybe its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeList[f_idx]);
            return -1;
        }
        nexons[f_idx] = sort_feature_info(&global_context, nfeatures[f_idx], loaded_features[f_idx], &geneid[f_idx], &start[f_idx], &stop[f_idx], &sorted_strand[f_idx], & block_end_index[f_idx], & block_min_start[f_idx], & block_max_end[f_idx], & merged_features[f_idx], f_idx);
        if(nexons[f_idx]<1){
            SUBREADprintf("Failed to merge feature type %s in file %s, maybe its format is incorrect.\n",nameFeatureTypeList[f_idx], argv[1]);
            return -1;
        }

        print_in_box(80,0,0,"   %ss after merge : %i", nameFeatureTypeList[f_idx], nexons[f_idx]);
        print_in_box(80,0,0,"   Genes : %d", global_context.gene_name_table[f_idx] -> numOfElements);
        print_in_box(80,0,0,"   Chromosomes : %d", global_context.exontable_nchrs[f_idx]);
        
        print_in_box(80,0,0,"");
        global_context.original_exons[f_idx] = nfeatures[f_idx];
        global_context.exontable_exons[f_idx] = nexons[f_idx];
    }

	char * tmp_pntr = NULL;
	char * files_used = malloc(strlen(argv[2])+1);
	strcpy(files_used, argv[2]);
	char * next_fn = strtok_r(files_used,";", &tmp_pntr);
    if(next_fn == NULL || strlen(next_fn) < 1)
    {
        SUBREADprintf("Failed to open the input file : %s\n", argv[2]);
        return -1;
    }
    
    unsigned int * column_numbers[featureTypeNum];
    for(f_idx = 0; f_idx < featureTypeNum; f_idx++)
    {
        column_numbers[f_idx] = calloc(nexons[f_idx], sizeof(unsigned int));
    }
    
    strcpy(global_context.input_file_name, next_fn);
    strcpy(global_context.raw_input_file_name, next_fn);
    
    int ret_int = readSummary_single_file(& global_context, column_numbers, nexons, geneid, start, stop, sorted_strand, block_end_index, block_min_start, block_max_end);
    
    if(ret_int != 0)
    {
        for(f_idx = 0; f_idx < featureTypeNum; f_idx++)
        {
            free(column_numbers[f_idx]);
        }
        return -1;
    }
    free(files_used);
    for(f_idx = 0; f_idx < featureTypeNum; f_idx ++)
    {
        char tmp_fname[330];
        sprintf(tmp_fname, "%s.%s.txt", argv[3], nameFeatureTypeList[f_idx]);
        fc_write_final_gene_results(&global_context, geneid[f_idx], start[f_idx], stop[f_idx], sorted_strand[f_idx], tmp_fname, nexons[f_idx], column_numbers[f_idx], argv[2], f_length, f_idx);
    }
    fc_write_final_counts(&global_context, argv[3], argv[2], argv[1]);

	print_FC_results(&global_context);
    
	KeyValuePair * cursor;
	int bucket;
    
    for(f_idx = 0; f_idx < featureTypeNum; f_idx++)
    {
        for(bucket=0; bucket < global_context.exontable_chro_table[f_idx] -> numOfBuckets; bucket++)
        {
            cursor = global_context.exontable_chro_table[f_idx] -> bucketArray[bucket];
            while (1)
            {
                if (!cursor) break;
                fc_chromosome_index_info * del_chro_info = cursor->value;
                free(del_chro_info->reverse_table_start_index);
                //free(del_chro_info->reverse_table_end_index);
                free((void *)cursor -> key);
                free(del_chro_info);
                cursor = cursor->next;
            }
        }
        free(column_numbers[f_idx]);
        HashTableDestroy(global_context.gene_name_table[f_idx]);
        free(global_context.gene_name_array[f_idx]);
        HashTableDestroy(global_context.exontable_chro_table[f_idx]);
        free(loaded_features[f_idx]);
        free(merged_features[f_idx]);
        free(geneid[f_idx]);
        free(start[f_idx]);
        free(sorted_strand[f_idx]);
        free(block_min_start[f_idx]);
        free(block_max_end[f_idx]);
        free(block_end_index[f_idx]);
        free(stop[f_idx]);
        free(global_context.unistr_buffer_space[f_idx]);
        free(nameFeatureTypeList[f_idx]);
    }
    
    free(global_context.longest_chro_name);
    free(global_context.exontable_nchrs);
    free(global_context.gene_name_array);
    free(global_context.gene_name_table);
    free(global_context.exontable_chro_table);
    free(global_context.unistr_buffer_space);
    free(global_context.unistr_buffer_used);
    free(global_context.unistr_buffer_size);
    free(global_context.original_exons);
    free(global_context.exontable_exons);
    free(global_context.read_counters);
    
	if(global_context.detail_output_fp) fclose(global_context. detail_output_fp);
    
    if(global_context.unprocessed_cnts)
    {
        int i = 0;
        for(i = 0; i < global_context.unprocessed_cnts; i++)
        {
            free(global_context.unprocessed_reads[i]);
        }
        free(global_context.unprocessed_reads);
    }
	return 0;
}

int readSummary_single_file(fc_thread_global_context_t * global_context, unsigned int ** column_numbers, int * nexons, int ** geneid, long ** start, long ** stop, unsigned char ** sorted_strand, long ** block_end_index, long ** block_min_start , long ** block_max_end)
{
	FILE *fp_in = NULL;
	int read_length = 0;
	int is_first_read_PE=0;
	//char * line = (char*)calloc(MAX_LINE_LENGTH, 1);
	char * file_str = "";

	if(strcmp(global_context->input_file_name,"STDIN")!=0)
	{
		int file_probe = is_certainly_bam_file(global_context->input_file_name, &is_first_read_PE);
        int alert_user = 0;
        
        // a Singel-end SAM/BAM file cannot be assigned as a PE SAM/BAM file;
        // but a PE SAM/BAM file may be assigned as a SE file if the user wishes to do so.
        if(is_first_read_PE==0)
        {
            global_context -> is_paired_end_data = 0;
            alert_user = 1;
        }
        

		if(file_probe == 1){
			global_context->is_SAM_file = 0;
		}
		else if(file_probe == 0) global_context->is_SAM_file = 1;

		global_context -> start_time = miltime();

		file_str = "SAM";
		if(file_probe == 1) file_str = "BAM" ;
		if(file_probe == -1)
        {
            file_str = "Unknown";
            SUBREADprintf("ERROR: The file does not exist, or its format is not supported by VERSE or it contains no read.\n");
            SUBREADprintf("No counts were generated for this file.\n");
            return -1;
        }

        print_in_box(80,0,0,"Process %s file %s...", file_str, global_context->input_file_name);
        if(is_first_read_PE)
        {
            if(global_context->is_paired_end_data)
                print_in_box(80,0,0,"   Quantifying in paired-end mode.");
            else {
                print_in_box(80,0,0,"   \x1B[31mALERT\x1b[0m : You specified single-end, but the file contains paired-end data.");
                print_in_box(80,0,0,"   Quantifying in single-end mode.");
            }
            
        }
        else {
            if(alert_user)
                print_in_box(80,0,0,"   \x1B[31mALERT\x1b[0m : You specified paired-end, but the data seem to be single-end. ");
            print_in_box(80,0,0,"   Quantifying in single-end mode.");
        }
	}

	int isInputFileResortNeeded = global_context->is_input_file_resort_needed;

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		FILE * exist_fp = fopen(global_context->input_file_name,"r");
		if(!exist_fp)
		{
			SUBREADprintf("ERROR: Failed to open file %s\n",  global_context->input_file_name);
			SUBREADprintf("No counts were generated for this file.\n");
			return -1;
		}
		fclose(exist_fp);
	}

	if(strcmp(global_context->input_file_name,"STDIN")!=0)
		if(warning_file_type(global_context->input_file_name, global_context->is_SAM_file?FILE_TYPE_SAM:FILE_TYPE_BAM))
			global_context->is_unpaired_warning_shown=1;
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		if(resort_input_file(global_context)) return -1;
	int isSAM = global_context->is_SAM_file;
	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.

	if(strcmp("STDIN",global_context->input_file_name)==0)
		fp_in = stdin;
	else
		fp_in = fopen(global_context->input_file_name,"r");

	// begin to load-in the data.
    if(global_context->is_paired_end_data)
        print_in_box(80,0,0,"   Assign read pairs to features...");
    else
        print_in_box(80,0,0,"   Assign reads to features...");
    
    if(global_context -> is_read_details_out)
    {
        char tmp_fname[350];
        sprintf(tmp_fname, "%s.detail.txt", global_context -> output_file_name);
        global_context -> detail_output_fp = fopen(tmp_fname, "w");
        if(!global_context -> detail_output_fp)
        {
            SUBREADprintf("ERROR: Unable to create file '%s'; the read assignment details are not written.\n", tmp_fname);
        }
    }
    else
        global_context -> detail_output_fp = NULL;
    

	int isPE = global_context->is_paired_end_data;

	SamBam_Reference_Info * sb_header_tab = NULL;
	
    char * chunk_in_buff[BLOCK_CHUNK_NUM];
    int i = 0;
    for(i = 0; i < BLOCK_CHUNK_NUM; i++) {
        chunk_in_buff[i] = malloc(65537);
    }

	char * binary_in_buff = malloc(65537 * BLOCK_CHUNK_NUM + 10000);
    int is_last_read = 0;

    char * ret = NULL;
    
    int remainder_read_data_len = 0;
    
    if(!isSAM) {
        // Load header
        
        PBum_load_header(fp_in, &sb_header_tab, binary_in_buff,  & remainder_read_data_len);
        //printf("RMD=%d\n", remainder_read_data_len);
        if(remainder_read_data_len && global_context->multithread_unzipping) {
            SUBREADprintf("Your BAM file is outdated and does not support multithreaded decompression.\n");
            SUBREADprintf("You can rerun without multithreaded decompression,\n");
            SUBREADprintf("or update your BAM using samtools [samtools view -b].\n");
            return -1;
        }
        
        global_context->sambam_chro_table = sb_header_tab;
    }
    
    fc_thread_start_threads(global_context, nexons, geneid, start, stop, sorted_strand, block_end_index, block_min_start , block_max_end, read_length);
    fc_thread_thread_context_t * one_thread_context = global_context->thread_contexts;
    
    if(isSAM && global_context->thread_number==1)
    {
        while (1){
            int is_second_read;
            
            for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
            {
                if(is_last_read)
                {
                    ret = NULL;
                    break;
                }
                if(one_thread_context->step_back)   // previous 2 reads sent to process is unpaired, the second read is now the first read
                {
                    is_second_read = 1;
                    strcpy(one_thread_context -> line_buffer1, one_thread_context -> line_buffer2);
                    one_thread_context -> step_back = 0;
                }
                char * lbuf = is_second_read?one_thread_context -> line_buffer2:one_thread_context -> line_buffer1;
                while(1)
                {
                    ret = fgets(lbuf, MAX_LINE_LENGTH, fp_in);
                    if(!ret) break;
                    if(lbuf[0] == '@')
                    {
                        int retlen = strlen(ret);
                        if(ret[retlen-1]!='\n')
                        {
                            while(1){
                                int nch = getc(fp_in);
                                if(nch == EOF || nch == '\n') break;
                            }
                        }
                    }
                    else break;
                }
                
                if(!ret)
                {
                    if(is_second_read) // paired-end, last read unpaired
                    {
                        char * fake_string = "THIS_IS_A_FAKE_READ_BY_VERSE\t89\t23\t35249027\t255\t0M100D\t*\t0\t0\tN\t#\tNH:i:1\tHI:i:1\n";
                        is_last_read = 1;
                        sprintf(one_thread_context -> line_buffer2, "%s", fake_string);
                    }
                    break;
                }
                
                if(read_length < 1)
                {
                    int tab_no = 0;
                    int read_len_tmp=0, read_cursor;
                    int curr_line_len = strlen(lbuf);
                    for(read_cursor=0; read_cursor<curr_line_len; read_cursor++)
                    {
                        if(lbuf[read_cursor] == '\t')
                            tab_no++;
                        else
                        {
                            if(tab_no == 9)	// SEQ
                                read_len_tmp++;
                        }
                    }
                    read_length = read_len_tmp;
                    global_context->read_length = read_length;
                }
            }
            
            if(!ret && !is_second_read) break;
            if(isPE) probe_pair_name(one_thread_context);
            process_line_buffer(global_context, one_thread_context);
        }
    }
    
    /**************************** BAM file handler ****************************/
    else if(!isSAM)
    {
        // Process the bulk data
        
        if(global_context -> multithread_unzipping) {
            write_compressed_block(global_context, fp_in, chunk_in_buff);
        } else {
            write_decompressed_chunk(global_context, fp_in, chunk_in_buff, remainder_read_data_len, binary_in_buff);
        }
    }
    else
    {
        SUBREADprintf("ERROR: Multithread SAM or other formats is currently not supported!!\n");
        SUBREADprintf("Please retry using one thread or convert the SAM to BAM using samtools.\n");
        return -1;
    }

    
    int j = 0;
    for(j = 0 ; j < BLOCK_CHUNK_NUM; j++) {
        free(chunk_in_buff[j]);
    }

	free(binary_in_buff);
	//free(preload_line);
	global_context->is_all_finished = 1;            // bulk finish
    if(isSAM || !isPE)
        global_context->is_really_finished = 1;     // if single-end or SAM, complete finish

	if(global_context->thread_number > 1 || !isSAM)
		fc_thread_wait_threads(global_context);

    // pool the unprocessed reads
    if(!isSAM && isPE)     // if paired-end BAM
    {
        combine_sort_unprocessed(global_context);
        //printf("Unprocessed cnts: %lli\n", global_context->unprocessed_cnts);
        if(global_context->unprocessed_cnts == 0)
            global_context->is_really_finished = 1;     // if there's no unprocessed
    }
    
    fc_thread_merge_results(global_context, column_numbers);

	fc_thread_destroy_thread_context(global_context);
    
    // restart a single thread to deal with the unassigned
    if(!isSAM && isPE && global_context -> unprocessed_cnts)
    {
        fc_thread_start_threads(global_context, nexons, geneid, start, stop, sorted_strand, block_end_index, block_min_start , block_max_end, read_length);
        
        fc_thread_thread_context_t * one_thread_context = global_context->thread_contexts;
        int is_second_read, read_idx = 0;
        
        while(read_idx < global_context -> unprocessed_cnts)
        {
            for(is_second_read=0; is_second_read < 2; is_second_read++)
            {
                char * lbuf = is_second_read?one_thread_context -> line_buffer2:one_thread_context -> line_buffer1;
                sprintf(lbuf, "%s\n", global_context -> unprocessed_reads[read_idx]);
                
                read_idx++;
                if (is_second_read == 0 && read_idx >= global_context -> unprocessed_cnts) // last read unpaired
                {
                    char * fake_string = "THIS_IS_A_FAKE_READ_BY_VERSE\t89\t23\t35249027\t255\t0M100D\t*\t0\t0\tN\t#\tNH:i:1\tHI:i:1\n";
                    global_context->missing_mates++;
                    sprintf(one_thread_context -> line_buffer2, "%s", fake_string);
                    break;
                }
            }
            
            probe_pair_name(one_thread_context);
            process_line_buffer(global_context, one_thread_context);
            
            if(one_thread_context -> step_back)
            {
                if(is_second_read)         //not last read unpaired
                {
                    read_idx--;
                    global_context->missing_mates++;
                    one_thread_context -> step_back = 0;
                }
            }
        }
        global_context -> is_really_finished = 1;
        fc_thread_merge_results(global_context, column_numbers);
        fc_thread_destroy_thread_context(global_context);
    }

	if(strcmp("STDIN",global_context->input_file_name)!=0)
		fclose(fp_in);
	if(sb_header_tab) free(sb_header_tab);
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		unlink(global_context->input_file_name);
	return 0;
}

int main(int argc, char ** argv)
{
    char * Rargv[38];
	char annot_name[300];
	char * out_name = malloc(300);
	int cmd_rebuilt_size = 200;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);
	char nameFeatureTypeColumn[200];
	char nameGeneIDColumn[66];
	int min_qual_score = 0;
	int min_dist = 50;
	int max_dist = 600;
	char debug_command[10];
	char min_dist_str[11];
	char max_dist_str[11];
	char min_qual_score_str[11];
    char run_mode_str[11];
	char feature_block_size_str[11];
	char Strand_Sensitive_Str[11];
	char * very_long_file_names;
	int is_Input_Need_Reorder = 0;
	int is_PE = 1;
	int is_SAM = 1;
	int is_Both_End_Mapped = 0;
	int feature_block_size = 14;
	int Strand_Sensitive_Mode = 0;
	int is_ReadDetail_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int is_Multi_Mapping_Allowed = 0;
	int is_Split_Alignment_Only = 0;
	int is_duplicate_ignored = 0;
	int reduce_5_3_ends_to_one = 0;
	int threads = 1;
	int isGTF = 1;
	char nthread_str[4];
	int option_index = 0;
	int c;
	int very_long_file_names_size = 200;
	int fiveEndExtension = 0, threeEndExtension = 0, minReadOverlap = 1, minReadLength = 1, maxReadNonoverlap = 0x7fff, minDifAmbiguous = 0;
	char strFiveEndExtension[11], strThreeEndExtension[11], strMinReadLength[11], strMinReadOverlap[11], strMaxReadNonoverlap[11], strMinDifAmbiguous[11];
    
    int run_mode = 0; // default mode false --> featureCounts
    int f_length = 0; // default not output feature length.
    int assignIndependently = 0; // assign hierarchically if 0, independently if 1. Default is hierarchical.
    int nonemptyModified = 0;
    int multithreadunzip = 0;
    
	very_long_file_names = malloc(very_long_file_names_size);
	very_long_file_names [0] = 0;

	debug_command[0] = 0;

	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
    annot_name[0]=0;out_name[0]=0;
	

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		if(strlen(cmd_rebuilt) + 300 > cmd_rebuilt_size)
		{
			cmd_rebuilt_size*=2;
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		sprintf(cmd_rebuilt+strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "g:t:T:o:a:d:D:L:Q:z:bs:SCBPRZlv?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'v':
				core_version_number("VERSE");
				return 0;
            case 'Z':
                print_run_mode();
                return 0;
			case 'Q':
				min_qual_score = atoi(optarg);
				break;
            case 'z':
                run_mode = atoi(optarg);
                break;
            case 'l':
                f_length = 1;
                break;
			case 't':
				strcpy(nameFeatureTypeColumn, optarg);
				break;
			case 'g':
				while((*optarg) == ' ') optarg++;
				strcpy(nameGeneIDColumn, optarg);
				break;
			case 'T':
				threads = atoi(optarg);
				break;
			case 'd':
				min_dist = atoi(optarg);
				break;
			case 'D':
				max_dist = atoi(optarg);
				break;
			case 'C':
				is_Chimeric_Disallowed = 1;
				break;
            case 'S':
                is_Input_Need_Reorder = 1;
                break;
			case 'P':
				is_PE_Dist_Checked = 1;
				break;
			case 'B':
				is_Both_End_Mapped = 1;
				break;
			case 'R':
				is_ReadDetail_Report = 1;
				break;
			case 's':
				Strand_Sensitive_Mode = atoi(optarg);
				break;
			case 'o':
				term_strncpy(out_name, optarg,299);
				break;
			case 'a':
				term_strncpy(annot_name, optarg,299);
				break;
			case 'L':
				feature_block_size = atoi(optarg);
				break;
			case 0 :	// long options
                if(strcmp("singleEnd", long_options[option_index].name)==0)
                {
                    is_PE = 0;
                }
                
				if(strcmp("readExtension5", long_options[option_index].name)==0)
				{
					fiveEndExtension = atoi(optarg);
					fiveEndExtension = max(0, fiveEndExtension);
				}

				if(strcmp("readExtension3", long_options[option_index].name)==0)
				{
					threeEndExtension = atoi(optarg);
					threeEndExtension = max(0, threeEndExtension);
				}

				if(strcmp("minReadOverlap", long_options[option_index].name)==0)
				{
					minReadOverlap = atoi(optarg);
				}
                
                if(strcmp("maxReadNonoverlap", long_options[option_index].name)==0)
                {
                    maxReadNonoverlap = atoi(optarg);
                }

				if(strcmp("debugCommand", long_options[option_index].name)==0)
				{
					strcpy(debug_command, optarg);
				}

				if(strcmp("ignoreDup", long_options[option_index].name)==0)
				{
					is_duplicate_ignored = 1 ;
				}

				if(strcmp("read2pos", long_options[option_index].name)==0)
				{
					if(optarg[0]=='3')
						reduce_5_3_ends_to_one = REDUCE_TO_3_PRIME_END;
					else if(optarg[0]=='5')
						reduce_5_3_ends_to_one = REDUCE_TO_5_PRIME_END;
				}				

				if(strcmp("countSplitAlignmentsOnly", long_options[option_index].name)==0)
				{
					is_Split_Alignment_Only = 1;
				}
                
                if(strcmp("minDifAmbiguous", long_options[option_index].name)==0)
                {
                    minDifAmbiguous = atoi(optarg);
                }
                if(strcmp("assignIndependently", long_options[option_index].name)==0)
                {
                    assignIndependently = 1;
                }
                if(strcmp("nonemptyModified", long_options[option_index].name)==0)
                {
                    nonemptyModified = 1;
                }
                if(strcmp("multithreadDecompress", long_options[option_index].name)==0)
                {
                    multithreadunzip = 1;
                }

				break;
			case '?':
			default :
                SUBREADprintf("Outputting documentation, VERSE is not running\n");
				print_usage();
                free(very_long_file_names);
                free(out_name);
                free(cmd_rebuilt);
				return -1;
				break;
		}


    if(out_name[0] == 0 || annot_name[0] == 0 || argc == optind || minReadOverlap < 1 || minReadLength < 1 || maxReadNonoverlap < 0 || maxReadNonoverlap > 0x7fff || minDifAmbiguous < 0 || (minDifAmbiguous > 0 && run_mode != 5) || (run_mode != 3 && nonemptyModified))
	{
        print_usage();
        free(very_long_file_names);
        free(out_name);
        free(cmd_rebuilt);
        return -1;
	}
    
    if(run_mode < 0 || run_mode > 5)
    {
        print_run_mode();
        free(very_long_file_names);
        free(out_name);
        free(cmd_rebuilt);
        return -1;
    }
    
    /*
    if(threads > 1 && is_ReadDetail_Report)
    {
        SUBREADprintf("\nDetail file cannot be written when using multiple threads, please use one thread instead.\n\n");
        free(very_long_file_names);
        free(out_name);
        free(cmd_rebuilt);
        return -1;
    }
     */

	for(; optind < argc; optind++)
	{
		int curr_strlen = strlen(very_long_file_names);
		if( very_long_file_names_size - curr_strlen <300)
		{
			very_long_file_names_size *=2;
			//printf("CL=%d ; NS=%d\n", curr_strlen , very_long_file_names_size);
			very_long_file_names=realloc(very_long_file_names , very_long_file_names_size);
		}

		strcat(very_long_file_names, argv[optind]);
		strcat(very_long_file_names, ";");
	}

	very_long_file_names[strlen(very_long_file_names)-1]=0;

	sprintf(strFiveEndExtension, "%d", fiveEndExtension);
	sprintf(strThreeEndExtension, "%d", threeEndExtension);
    sprintf(strMinReadLength, "%d", minReadLength);
	sprintf(strMinReadOverlap, "%d", minReadOverlap);
    sprintf(strMaxReadNonoverlap, "%d", maxReadNonoverlap);
	sprintf(nthread_str,"%d", threads);
	sprintf(min_dist_str,"%d",min_dist);
	sprintf(max_dist_str,"%d",max_dist);
	sprintf(min_qual_score_str,"%d", min_qual_score);
	sprintf(feature_block_size_str,"%d", feature_block_size);
	sprintf(Strand_Sensitive_Str,"%d", Strand_Sensitive_Mode);
    sprintf(strMinDifAmbiguous, "%d", minDifAmbiguous);
    
    sprintf(run_mode_str, "%d", run_mode);
	Rargv[0] = "Cverse";
	Rargv[1] = annot_name;
	Rargv[2] = very_long_file_names;
	Rargv[3] = out_name;
    Rargv[4] = run_mode_str;
	Rargv[5] = is_PE?"1":"0";
	Rargv[6] = min_dist_str;
	Rargv[7] = max_dist_str;
	Rargv[8] = is_SAM?"1":"0";
	Rargv[10] = nthread_str;
	Rargv[11] = isGTF?"1":"0";
	Rargv[12] = Strand_Sensitive_Str;
	Rargv[13] = is_ReadDetail_Report?"1":"0";
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	Rargv[19] = min_qual_score_str;
	Rargv[20] = is_Multi_Mapping_Allowed == ALLOW_PRIMARY_MAPPING?"2":(is_Multi_Mapping_Allowed == ALLOW_ALL_MULTI_MAPPING?"1":"0");
    Rargv[21] = assignIndependently?"1":"0";
	Rargv[22] = cmd_rebuilt;
	Rargv[23] = is_Input_Need_Reorder?"1":"0";
	Rargv[24] = feature_block_size_str;
	Rargv[25] = strFiveEndExtension;
	Rargv[26] = strThreeEndExtension;
    Rargv[27] = strMinReadLength;
	Rargv[28] = strMinReadOverlap;
    Rargv[29] = strMaxReadNonoverlap;
	Rargv[30] = is_Split_Alignment_Only?"1":"0";
	Rargv[31] = (reduce_5_3_ends_to_one == 0?"0":(reduce_5_3_ends_to_one==REDUCE_TO_3_PRIME_END?"3":"5"));
	Rargv[32] = debug_command;
	Rargv[33] = is_duplicate_ignored?"1":"0";
    Rargv[34] = strMinDifAmbiguous;
    Rargv[35] = f_length?"1":"0";
    Rargv[36] = nonemptyModified?"1":"0";
    Rargv[37] = multithreadunzip?"1":"0";

	readSummary(38, Rargv);

	free(very_long_file_names);
	free(out_name);
	free(cmd_rebuilt);

	return 0;

}


