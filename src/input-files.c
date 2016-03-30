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
  

#include <signal.h>
#include <dirent.h>

#include <unistd.h>
#include <assert.h>

#include "input-files.h"
#include "sambam-file.h"



char * _SAMSORT_SNP_delete_temp_prefix = NULL;
void SAM_SORT_SIGINT_hook(int param)
{
    int xk1, last_slash = -1;
    if(_SAMSORT_SNP_delete_temp_prefix != NULL)
    {
        char del2[300], del_suffix[200], del_name[400];
        SUBREADprintf("\n\nReceived a terminal signal. The temporary files were removed.\n");
        for(xk1=0; _SAMSORT_SNP_delete_temp_prefix[xk1]; xk1++)
        {
            if(_SAMSORT_SNP_delete_temp_prefix[xk1]=='/') last_slash = xk1;
            else if(_SAMSORT_SNP_delete_temp_prefix[xk1]=='\\')
            {
                SUBREADprintf("The file name is unknown.\n");
                return;
            }
        }
        if(last_slash>=0)
        {
            memcpy(del2, _SAMSORT_SNP_delete_temp_prefix, last_slash);
            del2[last_slash] = 0;
            strcpy(del_suffix , _SAMSORT_SNP_delete_temp_prefix + last_slash + 1);
        }
        else
        {
            strcpy(del2,".");
            strcpy(del_suffix , _SAMSORT_SNP_delete_temp_prefix);
        }
        
        if(strlen(del_suffix)>8)
        {
            DIR           *d;
            struct dirent *dir;
            
            d = opendir(del2);
            if (d)
            {
                while ((dir = readdir(d)) != NULL)
                {
                    if(strstr(dir->d_name, del_suffix))
                    {
                        //printf("%s\n", dir->d_name);
                        strcpy(del_name, del2);
                        strcat(del_name, "/");
                        strcat(del_name, dir->d_name);
                        unlink(del_name);
                    }
                }
            }
        }
        
    }
    
    exit(param);
}

void * old_sig_TERM = NULL, * old_sig_INT = NULL;
int sort_SAM_create(SAM_sort_writer * writer, char * output_file, char * tmp_path)
{
    char tmp_fname[MAX_FILE_NAME_LENGTH+40];
    memset(writer, 0, sizeof(SAM_sort_writer));
    
    old_sig_TERM = signal (SIGTERM, SAM_SORT_SIGINT_hook);
    old_sig_INT = signal (SIGINT, SAM_SORT_SIGINT_hook);
    
    sprintf(writer -> tmp_path, "%s/temp-sort-%06u-%08X-", tmp_path, getpid(), rand());
    _SAMSORT_SNP_delete_temp_prefix = writer -> tmp_path;
    
    sprintf(tmp_fname, "%s%s", writer -> tmp_path, "headers.txt");
    writer -> all_chunks_header_fp = fopen(tmp_fname,"w");
    if(!writer -> all_chunks_header_fp) return -1;
    fclose(writer -> all_chunks_header_fp);
    unlink(tmp_fname);
    
    writer -> out_fp = fopen(output_file,"w");
    if(!writer -> out_fp) return -1;
    
    return 0;
}

char * fgets_noempty(char * buf, int maxlen, FILE * fp)
{
    char * ret;
    while(1)
    {
        ret = fgets(buf, maxlen, fp);
        if(!ret)return NULL;
        if(ret[0]!='\n') return ret;
    }
}

int probe_file_type(char * fname, int * is_first_read_PE)
{
    FILE * fp = fopen(fname, "rb");
    if(!fp) return FILE_TYPE_NONEXIST;
    
    int ret = FILE_TYPE_UNKNOWN;
    int nch;
    char *test_buf=malloc(5000);
    
    nch = fgetc(fp);
    
    if(feof(fp))
    ret = FILE_TYPE_EMPTY;
    
    else
    {
        if(nch == '@')	// FASTQ OR SAM
        {
            char * rptr = fgets_noempty(test_buf, 4999, fp);
            int second_line_len = 0;
            if(rptr)
            {
                rptr = fgets_noempty(test_buf, 4999, fp);
                if(rptr)
                {
                    second_line_len = strlen(test_buf);
                    int tabs = 0, x1;
                    for(x1=0;x1<4999;x1++)
                    {
                        if(test_buf[x1]=='\n' || !test_buf[x1]) break;
                        if(test_buf[x1]=='\t'){
                            tabs++;
                            continue;
                        }
                        
                        if(tabs == 1)
                        if(!isdigit(test_buf[x1]))break;
                    }
                    if(rptr[0]=='@' || tabs>7)
                    ret = FILE_TYPE_SAM;
                }
            }
            if(ret == FILE_TYPE_UNKNOWN)
            {
                rptr = fgets_noempty(test_buf, 4999, fp);
                if(rptr[0] == '+')
                {
                    rptr = fgets_noempty(test_buf, 4999, fp);
                    if(rptr && second_line_len == strlen(test_buf))
                    ret = FILE_TYPE_FASTQ;
                }
            }
        }
        else if(nch == '>') // FASTA
        {
            char * rptr = fgets(test_buf, 4999, fp);
            int x1;
            if(rptr)
            {
                ret = FILE_TYPE_FASTA;
                for(x1=0;x1<4999;x1++)
                {
                    if(test_buf[x1]=='\n' || !test_buf[x1]) break;
                    nch = toupper(test_buf[x1]);
                    if(nch < ' ' || nch>127)
                    {
                        ret = FILE_TYPE_UNKNOWN;
                        break;
                    }
                }
                rptr = fgets(test_buf, 4999, fp);
                if(rptr && ret == FILE_TYPE_FASTA)
                {
                    for(x1=0;x1<4999;x1++)
                    {
                        if(test_buf[x1]=='\n' || !test_buf[x1]) break;
                        nch = toupper(test_buf[x1]);
                        if(nch == 'A' || nch == 'T' || nch == 'G' || nch == 'C' || nch == 'N' || nch == '.' || (nch >='0' && nch <= '3'))
                        ;
                        else
                        {
                            ret = FILE_TYPE_UNKNOWN;
                            break;
                        }
                    }
                    
                    if(x1==0) ret = FILE_TYPE_UNKNOWN;
                }
            }
        }
        else if(nch == 31) // BAM OR GZ_FASTQ
        {
            nch = fgetc(fp);
            if(nch == 139)
            {
                fclose(fp);
                fp=NULL;
                gzFile zfp = gzopen(fname, "rb");
                if(zfp)
                {
                    int rlen = gzread(zfp, test_buf,4);
                    if(rlen == 4 && memcmp(test_buf,"BAM\1",4)==0)
                    ret = FILE_TYPE_BAM;
                    if(rlen == 4 && test_buf[0]=='@')
                    ret = FILE_TYPE_GZIP_FASTQ;
                    if(rlen == 4 && test_buf[0]=='>')
                    ret = FILE_TYPE_GZIP_FASTA;
                    gzclose(zfp);
                }
            }
        }
        else if(nch >= 0x20 && nch <= 0x7f) // SAM without headers
        {
            int tabs = 0, x1;
            char * rptr = fgets(test_buf, 4999, fp);
            if(rptr)
            for(x1=0;x1<4999;x1++)
            {
                if(test_buf[x1]=='\n' || !test_buf[x1]) break;
                if(test_buf[x1]=='\t'){
                    tabs++;
                    continue;
                }
                if(tabs == 1)
                if(!isdigit(test_buf[x1]))break;
            }
            if(tabs>7)
            ret = FILE_TYPE_SAM;
            
        }
    }
    
    if(fp)fclose(fp);
    
    if(FILE_TYPE_BAM == ret || FILE_TYPE_SAM == ret)
    if(is_first_read_PE)
    {
        SamBam_FILE * tpfp = SamBam_fopen(fname, (FILE_TYPE_BAM  == ret)?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
        while(1)
        {
            char * tbr = SamBam_fgets(tpfp, test_buf, 4999, 0);
            if(!tbr){
                ret = FILE_TYPE_EMPTY;
                break;
            }
            if(tbr[0]=='@') continue;
            char * rname_str, *tmpstr;
            rname_str = strtok_r(tbr, "\t", &tmpstr);
            if(!rname_str)
            {
                ret = FILE_TYPE_UNKNOWN;
                break;
            }
            rname_str = strtok_r(NULL, "\t", &tmpstr);
            if((!rname_str)|| (!isdigit(rname_str[0])))
            {
                ret = FILE_TYPE_UNKNOWN;
                break;
            }
            
            int flags = atoi(rname_str);
            (*is_first_read_PE) = flags &1;
            break;
        }
        SamBam_fclose(tpfp);
    }
    
    free(test_buf);
    return ret;
}

int is_certainly_bam_file(char * fname, int * is_first_read_PE)
{
    
    int read_type = probe_file_type(fname, is_first_read_PE);
    if(read_type == FILE_TYPE_NONEXIST || read_type == FILE_TYPE_EMPTY || read_type == FILE_TYPE_UNKNOWN)
    return -1;
    if(read_type == FILE_TYPE_BAM)
    return 1;
    return 0;
}

void find_tag_out(char * read_line_buf, char * tag, char * hi_tag_out)
{
    int hi_tag = -1;
    char tag_str[10];
    sprintf(tag_str , "\t%s:i:", tag);
    char * hi_tag_str = strstr(read_line_buf, tag_str);
    if(hi_tag_str)
    {
        
        
        hi_tag = 0;
        int line_cursor;
        for(line_cursor=6; ; line_cursor++)
        {
            char nch = hi_tag_str[line_cursor];
            //								printf("HI:i=%s; nch [%d] ='%c'\n", hi_tag_str, line_cursor, nch);
            if(!isdigit(nch)) break;
            hi_tag = hi_tag*10 + (nch-'0');
        }
    }
    
    if(hi_tag >=0)
    {
        sprintf(hi_tag_out,"\t%s:i:%d", tag, hi_tag);
    }else hi_tag_out[0] = 0;
    
}

int sort_SAM_add_line(SAM_sort_writer * writer, char * SAM_line, int line_len)
{
    assert(writer -> all_chunks_header_fp);
    if(line_len<3) return 0;
    if(SAM_line[0]=='@')
    fputs(SAM_line, writer -> out_fp);
    else
    {
        char read_name[MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26];
        char chromosome_1_name[MAX_CHROMOSOME_NAME_LEN];
        char chromosome_2_name[MAX_CHROMOSOME_NAME_LEN];
        unsigned int pos_1, pos_2;
        int hi_tag,flags = 0, line_cursor = 0, field_cursor = 0, tabs=0;
        char * second_col_pos = NULL;
        
        chromosome_1_name[0]=0;
        chromosome_2_name[0]=0;
        pos_1 = 0;
        pos_2 = 0;
        hi_tag = -1;
        
        while(line_cursor < line_len)
        {
            char nch = SAM_line[line_cursor++];
            if(!nch)break;
            
            if(nch == '\t')
            {
                field_cursor = 0;
                tabs++;
                if(tabs == 1) second_col_pos = SAM_line + line_cursor;
                if(tabs>7) break;
            }
            else if(tabs == 0)
            {
                read_name[field_cursor++] = nch;
                if(MAX_READ_NAME_LEN<=field_cursor){
                    return -1;
                }
                read_name[field_cursor] = 0;
            }
            else if(tabs == 1)
            flags = flags*10+(nch-'0');
            else if(tabs == 2)
            {
                chromosome_1_name[field_cursor++] = nch;
                chromosome_1_name[field_cursor]=0;
                if(MAX_CHROMOSOME_NAME_LEN - 1 <= field_cursor) return -1;
            }
            else if(tabs == 3)
            pos_1 = pos_1 * 10 + (nch-'0');
            else if(tabs == 6)
            {
                chromosome_2_name[field_cursor++] = nch;
                chromosome_2_name[field_cursor] = 0;
                if(MAX_CHROMOSOME_NAME_LEN - 1 <= field_cursor) return -1;
            }
            else if(tabs == 7)
            pos_2 = pos_2 * 10 + (nch-'0');
            
        }
        if(tabs <= 7) return -1;
        
        char * hi_tag_str = strstr(SAM_line,"\tHI:i:");
        if(hi_tag_str)
        {
            hi_tag = 0;
            for(line_cursor=6; ; line_cursor++)
            {
                char nch = hi_tag_str[line_cursor];
                if(!isdigit(nch)) break;
                hi_tag = hi_tag*10 + (nch-'0');
            }
        }
        
        line_len = strlen(second_col_pos);
        sort_SAM_check_chunk(writer);
        
        for(field_cursor = 0; read_name[field_cursor] ; field_cursor++)
        if(read_name[field_cursor] == '/') read_name[field_cursor] = 0;
        
        if(chromosome_2_name[0]=='=')
        strcpy(chromosome_2_name, chromosome_1_name);
        
        
        // new read name format: OLD_READ_NAME\tCHR_R1:POS_R1:CHR_R2:POS_R2
        
        
        if(flags & SAM_FLAG_MATE_UNMATCHED)
        {
            if(chromosome_2_name[0] != '*')
            strcpy(chromosome_2_name , "*");
            pos_2 = 0;
        }
        
        
        if(flags & SAM_FLAG_UNMAPPED)
        {
            if(chromosome_1_name[0] != '*')
            strcpy(chromosome_1_name , "*");
            pos_1 = 0;
        }
        
        char hi_key [13];
        if(hi_tag >=0 && pos_1 && pos_2)
        sprintf(hi_key, ":%d", hi_tag);
        else
        hi_key[0]=0;
        
        if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
        sprintf(read_name+strlen(read_name), "\t%s:%u:%s:%u%s",chromosome_2_name, pos_2, chromosome_1_name, pos_1, hi_key);
        else
        sprintf(read_name+strlen(read_name), "\t%s:%u:%s:%u%s",chromosome_1_name, pos_1, chromosome_2_name, pos_2, hi_key);
        
        //if(memcmp("V0112_0155:7:1101:4561:132881", read_name, 27)==0)
        //	printf("RRN=%s\n", read_name);
        
        int read_name_len = strlen(read_name);
        unsigned long long int read_line_hash = sort_SAM_hash(read_name);
        
        int block_id = read_line_hash % SAM_SORT_BLOCKS;
        if(!writer -> current_block_fp_array[block_id])
        {
            char tmpfname[MAX_FILE_NAME_LENGTH+40];
            sprintf(tmpfname,"%sCHK%08d-BLK%03d.bin", writer -> tmp_path , writer -> current_chunk , block_id);
            writer -> current_block_fp_array[block_id] = fopen(tmpfname, "wb");
        }
        
        if(line_len < 2)
        {
            SUBREADprintf("unable to put the first read!\n");
            assert(0);
        }
        
        if(second_col_pos[0]==0 || second_col_pos[1]==0)
        {
            SUBREADprintf("unable to put the first read TEXT!\n");
            assert(0);
        }
        
        //		printf("WRNAME:%s\n", read_name);
        
        fwrite(&flags, 2, 1, writer -> current_block_fp_array[block_id]);
        fwrite(&read_name_len, 2, 1, writer -> current_block_fp_array[block_id]);
        fwrite(read_name, 1, read_name_len, writer -> current_block_fp_array[block_id]);
        fwrite(&line_len, 2, 1, writer -> current_block_fp_array[block_id]);
        fwrite(second_col_pos, 1, line_len, writer -> current_block_fp_array[block_id]);
        
        writer -> output_file_size += line_len;
        writer -> current_chunk_size += line_len;
        writer -> written_reads ++;
    }
    
    return 0;
}

void sort_SAM_finalise(SAM_sort_writer * writer)
{
    int x1_chunk, x1_block;
    int xk1;
    for(xk1=0;xk1<SAM_SORT_BLOCKS;xk1++)
    {
        if(writer -> current_block_fp_array[xk1])
        fclose(writer -> current_block_fp_array[xk1]);
    }
    memset(writer -> current_block_fp_array, 0, sizeof(FILE *)*SAM_SORT_BLOCKS);
    writer -> current_chunk_size = 0;
    writer -> current_chunk++;
    
    for(x1_block = 0; x1_block <SAM_SORT_BLOCKS; x1_block++){
        HashTable * first_read_name_table;
        first_read_name_table = HashTableCreate(SAM_SORT_BLOCK_SIZE / 100 );
        HashTableSetKeyComparisonFunction(first_read_name_table , fc_strcmp_chro);
        HashTableSetDeallocationFunctions(first_read_name_table , free, free);
        HashTableSetHashFunction(first_read_name_table, HashTableStringHashFunction);
        
        for(x1_chunk = 0; x1_chunk < writer -> current_chunk; x1_chunk++)
        {
            char tmpfname[MAX_FILE_NAME_LENGTH+40];
            sprintf(tmpfname, "%sCHK%08d-BLK%03d.bin", writer -> tmp_path, x1_chunk , x1_block);
            
            FILE * bbfp = fopen(tmpfname,"rb");
            if(!bbfp) continue;
            
            while(!feof(bbfp))
            {
                char * read_name = NULL;
                short flags;
                short read_name_len;
                short read_len;
                int ret = fread(&flags, 2,1 , bbfp);
                if(ret<1) break;
                fread(&read_name_len, 2,1 , bbfp);
                if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
                fseek(bbfp, read_name_len, SEEK_CUR);
                else
                {
                    read_name = malloc(read_name_len+1);
                    fread(read_name, 1, read_name_len, bbfp);
                    read_name[read_name_len] = 0;
                }
                fread(&read_len,2,1,bbfp);
                if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
                fseek(bbfp, read_len, SEEK_CUR);
                else
                {
                    char * new_line_mem = malloc(read_len+1);
                    fread(new_line_mem, 1, read_len, bbfp);
                    new_line_mem[read_len] = 0;
                    
                    if(read_len<2)
                    {
                        SUBREADprintf("Cannot determine read length from the tmp file!\n");
                        assert(0);
                    }
                    
                    
                    if( new_line_mem[0]==0 || new_line_mem[1]==0)
                    {
                        SUBREADprintf("Cannot load read part from the tmp file!\n");
                        assert(0);
                    }
                    
                    
                    char * old_line_mem = HashTableGet(first_read_name_table, read_name);
                    if(old_line_mem)
                    old_line_mem[0]=0xff;
                    else
                    HashTablePut(first_read_name_table, read_name, new_line_mem);
                    //if( first_read_name_table -> numOfElements<4)printf("RV=%s\n", read_name);
                }
            }
            
            fclose(bbfp);
        }
        
        //printf("BLK=%d; CKS=%d; READS=%llu\n", x1_block, x1_chunk, first_read_name_table -> numOfElements);
        unsigned long long int finished_second_reads = 0;
        
        for(x1_chunk = 0; x1_chunk < writer -> current_chunk; x1_chunk++)
        {
            char tmpfname[MAX_FILE_NAME_LENGTH+40];
            sprintf(tmpfname, "%sCHK%08d-BLK%03d.bin", writer -> tmp_path, x1_chunk , x1_block);
            
            //		printf("START_BLOCK: %s\n", tmpfname);
            
            FILE * bbfp = fopen(tmpfname,"rb");
            if(!bbfp) continue;
            
            char * read_line_buf = malloc(3000);
            char * read_name_buf = malloc(MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26);
            
            while(!feof(bbfp))
            {
                short flags;
                short read_name_len;
                short read_len;
                int ret = fread(&flags, 2,1 , bbfp);
                if(ret<1) break;
                
                fread(&read_name_len, 2,1 , bbfp);
                
                if(read_name_len>=MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26)
                SUBREADprintf("VERY_LONG_NAME(%d)\n", read_name_len);
                if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
                {
                    fread(read_name_buf, 1, read_name_len, bbfp);
                    read_name_buf[read_name_len] = 0;
                }
                else	fseek(bbfp, read_name_len, SEEK_CUR);
                fread(&read_len, 2,1 , bbfp);
                if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
                {
                    fread(read_line_buf, 1, read_len, bbfp);
                    read_line_buf[read_len] = 0;
                }
                else	fseek(bbfp, read_len, SEEK_CUR);
                
                
                if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
                {
                    //					printf("RRNAME:%s\n", read_name_buf);
                    
                    char * first_read_text = HashTableGet(first_read_name_table, read_name_buf);
                    strtok(read_name_buf,"\t");
                    if(first_read_text && first_read_text[0]!=(char)0xff)
                    {
                        fputs(read_name_buf, writer->out_fp);
                        putc('\t',  writer->out_fp);
                        fputs(first_read_text, writer->out_fp);
                        
                        fputs(read_name_buf, writer->out_fp);
                        putc('\t',  writer->out_fp);
                        fputs(read_line_buf, writer->out_fp);
                        
                        read_name_buf[strlen(read_name_buf)]='\t';
                        HashTableRemove(first_read_name_table, read_name_buf);
                        finished_second_reads ++;
                    }
                    else{
                        
                        int dummy_flags = 4 | 1, mate_flags = 0;
                        char * dummy_mate_chr = NULL;
                        char dummy_mate_chr_buf[120];
                        unsigned int dummy_mate_pos = 0, tmpi=0,dummy_char_strpos = 0;
                        int tabs = 0;
                        int read_cursor = 0;
                        
                        for(read_cursor = 0;; read_cursor++)
                        {
                            char nch = read_line_buf[read_cursor];
                            if(!nch) break;
                            if(nch == '\t')
                            {
                                if(tabs == 0){
                                    mate_flags = tmpi;
                                    dummy_mate_chr = read_line_buf+read_cursor+1;
                                }
                                else if(tabs == 1)
                                dummy_char_strpos = read_cursor;
                                else if(tabs == 2)
                                {
                                    dummy_mate_pos = tmpi;
                                    break;
                                }
                                tmpi=0;
                                tabs++;
                            }else{
                                if(tabs==0 || tabs == 2) tmpi = tmpi * 10 + (nch - '0');
                            }
                        }
                        
                        dummy_flags |= SAM_FLAG_FIRST_READ_IN_PAIR;
                        if(mate_flags & SAM_FLAG_UNMAPPED)  dummy_flags |= SAM_FLAG_MATE_UNMATCHED;
                        if(mate_flags & SAM_FLAG_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
                        if(mate_flags & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_REVERSE_STRAND_MATCHED;
                        
                        memcpy(dummy_mate_chr_buf, dummy_mate_chr, read_line_buf +dummy_char_strpos - dummy_mate_chr);
                        dummy_mate_chr_buf[read_line_buf +dummy_char_strpos - dummy_mate_chr]=0;
                        
                        char hi_tag_out[18];
                        char nh_tag_out[18];
                        
                        find_tag_out(read_line_buf, "HI", hi_tag_out);
                        find_tag_out(read_line_buf, "NH", nh_tag_out);
                        
                        // build a fake FIRST read for the mapped SECOND read.
                        // note that the TLEN, MATE_POS and MATE_CHAR are incorrect for general use.
                        fprintf(writer->out_fp, "%s\t%d\t*\t0\t0\t*\t%s\t%d\t0\tN\tI%s%s\n", read_name_buf, dummy_flags, dummy_mate_chr_buf, dummy_mate_pos, nh_tag_out, hi_tag_out);
                        fputs(read_name_buf, writer->out_fp);
                        putc('\t',  writer->out_fp);
                        fputs(read_line_buf, writer->out_fp);
                        writer -> unpaired_reads +=1;
                    }
                    
                    //else SUBREADprintf("WARNING: Unpaired read found in file:%s\n", read_name_buf);
                }
            }
            
            fclose(bbfp);
            unlink(tmpfname);
            free(read_name_buf);
            free(read_line_buf);
        }
        
        
        
        if(1)
        {
            writer -> unpaired_reads += first_read_name_table -> numOfElements;
            
            KeyValuePair * cursor;
            int bucket;
            
            // go through the hash table and write correct FIRST lines and dummy SECOND lines.
            for(bucket=0; bucket< first_read_name_table -> numOfBuckets; bucket++)
            {
                cursor = first_read_name_table -> bucketArray[bucket];
                while(1)
                {
                    if (!cursor) break;
                    char * first_read_text = (char *)cursor -> value;
                    char * first_read_name = (char *)cursor -> key;
                    
                    if(first_read_text[0]!=(char)0xff)
                    {
                        int dummy_flags = 4 | 1, mate_flags = 0;
                        char * dummy_mate_chr = NULL;
                        unsigned int dummy_mate_pos = 0, tmpi=0, dummy_char_strpos = 0;
                        int tabs = 0;
                        int read_cursor = 0;
                        
                        for(read_cursor = 0;; read_cursor++)
                        {
                            char nch = first_read_text[read_cursor];
                            if(!nch) break;
                            if(nch == '\t')
                            {
                                if(tabs == 0){
                                    mate_flags = tmpi;
                                    dummy_mate_chr = first_read_text+read_cursor+1;
                                }
                                else if(tabs == 1)
                                dummy_char_strpos = read_cursor;
                                else if(tabs == 2)
                                {
                                    dummy_mate_pos = tmpi;
                                    break;
                                }
                                tmpi=0;
                                tabs++;
                            }else{
                                if(tabs==0 || tabs == 2) tmpi = tmpi * 10 + (nch - '0');
                            }
                        }
                        
                        dummy_flags |= SAM_FLAG_SECOND_READ_IN_PAIR;
                        if(mate_flags & SAM_FLAG_UNMAPPED)  dummy_flags |= SAM_FLAG_MATE_UNMATCHED;
                        if(mate_flags & SAM_FLAG_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
                        if(mate_flags & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_REVERSE_STRAND_MATCHED;
                        
                        if((!first_read_text[0])||(!first_read_text[1]))
                        {
                            SUBREADprintf("unable to recover the first read! : '%s' , flags = %d\n", first_read_name, mate_flags);
                            assert(0);
                        }
                        
                        char nh_tag_out[18];
                        char hi_tag_out[18];
                        find_tag_out(first_read_text, "NH", nh_tag_out);
                        find_tag_out(first_read_text, "HI", hi_tag_out);
                        
                        strtok(first_read_name, "\t");
                        fputs(first_read_name, writer->out_fp);
                        putc('\t',  writer->out_fp);
                        fputs(first_read_text, writer->out_fp);
                        first_read_text[dummy_char_strpos] = 0;
                        fprintf(writer->out_fp, "%s\t%d\t*\t0\t0\t*\t%s\t%d\t0\tN\tI%s%s\n", first_read_name, dummy_flags, dummy_mate_chr, dummy_mate_pos, nh_tag_out,hi_tag_out);
                    }
                    cursor = cursor->next;
                }
            }
            
            
        }
        
        HashTableDestroy(first_read_name_table);
    }
    fclose(writer -> out_fp);
    signal (SIGTERM, old_sig_TERM);
    signal (SIGINT, old_sig_INT);
}

void sort_SAM_check_chunk(SAM_sort_writer * writer)
{
    if(writer -> current_chunk_size > SAM_SORT_BLOCK_SIZE * SAM_SORT_BLOCKS)
    {
        int xk1;
        for(xk1=0;xk1<SAM_SORT_BLOCKS;xk1++)
        {
            if(writer -> current_block_fp_array[xk1])
            fclose(writer -> current_block_fp_array[xk1]);
        }
        memset(writer -> current_block_fp_array, 0, sizeof(FILE *)*SAM_SORT_BLOCKS);
        writer -> current_chunk_size = 0;
        writer -> current_chunk++;
    }
}

unsigned long long int sort_SAM_hash(char * str)
{
    unsigned long long int hash = 5381;
    int c, xk1=0;
    
    while (1)
    {
        c = str[xk1++];
        if(!c)break;
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    return hash;
}

int is_pipe_file(char * fname)
{
    FILE * fp = fopen(fname,"r");
    if(!fp) return 0;
    
    int seeked = fseek(fp, 0, SEEK_SET);
    fclose(fp);
    
    return (seeked != 0);
}

int warning_file_type(char * fname, int expected_type)
{
    int ret_pipe_file = is_pipe_file(fname);
    if(ret_pipe_file)
    {
        SUBREADprintf("WARNING file '%s' is not a regular file.\n", fname);
        SUBREADprintf("        No alignment can be done on a pipe file.\n");
        SUBREADprintf("        If the FASTQ file is gzipped, please use gzFASTQinput option.\n");
        return 1;
    }
    
    int read_type = probe_file_type(fname, NULL);
    
    if(read_type == FILE_TYPE_NONEXIST)
    {
        SUBREADprintf("ERROR: unable to open file '%s'. File name might be incorrect, or you do not have the permission to read the file.\n", fname);
        return -1;
    }
    else if(read_type == FILE_TYPE_EMPTY)
    {
        SUBREADprintf("WARNING file '%s' is empty.\n", fname);
        return 1;
    }
    
    else if((expected_type == FILE_TYPE_FAST_ && (read_type!= FILE_TYPE_FASTQ && read_type!= FILE_TYPE_FASTA))||
            (expected_type == FILE_TYPE_GZIP_FAST_ && (read_type!= FILE_TYPE_GZIP_FASTQ && read_type!= FILE_TYPE_GZIP_FASTA)) ||
            ((  expected_type != FILE_TYPE_GZIP_FAST_ && expected_type != FILE_TYPE_FAST_) && expected_type != read_type))
    {
        char * req_fmt = "SAM";
        if(expected_type==FILE_TYPE_BAM) req_fmt = "BAM";
        else if(expected_type==FILE_TYPE_FAST_) req_fmt = "FASTQ or FASTA";
        else if(expected_type==FILE_TYPE_GZIP_FAST_) req_fmt = "gzip FASTQ or FASTA";
        
        char * real_fmt = "SAM";
        if(read_type==FILE_TYPE_BAM) real_fmt = "BAM";
        else if(read_type==FILE_TYPE_FASTA) real_fmt = "FASTA";
        else if(read_type==FILE_TYPE_FASTQ) real_fmt = "FASTQ";
        else if(read_type==FILE_TYPE_GZIP_FASTQ) real_fmt = "gzip FASTQ";
        else if(read_type==FILE_TYPE_GZIP_FASTA) real_fmt = "gzip FASTA";
        
        SUBREADprintf("WARNING format issue in file '%s':\n", fname);
        SUBREADprintf("        The required format is : %s\n", req_fmt);
        if(read_type == FILE_TYPE_UNKNOWN)
        SUBREADprintf("        The file format is unknown.\n");
        else
        SUBREADprintf("        The real format seems to be : %s\n", real_fmt);
        SUBREADprintf("A wrong format may result in wrong results or crash the program.\n");
        SUBREADprintf("Please refer to the manual for file format options.\n");
        SUBREADprintf("If the file is in the correct format, please ignore this message.\n");
        
        return 1;
    }
    return 0;
}
