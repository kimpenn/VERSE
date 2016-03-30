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
  
#include <string.h>
#include <stdarg.h>
#include <unistd.h>

#ifndef FREEBSD
#include <sys/timeb.h>
#endif

#include "HelperFunctions.h"


// original core functions
double miltime(){
    double ret;
#ifdef FREEBSD
    struct timeval tp;
    struct timezone tz;
    tz.tz_minuteswest=0;
    tz.tz_dsttime=0;
    
    gettimeofday(&tp,&tz);
    
    ret = tp.tv_sec+ 0.001*0.001* tp.tv_usec;
    
#else
    
    struct timeb trp;
    ftime(&trp);
    ret = trp.time*1.0+(trp.millitm*1.0/1000.0);
#endif
    
    return ret;
}

int fc_strcmp_chro(const void * s1, const void * s2)
{
    
    //	// we can compare the 3-th and 4-th bytes because we know that the buffers have enough lengths.
    //	if(((char *)s1)[4] != ((char *)s2)[4] || ((char *)s1)[3] != ((char *)s2)[3] )
    //		return 1;
    return strcmp((char*)s1, (char*)s2);
}

unsigned long fc_chro_hash(const void *key) {
    const unsigned char *str = (const unsigned char *) key;
    
    int xk1;
    unsigned long hashValue = 0;
    
    for(xk1=0; str[xk1]; xk1++)
    {
        unsigned long ch = str[xk1];
        hashValue += (ch + xk1) << (ch & 0xf);
    }
    
    
    return hashValue;
}

void print_in_box(int line_width, int is_boundary, int is_center, char * pattern,...)
{
    va_list args;
    va_start(args , pattern);
    char is_R_linebreak=0, * content, *out_line_buff;
    
    content= malloc(1000);
    out_line_buff= malloc(1000);
    out_line_buff[0]=0;
    vsprintf(content, pattern, args);
    int x1,content_len = strlen(content), state, txt_len, is_cut = 0, real_lenwidth;
    
    if(content_len>0&&content[content_len-1]=='\r'){
        content_len--;
        content[content_len] = 0;
        is_R_linebreak = 1;
    }
    
    if(content_len>0&&content[content_len-1]=='\n'){
        content_len--;
        content[content_len] = 0;
    }
    
    state = 0;
    txt_len = 0;
    real_lenwidth = line_width;
    for(x1 = 0; content [x1]; x1++)
    {
        char nch = content [x1];
        if(nch == CHAR_ESC)
            state = 1;
        if(state){
            real_lenwidth --;
        }else{
            txt_len++;
            
            if(txt_len == 80 - 6)
            {
                is_cut = 1;
            }
        }
        
        if(nch == 'm' && state)
            state = 0;
    }
    
    if(is_cut)
    {
        state = 0;
        txt_len = 0;
        for(x1 = 0; content [x1]; x1++)
        {
            char nch = content [x1];
            if(nch == CHAR_ESC)
                state = 1;
            if(!state){
                txt_len++;
                if(txt_len == 80 - 9)
                {
                    strcpy(content+x1, "\x1b[0m ...");
                    content_len = line_width - 4;
                    content_len = 80 - 4;
                    line_width = 80;
                    break;
                }
            }
            if(nch == 'm' && state)
                state = 0;
        }
    }
    
    if(content_len==0 && is_boundary)
    {
        strcat(out_line_buff,is_boundary==1?"//":"\\\\");
        for(x1=0;x1<line_width-4;x1++)
            strcat(out_line_buff,"=");
        strcat(out_line_buff,is_boundary==1?"\\\\":"//");
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);
        
        free(content);
        free(out_line_buff);
        return;
    }
    else if(is_boundary)
    {
        int left_stars = (line_width - content_len)/2 - 1;
        int right_stars = line_width - content_len - 2 - left_stars;
        strcat(out_line_buff,is_boundary==1?"//":"\\\\");
        for(x1=0;x1<left_stars-2;x1++) strcat(out_line_buff,"=");
        sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
        sprintf(out_line_buff+strlen(out_line_buff)," %s ", content);
        sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
        for(x1=0;x1<right_stars-2;x1++) strcat(out_line_buff,"=");
        strcat(out_line_buff,is_boundary==1?"\\\\":"//");
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);
        
        free(content);
        free(out_line_buff);
        return;
    }
    
    int right_spaces, left_spaces;
    if(is_center)
        left_spaces = (line_width - content_len)/2-2;
    else
        left_spaces = 1;
    
    right_spaces = line_width - 4 - content_len- left_spaces;
    
    char spaces[81];
    memset(spaces , ' ', 80);
    spaces[0]='|';
    spaces[1]='|';
    spaces[80]=0;
    
    //sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"||");
    
    //for(x1=0;x1<left_spaces;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
    
    spaces[left_spaces+2] = 0;
    strcat(out_line_buff,spaces);
    
    
    int col1w=-1;
    for(x1=0; content[x1]; x1++)
    {
        if(content[x1]==':')
        {
            col1w=x1;
            break;
        }
    }
    if(col1w>0 && col1w < content_len-1)
    {
        content[col1w+1]=0;
        strcat(out_line_buff,content);
        strcat(out_line_buff," ");
        sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
        strcat(out_line_buff,content+col1w+2);
        sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
    }
    else
        strcat(out_line_buff,content);
    //	for(x1=0;x1<right_spaces - 1;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
    
    memset(spaces , ' ', 80);
    spaces[79]='|';
    spaces[78]='|';
    
    right_spaces = max(1,right_spaces);

    sprintf(out_line_buff+strlen(out_line_buff)," %c[0m%s", CHAR_ESC , spaces + (78 - right_spaces + 1));
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, out_line_buff);
    free(out_line_buff);
    free(content);
}

void merge_sort_run(void * arr, int start, int items, int index, int compare (void * arr, int l, int r, int index), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2, int index))
{
    if(items > 11)
    {
        int half_point = items/2;
        merge_sort_run(arr, start, half_point, index, compare, exchange, merge);
        merge_sort_run(arr, start + half_point, items - half_point, index, compare, exchange, merge);
        merge(arr, start, half_point, items - half_point, index);
    }
    else
    {
        int i, j;
        for(i=start; i < start + items - 1; i++)
        {
            int min_j = i;
            for(j=i + 1; j< start + items; j++)
            {
                if(compare(arr, min_j, j, index) > 0)
                    min_j = j;
            }
            if(i!=min_j)
                exchange(arr, i, min_j);
        }
    }
}
void merge_sort(void * arr, int arr_size, int index, int compare (void * arr, int l, int r, int index), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2, int index))
{
    merge_sort_run(arr, 0, arr_size, index, compare, exchange, merge);
}

void print_verse_logo()
{
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m ===      === %c[0m%c[34m   %%%%  %%%%  %%%%%%%%%%%%  %%%%   %%%%  %%%%       %%%%%%%%   %%%%%%%%%%   ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m  ==      ==  %c[0m%c[34m   %%%% %%%%     %%%%    %%%%%% %%%%%%  %%%%      %%%%  %%%%  %%%%  %%%%  ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m   ==    ==   %c[0m%c[34m   %%%%%%%%      %%%%    %%%% %% %%%%  %%%%      %%%%%%%%%%%%  %%%%%%%%%%   ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m    ==  ==    %c[0m%c[34m   %%%% %%%%     %%%%    %%%%   %%%%  %%%%      %%%%  %%%%  %%%%  %%%%  ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m     ====     %c[0m%c[34m   %%%%  %%%%  %%%%%%%%%%%%  %%%%   %%%%  %%%%%%%%%%%%  %%%%  %%%%  %%%%%%%%%%   ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[40;33m      ==      %c[0m%c[34m  ................................................. %c[0m", CHAR_ESC, CHAR_ESC, CHAR_ESC, CHAR_ESC);
    char * spaces = "";
    if(strlen(VERSE_VERSION) == 8) spaces = "";
    else if(strlen(VERSE_VERSION) == 5) spaces = "  ";
    sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"	%sv%s",spaces,VERSE_VERSION);
}

int term_strncpy(char * dst, char * src, int max_dst_mem)
{
    int i;
    
    for(i=0; i<max_dst_mem; i++)
    {
        if(!src[i]) break;
        dst[i]=src[i];
        if(i == max_dst_mem-1)
            SUBREADprintf("String out of memory limit: '%s'\n", src);
    }
    if(i >= max_dst_mem) i = max_dst_mem-1;
    dst[i] = 0;
    
    return 0;
}

void core_version_number(char * program)
{
    VERSEprintf("\n%s v%s\n\n" , program, VERSE_VERSION);
}

// sublog
void remove_ESC_effects(char * txt)
{
    int x1;
    int ocur = 0;
    int state = 0;
    int trimmed = 0;
    
    //	return 0;
    
    for(x1=0;x1<1199;x1++)
    {
        if(!txt[x1])break;
        if((state == 0) && (txt[x1]==CHAR_ESC))
        {
            state = 1;
            trimmed = 1;
            continue;
        }
        if(state == 0){
            if(x1>ocur)
            txt[ocur] = txt[x1];
            ocur++;
        }
        
        if((state == 1) && (txt[x1]=='m'))
        state = 0;
    }
    if(trimmed)
    txt[ocur]=0;
}

int is_ESC_removed()
{
    return !isatty(fileno(stdout));    
}

void sublog_printf(int stage, int level, const char * pattern, ...)
{
    va_list args;
    va_start(args , pattern);
    if(level<MINIMUM_LOG_LEVEL) return;
    
    int to_remove_ESC = 1;
    to_remove_ESC = is_ESC_removed();
    
    if(to_remove_ESC)
    {
        char * vsbuf=malloc(1200);
        
        vsnprintf(vsbuf, 1199, pattern , args);
        remove_ESC_effects(vsbuf);
        
        fputs(vsbuf,stdout);
        fputs("\n", stdout);
        
        free(vsbuf);
    }
    else
    {
        vfprintf(stdout, pattern , args);
        va_end(args);
        fputs("\n", stdout);
        
        fflush(stdout);
    }
}

// original helper functions
int RSubread_parse_CIGAR_string(const char * CIGAR_Str, unsigned int * Staring_Points, unsigned short * Section_Length, int * is_N)
{
	unsigned int tmp_int=0;
	int cigar_cursor=0;
	unsigned short read_cursor=0;
	unsigned int chromosome_cursor=0;
	int ret=0;

	for(cigar_cursor=0; ; cigar_cursor++)
	{
		char ch = CIGAR_Str[cigar_cursor];

		if(ch >='0' && ch <= '9')
		{
			tmp_int=tmp_int*10+(ch - '0');
		}
		else
		{
			if(ch == 'M')
			{
				read_cursor += tmp_int;
				chromosome_cursor += tmp_int;
			}
			else if(ch == 'N' || ch == 'D' || ch == 0)
			{
                if(ch == 'N')
                    *is_N = 1;
				if(ret <6)
				{
					if(read_cursor>0)
					{
						Staring_Points[ret] = chromosome_cursor - read_cursor;
						Section_Length[ret] = read_cursor;
						ret ++;
					}
				}
				read_cursor = 0;

				if(ch == 'N' || ch == 'D') chromosome_cursor += tmp_int;
				else break;
			}
			//printf("C=%c, TV=%d, CC=%d, RC=%d\n", ch, tmp_int, chromosome_cursor, read_cursor);
			tmp_int = 0;
		}
		if(cigar_cursor>100) return -1;
	}

	return ret;
}


#define GECV_STATE_BEFORE 10
#define GECV_STATE_NAME 20
#define GECV_STATE_GAP 30
#define GECV_STATE_VALUE 40
#define GECV_STATE_QVALUE 50
#define GECV_STATE_QV_END 60
#define GECV_STATE_ERROR 9999

int GTF_extra_column_istoken_chr(char c)
{
	return (isalpha(c)||isdigit(c)||c=='_');
}

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size)
{
	int state = GECV_STATE_BEFORE;
	int col_cursor = 0, is_correct_name=0;
	char name_buffer[200];
	int name_cursor = 0, value_cursor=-1;

	while(1)
	{
		if(name_cursor>190) return -1;
		char nch = Extra_Col[col_cursor];
		if(nch == '\n' || nch == '\r') nch = 0;
		if(state == GECV_STATE_BEFORE)
		{
			if(GTF_extra_column_istoken_chr(nch))
			{
				name_buffer[0] = nch;
				name_cursor = 1;
				state = GECV_STATE_NAME;
			}
			else if(nch != ' ' && nch != 0)
			{
				state = GECV_STATE_ERROR;
			}
		}
		else if(state == GECV_STATE_NAME)
		{
			if(nch == ' ' || nch == '=')
			{
				state = GECV_STATE_GAP;
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				//printf("CORR=%d : '%s'\n", is_correct_name, name_buffer);
			}
			else if(nch == '"')
			{
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(GTF_extra_column_istoken_chr(nch))
				name_buffer[name_cursor++] = nch;
			else
			{
				state = GECV_STATE_ERROR;
				//printf("ERR2  : '%c'\n", nch);
			}
			
		}
		else if(state == GECV_STATE_GAP)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(nch != '=' && isgraph(nch))
			{
				state = GECV_STATE_VALUE;
				if(is_correct_name)
				{
					Target_Value[0]=nch;
					value_cursor = 1;
				}
			}
			else if(nch != ' ' && nch != '=')
				state = GECV_STATE_ERROR;
		}
		else if(state == GECV_STATE_VALUE)
		{
			if(nch == ';' || nch == 0)
			{
				state = GECV_STATE_BEFORE;
				if(is_correct_name)
				{
					Target_Value[value_cursor] = 0;
				}
				is_correct_name = 0;
			}
			else{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
					Target_Value[value_cursor++] = nch;
			}
		}
		else if(state == GECV_STATE_QVALUE)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QV_END;
				if(is_correct_name)
					Target_Value[value_cursor] = 0;
				is_correct_name = 0;
			}
			else
			{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
				{
					if(nch !=' ' || value_cursor>0)
						Target_Value[value_cursor++] = nch;
				}
			}
		}
		else if(state == GECV_STATE_QV_END)
		{
			if(nch == ';' || nch == 0)
				state = GECV_STATE_BEFORE;
			else if(nch != ' ')
				state = GECV_STATE_ERROR;
				
		}

		if (GECV_STATE_ERROR == state){
			Target_Value[0]=0;
			return -1;
		}
		if (nch == 0)
		{
			if(state == GECV_STATE_BEFORE && value_cursor>0)
			{
				int x1;
				for(x1 = value_cursor-1; x1>=0; x1--)
				{
					if(Target_Value[x1] == ' '){
						value_cursor --;
						Target_Value[x1]=0;
					}
					else break;
				}

				if(value_cursor>0)
					return value_cursor;
			}
			Target_Value[0]=0;
			return -1;
		}
		col_cursor++;
	}
}

char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; NULL != (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}

// rule: the string is ABC123XXXXXX...
// This is the priroity:
// First, compare the letters part.
// Second, compare the pure numeric part.
// Third, compare the remainder.
int strcmp_number(char * s1, char * s2)
{
	int x1 = 0;
	int ret = 0;

	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];

		if((!nch1) || !nch2){return nch2?1:(nch1?(-1):0);}
		if(isdigit(nch1) && isdigit(nch2))break;

		ret = nch1 - nch2;
		if(ret) return ret;
		x1++;
	}

	int v1 = 0, v2 = 0;
	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];
		if((!nch1) || !nch2){
			if(nch1 || nch2)
				return nch2?(-1):1;
			break;
		}
		int is_chr1_digit = isdigit(nch1);
		int is_chr2_digit = isdigit(nch2);

		if(is_chr1_digit || is_chr2_digit)
		{
			if(is_chr1_digit && is_chr2_digit)
			{
				v1 = v1*10+(nch1-'0');
				v2 = v2*10+(nch2-'0');
			}
			else
			{
				ret = nch1 - nch2;
				return ret;
			}
		}
		else break;
		x1++;
	}

	if(v1==v2)
		return strcmp(s1+x1, s2+x1);	
	else
	{
		ret = v1 - v2;
		return ret;
	}
}
