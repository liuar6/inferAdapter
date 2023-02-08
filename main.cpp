/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <map>
#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <htslib/sam.h>
#include <iomanip>
#include <algorithm>
using namespace std;
#include "tools.hpp"
struct read_sequence{
    char* seq;
    char* clip_seq;
    char* seed_seq;
    char* seq_end;
    bool unused;
};
int get_clip_len(bam1_t* bam){
    if (bam->core.flag & BAM_FUNMAP) return 0;
    if (bam_is_rev(bam)) return ((bam_cigar_op(bam_get_cigar(bam)[0])==4)?bam_cigar_oplen(bam_get_cigar(bam)[0]):0);
    else return ((bam_cigar_op(bam_get_cigar(bam)[bam->core.n_cigar-1])==4)?bam_cigar_oplen(bam_get_cigar(bam)[bam->core.n_cigar-1]):0);
}
char* get_seq(bam1_t* bam, char* seq, int seq_len) {
    int i;
    for (i = 0; i < bam->core.l_qseq && i < seq_len - 1; i++)
        seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(bam), i)];
    seq[i] = '\0';
    if (bam_is_rev(bam)) _reverse(_complement(seq));
    return seq;
}
int seed_len=8;
bool seed_comp(struct read_sequence *a, struct read_sequence *b){
    return(strncmp(a->seed_seq, b->seed_seq, seed_len)>0);
}
struct parameter_info{
    const char *bam_file;
    const char *out_file;
    uint32_t bin_size;
    uint32_t seed_step;
    uint32_t min_seed_count;
    uint32_t min_extend_count;
    uint32_t min_adapter_len;
    double major_base_proportion_cutoff;
    uint32_t mask_min_match;
    int32_t mask_search_start_offset;
    int32_t mask_search_end_offset;
    uint32_t max_try;
    uint32_t max_hit;

    FILE* out;
    FILE* out1;
};
char get_extended_base(vector<read_sequence*> *seqs, int position, struct parameter_info *parameters){
    double major_base_proportion=parameters->major_base_proportion_cutoff;
    int valid_count=0, total_count=0;
    int count_matrix[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    char base;
    for (int i = 0; i < seqs->size(); ++i) {
        if (seqs->at(i)->unused || seqs->at(i)->clip_seq + position >= seqs->at(i)->seq_end|| *(seqs->at(i)->clip_seq + position)=='N') continue;
            count_matrix[seq_nt16_table[*(seqs->at(i)->clip_seq + position)]] += 1;
            valid_count++;
    }
    fprintf(parameters->out, "%d\t%d\t%d\t%d\n", count_matrix[1], count_matrix[2], count_matrix[4], count_matrix[8]);
    if (valid_count>=parameters->min_extend_count){
        total_count=valid_count;
            if ((double)count_matrix[1]/total_count>=major_base_proportion) base='A';
            else if ((double)count_matrix[2]/total_count>=major_base_proportion) base='C';
            else if ((double)count_matrix[4]/total_count>=major_base_proportion) base='G';
            else if ((double)count_matrix[8]/total_count>=major_base_proportion) base='T';
            else base='N';
    } else base='?';
    return base;
};
void search_adapter(vector<read_sequence*> *seqs, char* adapter,  struct parameter_info *parameters){
    char best_seed_seq[200]="\0", seed_seq[200]="\0";
    int best_seed_seq_count=0, seed_seq_count=0;
    int best_seed_start=0, seed_start=0;
    for (int i=0; i<seqs->size(); ++i) seqs->at(i)->seed_seq=seqs->at(i)->clip_seq; /* reset seed position*/
    for (int i=0; i<1000; ++i) adapter[i]='\0'; /* reset adapter */
    int valid_count=parameters->min_seed_count;
    while(valid_count>=parameters->min_seed_count){
        valid_count=0;
        sort(seqs->begin(), seqs->end(), seed_comp);
        for (int i=0; i<seqs->size(); ++i){
            if(seqs->at(i)->seq_end-seqs->at(i)->seed_seq < seed_len) continue;
            if (strncmp(seed_seq, seqs->at(i)->seed_seq, seed_len)!=0){
                if (seed_seq_count>best_seed_seq_count && strchr(seed_seq, 'N')==NULL){
                    strncpy(best_seed_seq, seed_seq, seed_len);
                    best_seed_seq_count=seed_seq_count;
                    best_seed_start=seed_start;
                }
                strncpy(seed_seq, seqs->at(i)->seed_seq, seed_len);
                seed_seq_count=0;
            }
            seed_seq_count++;
            valid_count++;
        }
        if (seed_seq_count>best_seed_seq_count && strchr(seed_seq, 'N')==NULL){
            strncpy(best_seed_seq, seed_seq, seed_len);
            best_seed_seq_count=seed_seq_count;
            best_seed_start=seed_start;
        }
        seed_seq[0]='\0';
        seed_seq_count=0;
        for (int i=0; i<seqs->size(); ++i) seqs->at(i)->seed_seq+=parameters->seed_step;
        seed_start+=parameters->seed_step;
    }
    fprintf(parameters->out, "#seed status\nseed:      \t%s\nseed count:\t%d\nseed pos:  \t%d\n", best_seed_seq, best_seed_seq_count, best_seed_start);
    if (best_seed_seq_count<parameters->min_seed_count) {
        if (parameters->out1) fprintf(parameters->out1, "#stop for seed count too small.\n");
        return;
    }
    /* extended the seed to get full length adapter */
    int adapter_start=best_seed_start;
    int adapter_end=best_seed_start+seed_len-1;
    strcpy(adapter+adapter_start, best_seed_seq);
    for (int i=0; i<seqs->size(); ++i) seqs->at(i)->unused=(seqs->at(i)->clip_seq+adapter_end >= seqs->at(i)->seq_end || strncmp(seqs->at(i)->clip_seq+adapter_start, adapter+adapter_start, adapter_end-adapter_start+1)!=0);
    if (adapter_start>0) fprintf(parameters->out, "\n#upstream extension\nA\tC\tG\tT\n");
    while(adapter_start > 0) {
        char base=get_extended_base(seqs, --adapter_start, parameters);
        adapter[adapter_start]=base;
        if (base!='N') for (int i=0; i<seqs->size(); ++i)
            if (!seqs->at(i)->unused && *(seqs->at(i)->clip_seq+adapter_start)!=base) seqs->at(i)->unused=true;
    }
    fprintf(parameters->out, "\n#downstream extension\nA\tC\tG\tT\n");
    while(adapter_end < 1000){
        char base=get_extended_base(seqs, ++adapter_end, parameters);
        adapter[adapter_end]=base;
        if (base!='N') for (int i=0; i<seqs->size(); ++i)
            if (!seqs->at(i)->unused && *(seqs->at(i)->clip_seq+adapter_end)!=base) seqs->at(i)->unused=true;
        if (base=='?') break;
    }
    /* clip trailing 'N's and '?'*/
    for (char *p=adapter+strlen(adapter)-1; p>adapter; --p) if ((*p=='N' || *p=='?') && (*(p-1)=='N' || *(p-1)=='?')) *p='\0'; else break;
}
/* mask adapter sequences with 'N', currently mismatch is not supported.*/
int mask_adapter(vector<read_sequence*> *seqs, char* _adapter, int search_start_offset, int search_end_offset, long min_match, struct parameter_info *parameters){
    char bases[1000];
    char *adapter=bases;
    strcpy(adapter, _adapter);
    char *p=adapter+strlen(adapter)-1;
    /* ignore leading and trailing special characters*/
    while (*p=='N' || *p=='?') *(p--)='\0';
    while (*adapter=='N' || *adapter=='?') ++adapter;
    int adapter_len=strlen(adapter);
    int search_start=0, search_end=0;
    int match_count=0, count=0;
    int i, j, k;
    for (i=0; i<seqs->size(); ++i){
        char *seq=seqs->at(i)->seq;
        int seq_len=seqs->at(i)->seq_end-seqs->at(i)->seq;
        /* search start and search end, the first and last index to start a search for adapter */
        search_start=max(0L, seqs->at(i)->clip_seq-seqs->at(i)->seq+search_start_offset);
        search_end=min(seq_len-min_match, seqs->at(i)->clip_seq-seqs->at(i)->seq+search_end_offset);
        for (j=search_start; j<=search_end; ++j){
            match_count=0;
            for (k=0; j+k<seq_len && k<adapter_len; ++k){
                if (adapter[k]!='N'){
                    if (seq[j+k]==adapter[k]) match_count++;
                    else break;
                }
            }
            if ((j+k==seq_len || k==adapter_len) && match_count>=min_match){
                /* currently, bases are not to be masked if the corresponding position in adapter is 'N'. */
                for(int l=0 ; l<k; ++l) if (adapter[l]!='N') seq[j+l]='N';
                ++count;
                break;
            }
        }
    }
    return count;
}
void search_many_adapters(vector<read_sequence*> *seqs, struct parameter_info *parameters){
    char adapter[1000];
    int adapter_count=0;
    int i=0, hit=0;
    uint32_t n_base=0;
    for (i=0; i<parameters->max_try && hit<parameters->max_hit; ++i) {
        fprintf(parameters->out, "== Round %d ==\n", i + 1);
        search_adapter(seqs, adapter, parameters);
        if (strlen(adapter) == 0) {
            fprintf(parameters->out, "\n\n");
            break;
        }
        n_base = 0;
        for (int j = 0; adapter[j] != '\0'; ++j) if (adapter[j] != 'N' && adapter[j] != '?') ++n_base;
        /* is it nessesary to fix the search position at where the adapter is found?*/
        adapter_count = mask_adapter(seqs, adapter, parameters->mask_search_start_offset,
                                     parameters->mask_search_end_offset, min(n_base, parameters->mask_min_match),
                                     parameters);
        fprintf(parameters->out, "\n#result\n");
        if (n_base >= parameters->min_adapter_len) {
            hit++;
            fprintf(parameters->out, "seq     :\t%s\n", adapter);
            fprintf(parameters->out, "count   :\t%d\n", adapter_count);
            fprintf(parameters->out, "prop    :\t%.4f\n", (double)adapter_count/(double)seqs->size());
        } else fprintf(parameters->out, "the suggested adapter is too short.\n");
        fprintf(parameters->out, "\n\n");
    }
    if (i>=parameters->max_try && parameters->max_try>1) if (parameters->out1) fprintf(parameters->out1, "#stop search for reaching predefined times of tries.\n");
    else if (hit>=parameters->max_hit && parameters->max_hit>1) if (parameters->out1) fprintf(parameters->out1, "#stop search for reaching predefined number of hits.\n");
}
void usage(){
    fprintf(stderr, "%s", "inferAdapter: infer sequencing adapters from alignment file.\n\
Usage:  inferAdapter [options] -i <alignment file>\n\
[options]\n\
-i/--input                     : input bam file.[required]\n\
-o/--output                    : output text file.[required]\n\
-v/--version                   : show version.\n\
-h/--help                      : show help informations.\n\
-b/--bin-size                  : for each 2^(bin size) genomic region, one alignment will be used for search, default 8.\n\
-l/--min-adapter-len           : minimal adapter length required, default 12.\n\
-p/--base-proportion           : minimal base proportion for determining the base, default 0.9.\n\
--seed-len                     : length of seed to look for initial sequence for extension, default 8.\n\
--seed-step                    : shift of position to look for seed, default 3.\n\
--min-seed-count               : minimal count of seed sequences to start extension step, default 1000. \n\
--min-extend-count             : minimal count of available reads to continue extension step, default 300. \n\
--max-try                      : max times of tries to search for adapters, default 5.\n\
--max-hit                      : number of adapter required, default 3.\n\
");
    exit(1);
}
void parse_args(int argc, char *argv[], struct parameter_info *parameters){

    parameters->bam_file=nullptr;
    parameters->out_file=nullptr;
    parameters->bin_size=8;
    parameters->seed_step=3;
    seed_len=8;
    parameters->min_seed_count=1000;
    parameters->min_extend_count=300;
    parameters->min_adapter_len=12;
    parameters->major_base_proportion_cutoff=0.9;
    parameters->mask_min_match=8;
    parameters->mask_search_start_offset=-10;
    parameters->mask_search_end_offset=200;
    parameters->max_try=5;
    parameters->max_hit=3;

    /* special case */
    if (argc>=2 && argv[1][0]!='-'){
        parameters->bam_file=argv[1];
        --argc;
        ++argv;
    } else if (argc==1) usage();

    char c;
    int showHelp=0;
    int showVersion=0;
    const char* version="1.0.0";
    const char *shortOptions = "hvo:i:b:s:e:l:p:t:m:A:B:C:D:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "output" , required_argument , NULL, 'o' },
                    { "input" , required_argument, NULL, 'i' },
                    { "bin-size" , required_argument, NULL, 'b' },
                    { "seed-len" , required_argument, NULL, 'A' },
                    { "seed-step" , required_argument, NULL, 'B' },
                    { "min-seed-count" , required_argument, NULL, 's' },
                    { "min-extend-count" , required_argument, NULL, 'e' },
                    { "min-adapter-len" , required_argument, NULL, 'l' },
                    { "base-proportion" , required_argument, NULL, 'p' },
                    { "mask-min-match" , required_argument, NULL, 'C' },
                    { "mask-search-start" , required_argument, NULL, 'D' },
                    { "max-try" , required_argument, NULL, 't' },
                    { "max-hit" , required_argument, NULL, 'm' },
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };

    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                showHelp = 1;
                break;
            case 'v':
                showVersion = 1;
                break;
            case 'o':
                parameters->out_file=optarg;
                break;
            case 'i':
                parameters->bam_file=optarg;
                break;
            case 'b':
                parameters->bin_size=max(0L, strtol(optarg, nullptr, 10));
                break;
            case 'A':
                seed_len=max(1L, strtol(optarg, nullptr, 10)); /* global variable*/
                break;
            case 'B':
                parameters->seed_step=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 's':
                parameters->min_seed_count=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 'e':
                parameters->min_extend_count=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 'l':
                parameters->min_adapter_len=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 'p':
                parameters->major_base_proportion_cutoff=max(0.5, strtod(optarg, nullptr));
                break;
            case 'C':
                parameters->mask_min_match=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 'D':
                parameters->mask_search_start_offset=strtol(optarg, nullptr, 10);
                break;
            case 't':
                parameters->max_try=max(1L, strtol(optarg, nullptr, 10));
                break;
            case 'm':
                parameters->max_hit=max(1L, strtol(optarg, nullptr, 10));
                break;
            case '?':
                showHelp = 1;
                break;
            default:
                usage();
        }
    }

    if (argc != optind) usage();
    if (showVersion)
    {
        fprintf(stderr, "inferAdapter-%s\n", version);
        exit(1);
    }

    if (showHelp)
    {
        usage();
        exit(1);
    }
    if (parameters->out_file==nullptr || strcmp(parameters->out_file, "-")==0){
        parameters->out=stdout;
        parameters->out1=nullptr;
    } else {
        parameters->out=fopen(parameters->out_file, "w");
        parameters->out1=stderr;
    }
    if (parameters->out==nullptr){
        fprintf(stderr, "#Error, can not open the bam file\n");
        exit(0);
    }
}
int main(int argc, char *argv[]){
    auto parameters=new struct parameter_info;
    parse_args(argc, argv, parameters);
    BGZF *fp;
    bam_hdr_t *header;
    bam1_t *bam;
    fp=bgzf_open(parameters->bam_file, "r");
    header=bam_hdr_read(fp);
    bam=bam_init1();
    int clip_len;
    int index;
    int count=0;
    char seq[2000];
    unordered_map<uint64_t, struct read_sequence*> keeps;
    vector<read_sequence*> seqs;
    while(bam_read1(fp, bam)>0){
        if (bam->core.flag & (uint64_t) BAM_FUNMAP) continue;
        if ((clip_len=get_clip_len(bam)) < seed_len) continue;
        if (parameters->bin_size >= 0) index=((uint64_t)(bam->core.tid)<<32u)+(((uint64_t)(bam->core.pos)>>parameters->bin_size)<<2u)+
                    (((bam->core.flag & (uint64_t) BAM_FREAD1)>0))+
                    (((bam->core.flag & (uint64_t) BAM_FREAD2)>0)<<1u);
        else index=count;
        if (keeps.find(index)!=keeps.end()){
            auto read_seq=keeps[index];
            if (strlen(read_seq->clip_seq) < clip_len){
                delete [] read_seq->seq;
                read_seq->seq=newstr(get_seq(bam, seq, 2000));
                read_seq->clip_seq=read_seq->seq+strlen(read_seq->seq)-clip_len;
                read_seq->seed_seq=read_seq->clip_seq;
                read_seq->seq_end=read_seq->clip_seq+strlen(read_seq->clip_seq);

            }
        } else {
            auto read_seq = (keeps[index]=new struct read_sequence);
            read_seq->seq=newstr(get_seq(bam, seq, 2000));
            read_seq->clip_seq=read_seq->seq+strlen(read_seq->seq)-clip_len;
            read_seq->seed_seq=read_seq->clip_seq;
            read_seq->seq_end=read_seq->clip_seq+strlen(read_seq->clip_seq);
            ++count;
        }
    }
    if (fp) bgzf_close(fp);
    if (header) bam_hdr_destroy(header);
    if (bam) bam_destroy1(bam);
    for (auto it: keeps) if ((it.first & 3u)==0) seqs.push_back(it.second);
    if (!seqs.empty()) {
        if (parameters->out1) fprintf(parameters->out1, "#search adapters.\n");
        fprintf(parameters->out, "=== Adatpers ===\n");
        search_many_adapters(&seqs, parameters);
        seqs.clear();
    }
    for (auto it: keeps) if (it.first & 1u) seqs.push_back(it.second);
    if (!seqs.empty()) {
        if (parameters->out1) fprintf(parameters->out1, "#search adapters for read 1.\n");
        fprintf(parameters->out, "=== Adatpers for Read 1 ===\n");
        search_many_adapters(&seqs, parameters);
        seqs.clear();
    }
    for (auto it: keeps) if (it.first & 2u) seqs.push_back(it.second);
    if (!seqs.empty()) {
        if (parameters->out1) fprintf(parameters->out1, "#search adapters for read 2.\n");
        fprintf(parameters->out, "=== Adapters for Read 2 ===\n");
        search_many_adapters(&seqs, parameters);
        seqs.clear();
    }
    if(parameters->out1) fprintf(parameters->out1, "#finished!\n");
    fclose(parameters->out);
    delete parameters;
    return 0;
}
