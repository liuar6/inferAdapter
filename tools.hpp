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

#include<unordered_map>
char * newstr(const char *str){
    char *cstr=new char[strlen(str)+1];
    strcpy(cstr, str);
    return cstr;
}
char ** split(char * line, char ** results, int length,char c='\t'){
    char *start=line;
    char *end=nullptr;
    int i=0;
    while ((end=strchr(start, c))!=nullptr && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (;i<length;++i) results[i]=nullptr;
    return results;
}
char* _complement(char *seq) {
    char *end=seq+strlen(seq);
    for (char *p=seq; p!=end; ++p) {
        switch (*p) {
            case 'A':  *p = 'T';break;
            case 'U':  *p = 'A';break;
            case 'T':  *p = 'A';break;
            case 'C':  *p = 'G';break;
            case 'G':  *p = 'C';break;
            case 'a':  *p = 't';break;
            case 'u':  *p = 'a';break;
            case 't':  *p = 'a';break;
            case 'c':  *p = 'g';break;
            case 'g':  *p = 'c';break;
        }
    }
    return seq;
}
char* _reverse(char *seq)
{
    int length=strlen(seq);
    int halfLen=length>>1;
    char *p=seq+length;
    char c;
    while (--halfLen>=0)
    {
        c=*seq;
        *seq++=*--p;
        *p=c;
    }
    return seq;
}
typedef struct faiItem {
    int lineCharNum;
    int lineBaseNum;
    long chromLen;
    long offset;
} faiItem;

class genomeReader{
private:
    ifstream genome;
    unordered_map <string, faiItem> genomeIndex;
public:
    ~genomeReader(){
        genome.close();
    }
    void loadGenome(const char* genomeFile){
        genome.open(genomeFile);
    }
    void loadGenomeIndex(const char* genomeIndexFile){
        char line[200];
        char *items[5];
        ifstream input(genomeIndexFile);
        input.getline(line, 200);
        while (!input.eof()){
            split(line, items, 5);
            genomeIndex[items[0]].chromLen=atol(items[1]);
            genomeIndex[items[0]].offset=atol(items[2]);
            genomeIndex[items[0]].lineBaseNum=atoi(items[3]);
            genomeIndex[items[0]].lineCharNum=atoi(items[4]);
            input.getline(line, 200);
        }
        input.close();
    }
    int load(const char* genomeFile, const char* genomeIndexFile){
        loadGenome(genomeFile);
        loadGenomeIndex(genomeIndexFile);
        return 1;
    }
    char* getSequence(const char *chrom, int chromStart,int chromEnd, char strand, bool toUpper=false){
        auto const &index=genomeIndex[chrom];
        int sepNum=index.lineCharNum-index.lineBaseNum;
        auto seq=new char[chromEnd-chromStart+sepNum];
        auto p=seq;
        int startBase=chromStart%index.lineBaseNum;
        int startLine=chromStart/index.lineBaseNum;
        int endBase=(chromEnd-1)%index.lineBaseNum;
        int endLine=(chromEnd-1)/index.lineBaseNum;
        long offset=index.offset+startLine*index.lineCharNum+startBase;
        int lineNum=endLine-startLine-1;
        genome.seekg(offset, ios::beg);
        if (genome.tellg()==-1) {cout<<"seekg error"<<endl; exit(1);}
        if (lineNum<0) genome.read(p, chromEnd-chromStart);
        else {
            genome.read(p, index.lineCharNum-startBase);
            p+=index.lineBaseNum-startBase;
            for (int i=0; i<lineNum; ++i) {
                genome.read(p, index.lineCharNum);
                p+=index.lineBaseNum;
            }
            genome.read(p, endBase+1);
        }
        seq[chromEnd-chromStart]='\0';
        if (strand=='-'){
            _complement(seq);
            _reverse(seq);
        }
        if (toUpper){
            char *pend=seq+chromEnd-chromStart;
            for (p=seq; p!=pend; ++p)
                *p=toupper(*p);
        }
        return seq;
    }
    int getLength(const char *chrom){
        return genomeIndex[chrom].chromLen;
    }
    bool exist(const char* chrom){
        return genomeIndex.find(chrom)!=genomeIndex.end();
    }
    unordered_map<string, faiItem>::iterator find(const char *chrom){
        return genomeIndex.find(chrom);
    }
    unordered_map<string, faiItem>::iterator begin(){
        return genomeIndex.begin();
    }
    unordered_map<string, faiItem>::iterator end(){
        return genomeIndex.end();
    }
};


#include <unordered_map>
#include <map>
#include <vector>
#include <string.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
using namespace std;

#define isUnmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define isReverse(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define getStrand(b) ((isReverse(b)?'-':'+'))
#define getName(b) ((char*)(b)->data)
#define getReferenceID(b) ((b)->core.tid)
#define getReferenceName(b, h) (h->target_name[(b)->core.tid])
#define getPosition(b) ((b)->core.pos)

#define getCigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define getCigarNum(b) ((b)->core.n_cigar)
#define getCigarType(b, i) (BAM_CIGAR_STR[((uint32_t*)((b)->data + (b)->core.l_qname))[i] & 0xF])
#define getCigarValue(b, i) ((((uint32_t*)((b)->data + (b)->core.l_qname))[i] & 0xFFFFFFF0)>>4)

#define getQueryLength(b) ((b)->core.l_qseq)
#define getQuery(b) ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define getQuality(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define getBase(s, i) (seq_nt16_str[(s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf])


class bamReader{
public:
    BGZF *fp;
    bam_hdr_t *header;
    bam1_t *b;
    unordered_map<string, int> chrom2id;
    unordered_map<int, uint64_t> offset;

    explicit bamReader(const char *fn){
        fp=nullptr;
        header=nullptr;
        b=nullptr;
        open(fn);
    }
    explicit bamReader(){
        fp=nullptr;
        header=nullptr;
        b=nullptr;
    }
    int open(const char *fn){
        fp=bgzf_open(fn, "r");
        header=bam_hdr_read(fp);
        b=bam_init1();
        for (int i=0; i<header->n_targets; ++i) chrom2id[header->target_name[i]]=i;
        return 1;
    }
    int reopen(){
        if (header) bam_hdr_destroy(header);
        if (b) bam_destroy1(b);
        bgzf_seek(fp, 0, SEEK_SET);
        header=bam_hdr_read(fp);
        b=bam_init1();
        return 1;
    }
    int jumpToChrom(const char* chrom){
        return bgzf_seek(fp, offset[chrom2id[chrom]], SEEK_SET);
    }
    bam1_t* next(){
        if (b==nullptr){
            b=bam_init1();
        }
        if (bam_read1(fp, b)>0) return b;
        else return nullptr;
    }
    static int return_with_error(const char* message, int ret){
        if (ret){
            cout<<message<<endl;
            return ret;
        }
    }
    int parse(){
        reopen();
        int last_coor=0;
        int last_tid=-1;
        u_int64_t last_offset=bgzf_tell(fp);
        int tid;
        bool no_coor=false;
        while (next()){
            tid=b->core.tid;
            if (isUnmapped(b)) no_coor=true;
            if (tid>=0 && last_tid!=tid){
                if (no_coor) return_with_error("NO_COOR reads not in a single block at the end", -1);
                if (offset.find(tid)!=offset.end()) return_with_error("Chromosome blocks not continuous", -1);
                else offset[tid]=last_offset;
            }
            else if (tid >= 0 && last_coor > (b->core.pos)) return_with_error("test", 0);
            last_tid=tid;
            last_offset=bgzf_tell(fp);
            last_coor=b->core.pos;
        }
        reopen();
        return 1;
    }
    ~bamReader(){
        if (fp) bgzf_close(fp);
        if (header) bam_hdr_destroy(header);
        if (b) bam_destroy1(b);
    }
};