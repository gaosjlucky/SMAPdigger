#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<vector>
#include<string>
#include<cstdio>
#include<set>
#include<algorithm>
#include<cmath>

#include "param.h"
#include "reads.h"
#include "dbseq.h"

using namespace std;

const unsigned int FIXSIZE=SEGLEN*FIXELEMENT; 
typedef gHit HitArray[MAXHITS+1];
typedef HitArray HitMatrix[MAXSNPS+1];

typedef bit32_t SegArray[FIXELEMENT*2];
typedef bit32_t SeedArray[16];

extern Param param;
extern char rev_char[];
extern char chain_flag[];

class SingleAlign {
public:
	SingleAlign();
	~SingleAlign();
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a);
	bit32_t CountNs();
	void set_RRBS_start(); //by yxi
	int TrimLowQual();
	void ConvertBinaySeq();
	inline void GenerateSeeds(int n, int start);
    int CountSeeds(RefSeq &ref, int n, bit32_t start);
	
	inline void CountMismatch(bit64_t *q, bit64_t *s);
	void ClearHits();
	int RunAlign(RefSeq &ref);
	int FilterReads();
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void Reverse_Seq();
	void s_OutHit(int chain, int n, bit8_t nspsn, gHit *hit, int insert_size, RefSeq &ref, string &os);
    
    //by yxi
    void SnpAlign(RefSeq &ref, bit32_t mode);
    int TrimAdapter();
    void SortHits4PE(int n);
    void Fix_Unpaired_Short_Fragment(RefSeq &ref);
    void ReorderSeed(RefSeq &ref);
    bit32_t GetTotalSeedLoc(RefSeq &ref, bit32_t start);
    void AdjustSeedStartArray(RefSeq &ref);
    inline bit32_t MismatchPattern0(bit64_t *q, bit64_t *s, bit32_t seglen_offset);
    inline bit32_t MismatchPattern1(bit64_t *q, bit64_t *s, bit32_t gap_index, bit32_t seglen_offset);
    bit32_t GapAlign(RefSeq &ref, bit32_t readlen, bit32_t mode, bit32_t seed_pos);
    //int MatchGap(bit32_t mmi1[], bit32_t mmi2[], bit32_t shift);    
    inline bit32_t AddHit(RefSeq &ref, bit32_t w, bit32_t mode);
    void Process_3nt();

public:
    //Hot data section
	//cache line 1
	bit32_t tmp_snp __attribute__((aligned(64))); //keep in the same cache line, 4 byte
	bit32_t snp_thres, map_readlen, raw_readlen; // 16 byte
	bit32_t read_chain_index, ref_chain_index, _seed, cseed_offset; //32 byte
	bit32_t rand_rSeed, n_aligned, read_max_snp_num, seedseg_num; // 48 byte
	bit32_t xflag_chain[2], xseed_start_offset[2]; // 64 byte

	//cache line 2
	gHit _hit, _ghit; // 16 byte
	bit32_t *_pro, *_refloc2; // 32 byte
	Hit *_refloc; //40 byte
    HitMatrix *xhits; // 48 byte
    HitArray *hits, *chits; //64 byte

    //cache line 3
    bit32_t N_count, tmp_snp0;

    SegArray xseq[2][SEGLEN] __attribute__((aligned(64))), *bseq, *cbseq;

	SeedArray xseeds[2][MAXSNPS+1];
	bit32_t xseed_array[2][FIXSIZE-SEGLEN];	
	bit32_t x_cur_n_hit[2][MAXSNPS+1], *_cur_n_hit, *_cur_n_chit;
	bit32_t xseed_start_array[2][MAXSNPS+1];

    set<ref_loc_t> *hitset, *ghitset; //, *chitset; 
    vector<pair<int,int> > xseedindex[2];
    //bit32_t read_chain_index; //forward or reverse reads, bseq or cbseq
    //bit32_t ref_chain_index; //Watson or Crick, ref.chr%2
	
	bit32_t mm_index[2*MAXGAPS+1][MAXSNPS+1];
	char cigar[16];
	bit64_t total_candidates, total_reads, total_seeds;

    //cold data section
	vector<ReadInf>::iterator _pread;
	string _outseq[2], _outqual[2];
	bit32_t num_reads;
	vector<ReadInf> mreads;
	string _str_align;   //align results, prepare for output
    vector<string>::iterator read_motif_iter;

	//local variables
	string::iterator _sp;
	string::reverse_iterator _sq;
	
	char _ch[1024];
    pair<ref_loc_t,bit32_t> seg_info;
    char _mapseq[256];
    string::iterator _readnt, _adapternt;
};

/*n=0: ab; 1: cd; 2: bc; 3: ac; 4: bd; 5: ad*/
//n<3: ab, cd, bc	
inline void SingleAlign::GenerateSeeds(int n, int start)
{
	bit32_t i;
 	//cout<<"cseed_offset="<<cseed_offset<<endl;
    if(param.RRBS_flag)	xseeds[read_chain_index][n][0]=xseed_array[read_chain_index][param.profile[n][0]+start+cseed_offset*read_chain_index];
    else for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) xseeds[read_chain_index][n][i]=xseed_array[read_chain_index][(*_pro)+start-i];
}

inline void SingleAlign::CountMismatch(register bit64_t *q, register bit64_t *s) {
/*
	if((tmp_snp=param.XM64((q[1]&param.XC64(s[1])^s[1])&q[6]))>snp_thres) return;
    if((tmp_snp+=param.XM64((q[0]&param.XC64(s[0])^s[0])&q[5]))>snp_thres) return;
    if((tmp_snp+=param.XM64((q[2]&param.XC64(s[2])^s[2])&q[7]))>snp_thres) return;
    tmp_snp+=param.XM64X2((q[3]&param.XC64(s[3])^s[3])&q[8],(q[4]&param.XC64(s[4])^s[4])&q[9]);
*/
	tmp_snp=N_count;
    if(param.nt3) {
    	for(bit32_t i=0;i<5;i++) {
    		if((tmp_snp+=param.XM64((param.XT64(q[i])^param.XT64(s[i]))&q[i+5]))>snp_thres) return;
    	}
    	tmp_snp0=N_count;
    	for(bit32_t i=0;i<5;i++) tmp_snp0+=param.XM64((q[i]&param.XC64(s[i])^s[i])&q[i+5]);
    }
    else {
    	for(bit32_t i=0;i<5;i++) {
    		if((tmp_snp+=param.XM64((q[i]&param.XC64(s[i])^s[i])&q[i+5]))>snp_thres) return;
		}
    }
}	

inline bit32_t SingleAlign::MismatchPattern0(bit64_t *q, bit64_t *s, bit32_t seglen_offset) {
    //assuming little endianness for x86, need to be re-written for portability
    bit64_t tmp; 
    bit32_t *mm_array; int i,ii, j, ss=0, end_element, end_offset, left, right, SEGLEN64=32, maxj;
    end_element=(map_readlen+seglen_offset-1)/SEGLEN64; left=SEGLEN64-seglen_offset; 
    right=(map_readlen+seglen_offset-1)%SEGLEN64+1; end_offset=SEGLEN64-right;
    //cout<<"in MMpattern()  gap_index:"<<gap_index<<endl;
    //cout<<"shift="<<(1-((int)gap_index%2)*2)*((int)gap_index+1)/2<<" offset="<<seglen_offset<<" end_element="<<end_element<<endl;
    //cout<<"end_offset="<<end_offset<<endl;
    mm_array=mm_index[0];
	tmp=param.swap_endian64(q[0]&param.XC64(s[0])^s[0])<<(seglen_offset*2);
	for(j=0;j<left;j++,tmp<<=2) if(tmp&0xC000000000000000ULL) {
		mm_array[ss++]=j;
		if(ss>(int)snp_thres-2) return j;
	}
	//cout<<"ss="<<ss<<" tmp_snp="<<tmp_snp<<" snp_thres="<<snp_thres<<endl;
	for(i=1,ii=left;i<=end_element;i++,ii+=SEGLEN64) {
		//disp_bfa(*s); 
		tmp=param.swap_endian64(q[i]&param.XC64(s[i])^s[i]);
		maxj=SEGLEN64-end_offset*(i==end_element);
        for(j=0;j<maxj;j++,tmp<<=2) if(tmp&0xC000000000000000ULL) {
			mm_array[ss++]=ii+j;
			if(ss>(int)snp_thres-2) return ii+j;
		}
		//cout<<"ss="<<ss<<" tmp_snp="<<tmp_snp<<" snp_thres="<<snp_thres<<endl;
    }
    for(;ss<=(int)snp_thres-2;ss++) mm_array[ss]=map_readlen;
	return map_readlen;
}

inline bit32_t SingleAlign::MismatchPattern1(bit64_t *q, bit64_t *s, bit32_t gap_index, bit32_t seglen_offset) {
    //assuming little endianness for x86, need to be re-written for portability
    bit64_t tmp;
    bit32_t *mm_array; int i,ii, j, ss=0, end_element, right, SEGLEN64=32, maxj, max0;
    end_element=(map_readlen+seglen_offset-1)/SEGLEN64; 
    right=(map_readlen+seglen_offset-1)%SEGLEN64+1; 
    //cout<<"in MMpattern()  gap_index:"<<gap_index<<endl;             
    //cout<<"shift="<<(1-((int)gap_index%2)*2)*((int)gap_index+1)/2<<" offset="<<seglen_offset<<" end_element="<<end_element<<endl;
    //cout<<"end_offset="<<end_offset<<endl;  
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(q[i])); cout<<endl;
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(s[i])); cout<<endl;
    mm_array=mm_index[gap_index];
	tmp=param.swap_endian64(q[end_element]&param.XC64(s[end_element])^s[end_element])>>((SEGLEN64-right)*2);
	max0=mm_index[0][snp_thres-2];
	if((param.XM64(tmp)>snp_thres-2)&&(right+max0+param.gap<map_readlen)) return 0;
	for(j=0;j<right;j++,tmp>>=2) if(tmp&0x3ULL) {
		mm_array[ss++]=j;
		if(ss>(int)snp_thres-2) return 1;
	}
	for(i=end_element-1,ii=right;i>=0;i--,ii+=SEGLEN64) {
        //disp_bfa(*s); 
		tmp=param.swap_endian64(q[i]&param.XC64(s[i])^s[i]);
		maxj=SEGLEN64-seglen_offset*(i==0);
		if((param.XM64(tmp)+ss>snp_thres-2)&&(ii+maxj+max0+param.gap<map_readlen)) return 0;
		for(j=0;j<maxj;j++,tmp>>=2) if(tmp&0x3ULL) {
			mm_array[ss++]=ii+j;
            if(ss>(int)snp_thres-2) return 1;
        }
	}
    for(;ss<=(int)snp_thres-2;ss++) mm_array[ss]=map_readlen;
    return map_readlen;
    /*
    for(i=0;i<=end_element;i++) disp_bfa(s[i]); cout<<endl;
    for(i=0;i<=end_element;i++) disp_bfa(q[i]); cout<<endl;
    for(i=0;i<=read_max_snp_num-2;i++) cout<<(int)mm_index[gap_index][i]<<"\t"; cout<<endl;
    */
    /*
    cout<<endl;
    for(i=0;i<FIXELEMENT;i++){
    	for(j=0;j<16;j++) cout<<(int)gap_pattern[shift_index][i*16+j];
    	cout<<" ";
    }
    cout<<endl;
    */
}	


inline bit32_t SingleAlign::AddHit(RefSeq &ref, bit32_t w, bit32_t mode) {
	/*
	if(_ghit.gap_size){
	bit32_t i; 
	cout<<"chain_index:"<<chain_index<<" w="<<w<<endl;
	for(i=0;i<=read_max_snp_num;i++) cout<<_cur_n_hit[i]<<"\t"; cout<<endl;
	for(i=0;i<=read_max_snp_num;i++) cout<<_cur_n_chit[i]<<"\t"; cout<<endl;
	}
	*/
	if((int)_ghit.loc<0) return 0; //underflow the start of refseq
  	if(_ghit.loc+map_readlen>ref.title[_ghit.chr].size) return 0; //overflow the end of refseq
    if(_ghit.gap_size) {
    	if(!ghitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0;
    }
    else {
    	if(!hitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0; //hit already exist
	}
	//cout<<"###HIT### "<<ref.title[_ghit.chr].name<<":"<<_ghit.loc<<" mis:"<<tmp_snp<<" snp_thres:"<<snp_thres<<endl;
	if(param.nt3) {
		if(tmp_snp0>31) tmp_snp0=31;
		_ghit.chr|=tmp_snp0<<11;
	}
    xhits[read_chain_index][w][x_cur_n_hit[read_chain_index][w]++]=_ghit;
    if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return 1;
    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) {
        if(w==0) return 1; 
        else snp_thres=w-1;
	}
    return 0;
}

inline void SingleAlign::Reverse_Seq() {
	_outseq[0]=_pread->seq; _outseq[1]=_pread->seq;
	reverse(_outseq[1].begin(), _outseq[1].end());
	for(string::iterator p=_outseq[1].begin(); p!=_outseq[1].end(); ++p) *p=rev_char[(unsigned char)*p];
	_outqual[0]=_pread->qual; _outqual[1]=_pread->qual;
	reverse(_outqual[1].begin(), _outqual[1].end());
}
#endif //_ALIGN_H_
