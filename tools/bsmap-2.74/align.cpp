#include "align.h"


extern bit8_t alphabet[];
extern bit8_t reg_alphabet[];
extern bit8_t rev_alphabet[];
extern char rev_char[];
extern char nt_code[];
extern char revnt_code[];

//create seed profile
SingleAlign::SingleAlign() {
	n_aligned=0;
    _cur_n_hit=x_cur_n_hit[0];
    _cur_n_chit=x_cur_n_hit[1];

    bseq=xseq[0]; cbseq=xseq[1];

    xhits = new HitMatrix[2];
    hits = xhits[0]; chits = xhits[1];
	
	hitset= new set<ref_loc_t>[param.total_ref_seq];
	ghitset = new set<ref_loc_t>[param.total_ref_seq];

	total_candidates=total_reads=total_seeds=0;	
    if(!param.pairend) _str_align.reserve(BatchNum*400);
	//chitset= new set<ref_loc_t>[param.total_ref_seq];
	rand_rSeed=getpid()*time(NULL);
	N_count=0;
	//for(int i=0;i<256;i++) nothing[i]=0;
	//cout<<"cache line 1:"<<&tmp_snp<<endl;
	//cout<<"cache line 2:"<<&_hit<<endl;
}

SingleAlign::~SingleAlign() {
	delete [] hitset;
	delete [] ghitset;
    delete [] xhits;
    //for(int i=0;i<60;i++) cout<<i<<"\t"<<nothing[i+100]<<"\t"<<nothing[i]<<"\t"<<(nothing[i+100]/(nothing[i]+1))<<endl; cout<<endl;
	//delete [] chitset;

    //cerr<<"Total candidates:"<<total_candidates<<" \ttotal_seeds:"<<total_seeds<<" \ttotal_reads:"<<total_reads<<endl;
    //cerr<<total_candidates/total_reads<<"candidates/reads\n";
    //cerr<<total_seeds/total_reads<<"seeds/reads\n";
    //cerr<<total_candidates/total_seeds<<"candidates/seeds\n";
}

void SingleAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a) {
	num_reads=n;
	mreads=a;
}

bit32_t SingleAlign::CountNs() {
	int n=0;
	for(_sp=_pread->seq.begin(); _sp!=_pread->seq.end(); _sp++) {
		if(!reg_alphabet[(unsigned char)*_sp]) n++;
    }
    if(param.N_mis) N_count=n; 
	return n;	
}


//trim low quality at 3'-end, cut at 3bp continuous high-quality bps
int SingleAlign::TrimLowQual() { //return 0: normal trim, 1: read too short after trimming
	bit32_t i; bit8_t qual_thres=param.zero_qual+param.qual_threshold;

    if(_pread->seq.size()!=_pread->qual.size())
        _pread->qual=string(_pread->seq.size(), param.zero_qual+param.default_qual);

    qual_thres=param.zero_qual+param.qual_threshold;
    if(param.out_sam&&param.zero_qual!='!'){
        for(_sp=_pread->qual.begin(); _sp!=_pread->qual.end();++_sp) {
            //cout<<*_sp<<" "<<(int)*_sp<<" ";
            *_sp-=(param.zero_qual-'!');
            //cout<<(int)*_sp<<" \t";
        }
        qual_thres-=(param.zero_qual-'!');
    }

    if(param.qual_threshold==0) return 0;
    for(_sq=_pread->qual.rbegin(),i=_pread->qual.size(); _sq!=_pread->qual.rend(); i--,++_sq) {
        //cout<<"i="<<i<<" *_sq:"<<*_sq<<' '<<(int) *_sq<<"\t\t";
        if(*_sq>qual_thres) break;
    }
    
	if(i<param.seed_size+param.index_interval-1) return 1;
	_pread->qual.erase(i); _pread->seq.erase(i);
	return 0;
}

//shift 1bp to left
inline void RightShiftBinSeq(bit32_t *orib, bit32_t *newb) {
/*
	newb[0]=orib[0]>>2; newb[1]=orib[1]>>2|orib[0]<<(SEGLEN*2-2);
	newb[4]=orib[4]>>2|orib[1]<<(SEGLEN*2-2); newb[5]=orib[5]>>2|orib[4]<<(SEGLEN*2-2); 
	newb[8]=orib[8]>>2|orib[5]<<(SEGLEN*2-2); newb[9]=orib[9]>>2|orib[8]<<(SEGLEN*2-2);
	newb[12]=orib[12]>>2|orib[9]<<(SEGLEN*2-2); newb[13]=orib[13]>>2|orib[12]<<(SEGLEN*2-2);
	newb[16]=orib[16]>>2|orib[13]<<(SEGLEN*2-2); newb[17]=orib[17]>>2|orib[16]<<(SEGLEN*2-2);
*/
	bit32_t i; newb[0]=orib[0]>>2;
    for(i=1;i<FIXELEMENT;i++) newb[i]=(orib[i]>>2)|(orib[i-1]<<(SEGLEN*2-2));
}

//convert string seq to binary type
void SingleAlign::ConvertBinaySeq() {
	bit32_t i, h, _a, _b, s;
	Reverse_Seq();
    xflag_chain[0]=param.chains||(_pread->readset<2);
    xflag_chain[1]=param.chains||(_pread->readset==2);
    if(xflag_chain[0]) {	//direct chain
    	h=_a=_b=s=0;
    	for(_sp=_pread->seq.begin(),i=1; _sp!=_pread->seq.end(); _sp++,i++) {
    		_a<<=2; _b<<=2;
    		_a|=alphabet[(unsigned char)*_sp];
    		_b|=reg_alphabet[(unsigned char)*_sp];
    		if(i>param.seed_size){
    		    s<<=2; s|=_a&0x3; 
    		    xseed_array[0][i-param.seed_size]=param.XT(s&param.seed_bits);
    		}
    		else if(i==param.seed_size){s=_a; xseed_array[0][0]=param.XT(s);}
    		if(0==i%SEGLEN) {
    			bseq[0][h]=_a;
    			bseq[0][h+FIXELEMENT]=_b;
    			_a=_b=0; h++;
    		}		
    	}
    	for(; i!=FIXSIZE+1; i++) {
      		_a<<=2; _b<<=2;
    		if(0==i%SEGLEN) {
    			bseq[0][h]=_a;
    			bseq[0][h+FIXELEMENT]=_b;
    			_a=_b=0; h++;
    		}		
    	}
    	for(i=1; i!=SEGLEN; i++) {
    		RightShiftBinSeq(bseq[i-1], bseq[i]);
    		RightShiftBinSeq(bseq[i-1]+FIXELEMENT, bseq[i]+FIXELEMENT);
    	}
	}
/*
	int ii;
	for(i=0; i!=SEGLEN; i++) {
		cout<<"bin seq: "<<i<<endl;
		for(ii=0;ii<FIXELEMENT;ii++) disp_bfa(bseq[i][ii]); cout<<endl;
		for(;ii<FIXELEMENT*2;ii++) disp_bfa(bseq[i][ii]); cout<<endl;
	}
*/	
    if(xflag_chain[1]) { //reverse seq
    	h=_a=_b=s=0;
    	for(_sq=_pread->seq.rbegin(),i=1; _sq!=_pread->seq.rend(); _sq++,i++) {
    		_a<<=2; _b<<=2;
    		_a|=rev_alphabet[(unsigned char)*_sq];
    		_b|=reg_alphabet[(unsigned char)*_sq];
    		if(i>param.seed_size){
    		    s<<=2; s|=_a&0x3; 
    		    xseed_array[1][i-param.seed_size]=param.XT(s&param.seed_bits);
    		}
    		else if(i==param.seed_size){s=_a; xseed_array[1][0]=param.XT(s);}
    		if(0==i%SEGLEN) {
    			cbseq[0][h]=_a;
    			cbseq[0][h+FIXELEMENT]=_b;
    			_a=_b=0; h++;
    		}		
    	}
    	for(; i!=FIXSIZE+1; i++) {
    		_a<<=2; _b<<=2;
    		if(0==i%SEGLEN) {
    			cbseq[0][h]=_a;
    			cbseq[0][h+FIXELEMENT]=_b;
    			_a=_b=0; h++;
    		}		
    	}
    	
    	for(i=1; i!=SEGLEN; i++) {
    		RightShiftBinSeq(cbseq[i-1], cbseq[i]);
    		RightShiftBinSeq(cbseq[i-1]+FIXELEMENT, cbseq[i]+FIXELEMENT);
    	}
	}			
}


void SingleAlign::SnpAlign(RefSeq &ref, bit32_t mode) {
   	bit32_t i,j,m, modeindex, cmodeindex, mc, h, read_chain_index_mask, _z;
   	Hit prefetch_hit, _hHit;

    if(param.RRBS_flag){ //RRBS mode
    	for(read_chain_index=0; read_chain_index<2;read_chain_index++){
            if(!xflag_chain[read_chain_index]) continue;
            cmodeindex=modeindex=xseedindex[read_chain_index][mode].second;
            if(read_chain_index) cmodeindex=map_readlen/param.seed_size-1-modeindex;
       		_seed=xseeds[read_chain_index][modeindex][0];
    		//cout<<"mode:"<<mode<<" m="<<ref.index[_seed].n1<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<' '<<_seed<<endl;
    		if((m=ref.index[_seed].n1)==0) continue; 
 	    	_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
       		h=param.profile[modeindex][0]+cseed_offset*read_chain_index;
       		read_chain_index_mask=read_chain_index<<24;
       		for(j=0; j!=m; j++) {
   	    		_hHit=_refloc[j];
   	    		if(((_hHit.chr^read_chain_index_mask)>>16)!=cmodeindex) continue; // mode or strand not match
                if(j+PREFETCH_LOOP<m){      
                    prefetch_hit=_refloc[j+PREFETCH_LOOP];
                    //cache_ptr=ref.bfa[prefetch_hit.chr&0xffff].s+(prefetch_hit.loc-h)/SEGLEN;                            
                    //__builtin_prefetch(cache_ptr,0,0); __builtin_prefetch(cache_ptr+1,0,0);
                    __builtin_prefetch(ref.bfa[prefetch_hit.chr&0xffff].s+(prefetch_hit.loc-h)/SEGLEN,0,0);
                }                                                                                        
   	    		_hit.chr=_hHit.chr&0xffff;
   	    		if(_hHit.loc<h) continue; //underflow the start of refseq
   	    		_hit.loc=_hHit.loc-h; _z=_hit.loc%SEGLEN;
                CountMismatch((bit64_t*)xseq[read_chain_index][_z], (bit64_t*)(ref.bfa[_hit.chr].s+_hit.loc/SEGLEN));
   	    		//cout<<" j="<<j<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<endl;
   	    	 	if(tmp_snp<=snp_thres) {
                    _ghit=_hit;
       	    	 	if(_ghit.chr&1) _ghit.loc=ref.title[_ghit.chr].rc_offset-map_readlen-_ghit.loc;
                    if((!param.pairend)&&(read_chain_index==0)){
                        seg_info=ref.CCGG_seglen(_ghit.chr, _ghit.loc, map_readlen); //get fragment information
                        //cout<<"seg1:"<<seg_info.first<<" seg2:"<<seg_info.second<<endl;
                        if(seg_info.second>param.max_insert) continue; // fragent too large
                        if(seg_info.second<param.min_insert) continue; // fragment too small
                    }
                    _ghit.gap_size=0;
                    if(param.nt3) _ghit.chr|=tmp_snp0<<11;                    
          			if(AddHit(ref, tmp_snp, mode)) return;
                }
                if(param.gap>0) {
                    _hit.loc+=ref.ref_anchor[_hit.chr/2];
                    if(GapAlign(ref, map_readlen*ref_chain_index, mode, h)) return;
                }                             
    		}
    	}
    }
    else{ //WGBS mode
        for(read_chain_index=0; read_chain_index<2;read_chain_index++){
            if(!xflag_chain[read_chain_index]) continue;
    	    modeindex=xseedindex[read_chain_index][mode].second;
    	    //cout<<"read_chain:"<<read_chain_index<<" modeindex:"<<modeindex<<" mode:"<<mode<<endl;
        	for(i=0; i!=param.index_interval; i++) {
        		_seed=xseeds[read_chain_index][modeindex][i];
        		m=ref.index2[_seed].n[0]; 
        		//cout<<" i="<<i<<" m="<<m<<" mc="<<mc<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<endl;
        		total_seeds++;
				if(m==0||m>param.max_kmer_num) continue;
				total_candidates+=m;
				mc=ref.index2[_seed].n[1]-1;
        		h=param.profile[modeindex][i]+xseed_start_array[read_chain_index][modeindex]-i;
        		//cout<<" h="<<h<<" "<<xseed_start_array[read_chain_index][modeindex]<<" "<<endl;;
        		for(j=0,_refloc2=ref.index2[_seed].loc1; j!=m; j++,_refloc2++) { 
					//cout<<"j="<<j<<"prefetch:"<<_refloc2+PREFETCH_LOOP<<endl;
					//cout<<"content:"<<*(_refloc2+PREFETCH_LOOP)<<endl;
					//if(j+PREFETCH_LOOP<m) {
						__builtin_prefetch(ref.xref[(bit32_t)((int)mc-(int)j-(int)PREFETCH_LOOP)>>31]+(*(_refloc2+PREFETCH_LOOP)-h)/SEGLEN,0,0);
						//_c=(*(_refloc2+PREFETCH_LOOP)-h)/SEGLEN;
						//cache_ptr=ref.xref[(j+PREFETCH_LOOP>=mc)]+(*(_refloc2+PREFETCH_LOOP)-h)/SEGLEN;
						//__builtin_prefetch(cache_ptr,0,0); __builtin_prefetch(cache_ptr+1,0,0);
					//}
        			_hit.loc=(*_refloc2)-h; 
        			ref_chain_index=(bit32_t)((int)mc-(int)j)>>31; 
        			//cerr<<"j="<<j<<" mc="<<mc<<" chain="<<ref_chain_index<<endl;
                    //_ghit=ref.int2hit(_hit,map_readlen*ref_chain_index,0,0);
					//cout<<" j="<<j<<" ref_chain:"<<ref_chain_index<<" refloc:"<<(*_refloc2)<<" h:"<<h<<" loc:"<<ref.int2hit(_hit,map_readlen*ref_chain_index,0,0).loc<<endl;
                    CountMismatch((bit64_t*)xseq[read_chain_index][_hit.loc%SEGLEN], (bit64_t*)(ref.xref[ref_chain_index]+(_hit.loc/SEGLEN)));
					//cout<<"id:"<<_hit.loc%SEGLEN<<" mis:"<<w<<" snp_thres:"<<snp_thres<<endl;
                    if(tmp_snp<=snp_thres) {
                        _ghit=ref.int2hit(_hit,map_readlen*ref_chain_index,0,0);
                        //cout<<"###HIT### "<<ref.title[_ghit.chr].name<<":"<<_ghit.loc<<" mis:"<<tmp_snp<<" snp_thres:"<<snp_thres<<endl;                    
                        if(AddHit(ref, tmp_snp, mode)) return;
                    }
                    if(param.gap>0) 
                        if(GapAlign(ref, map_readlen*ref_chain_index, mode, h)) return;
        		}
        	}
    	}
    }
}

bit32_t SingleAlign::GapAlign(RefSeq &ref, bit32_t readlen, bit32_t mode, bit32_t seed_pos) { 
    //chain_index: read sequence direction, 0: original seq (forward), 1: RC seq (reverse strand)
    //readlen: reference sequence direction, 0: Watson strand, >0: Crick strand (value=read length)
    bit32_t ghit_loc=_hit.loc, t=0, tt, ghit_loc1; int shift, shift1, clip;
    //gHit hh=_hit; ref.int2hit(hh,readlen); cout<<ref.title[hh.chr].name.c_str()<<":"<<hh.loc<<endl;
    SegArray *readseq; bit32_t *refseq, i, j, gap_snp, gap_pos, m2, rl; bit32_t *mmi1, *mmi2; 
    if(snp_thres<2) return 0;
	//for(i=0;i<=snp_thres-2;i++) cout<<" i="<<i<<" mmi0="<<mm_index[0][i]; cout<<endl; 
    refseq=ref.xref[readlen!=0]; readseq=xseq[read_chain_index];
    //cout<<"ghit_loc="<<ghit_loc<<" "<<(int)ghit_loc<<" readlen="<<readlen<<" refseq="<<refseq<<endl;
    if(MismatchPattern0((bit64_t*)readseq[ghit_loc%SEGLEN], (bit64_t*)(refseq+(ghit_loc/SEGLEN)), ghit_loc%SEGLEN)<seed_pos+param.seed_size) return 0;
    mmi1=mm_index[0]; 
    for(tt=1;tt<=param.gap*2;tt++) {
        //cout<<"tt="<<tt<<" t="<<t<<endl;
        t=(tt+1)/2; shift=(1-(tt%2)*2)*t; shift1=shift*(shift<0);
        if(snp_thres<1+t) break;
        ghit_loc1=ghit_loc+shift; 
        if(MismatchPattern1((bit64_t*)readseq[(ghit_loc1)%SEGLEN], (bit64_t*)(refseq+((ghit_loc1)/SEGLEN)), tt, (ghit_loc1)%SEGLEN)==0) continue;
		//cout<<"shift="<<shift; for(i=0;i<=snp_thres-2;i++) cout<<" i="<<i<<" mmi0="<<mm_index[tt][i]; cout<<endl;
        rl=map_readlen-t-1; mmi2=mm_index[tt];
        for(i=0;i<snp_thres-t;i++) {
            gap_pos = mmi1[i]; gap_snp=0;
            if(gap_pos<param.gap_edge||gap_pos>=rl) continue;
            for(j=0;j<snp_thres-t-i;j++) {
                m2=mmi2[j];
                if(m2<param.gap_edge||m2>=rl) continue;
                if((int)gap_pos+(int)m2-shift1<(int)map_readlen) continue;
                gap_snp=i+j+1+t;
                clip=(int)gap_pos+(int)param.gap_edge-(int)map_readlen;
                clip-=shift1;
                if(clip>0) gap_pos-=clip;
                //cout<<"clip="<<clip<<" gap_pos="<<gap_pos<<" gap_snp="<<gap_snp<<endl;
                _ghit=ref.int2hit(_hit,readlen,shift,gap_pos);
                //cout<<"i1="<<i<<" m1="<<gap_pos<<" i2="<<j<<" m2="<<m2<<" shift="<<shift<<" gap_snp="<<gap_snp<<" clip="<<clip<<" hitloc="<<_ghit.loc<<endl;
                //cout<<"gap="<<(int)_hit.gap_size<<" gap_pos="<<(int)_hit.gap_pos<<" gap_snp="<<gap_snp<<" "<<ref.title[_hit.chr].name.c_str()<<":"<<_hit.loc<<endl;
                return AddHit(ref, gap_snp, mode);
            }
            //if(gap_snp) break;
        }
        /*         
        if((gg=MatchGap(mm_index[0], mm_index[tt], tt))>0) {
            _ghit=ref.int2hit(_hit,readlen,shift,gg);
            //cout<<"gap="<<(int)_hit.gap_size<<" gap_pos="<<(int)_hit.gap_pos<<" gap_snp="<<gap_snp<<" "<<ref.title[_hit.chr].name.c_str()<<":"<<_hit.loc<<endl;
            //cout<<endl<<_pread->name.c_str()<<"#"<<_pread->readset<<" \t"<<_pread->seq<<endl;
            if(AddHit(ref, chain_index, gap_snp, mode)) return 1; 
            //break;
        }
        */
    }
    return 0;

    /*
    bit32_t i, rev_flag;
    if(gg>0&&t==param.gap&&gap_snp<=snp_thres&&readlen>0){
        cout<<"read_dir="<<chain_index<<" ref_dir="<<(int)(readlen!=0)<<endl;
        cout<<"gap="<<(int)_hit.gap_size<<" gap_pos="<<(int)_hit.gap_pos<<" gap_snp="<<gap_snp<<" "<<ref.title[_hit.chr].name.c_str()<<":"<<_hit.loc<<endl;
        //for(i=0;i<=read_max_snp_num-2;i++) cout<<(int)mm_index[0][i]<<"\t"; cout<<endl;
        rev_flag=chain_index^(readlen!=0); gg=(bit32_t) _hit.gap_pos;
        if(tt%2) { //insert on reference
            for(i=0;i<gg;i++) cout<<param.useful_nt[(ref.bfa[_hit.chr&0xfffeU].s[(_hit.loc+i)/SEGLEN]>>((SEGLEN-1-(_hit.loc+i)%SEGLEN)*2))&0x3];
            for(i=0;i<t;i++) cout<<"-";
            for(i=gg;i<map_readlen-t;i++) cout<<param.useful_nt[(ref.bfa[_hit.chr&0xfffeU].s[(_hit.loc+i)/SEGLEN]>>((SEGLEN-1-(_hit.loc+i)%SEGLEN)*2))&0x3];
            cout<<endl;
            for(i=0;i<map_readlen;i++) cout<<_outseq[rev_flag][i]; cout<<endl;
        }
        else { //insert on read
            for(i=0;i<map_readlen+t;i++) cout<<param.useful_nt[(ref.bfa[_hit.chr&0xfffeU].s[(_hit.loc+i)/SEGLEN]>>((SEGLEN-1-(_hit.loc+i)%SEGLEN)*2))&0x3];
            cout<<endl;
            for(i=0;i<gg;i++) cout<<_outseq[rev_flag][i];
            for(i=0;i<t;i++) cout<<"-";
            for(i=gg;i<map_readlen;i++) cout<<_outseq[rev_flag][i]; cout<<endl;
        }
        //for(i=0;i<=read_max_snp_num-2;i++) cout<<(int)mm_index[tt][i]<<"\t"; cout<<endl;  
        if(++nothing==10) exit(0);
    }
    return 0;
    */
}                                            

void SingleAlign::SortHits4PE(int n) {
    //cout<<"SORT "<<_pread->name.c_str()<<" "<<n<<"\t"<<_cur_n_hit[n]<<" "<<_cur_n_chit[n]<<endl;
	sort(hits[n], hits[n]+_cur_n_hit[n], HitComp);
	sort(chits[n], chits[n]+_cur_n_chit[n], HitComp);
}


int SingleAlign::TrimAdapter() {
    bit32_t i,k,m,m0, pos;
    raw_readlen=_pread->seq.size();
    if(param.RRBS_flag){
        for(i=0; i<param.n_adapter; i++){
            for(pos=param.seed_size+param.index_interval-1;pos<_pread->seq.size()-5;pos++){
                m0=0;
                for(k=0,_readnt=_pread->seq.begin()+pos,_adapternt=param.adapter[i].begin();k<param.adapter[i].size()&&k<15&&_readnt!=_pread->seq.end();k++,++_readnt,++_adapternt){
    		        if((m0+=(*_adapternt!=*_readnt))>4) break;
                }
                if(k<m0*5) continue;
    	        m=m0;
                for(_adapternt=param.digest_site.begin(),_readnt=_pread->seq.begin()+pos-param.digest_site.size()+param.digest_pos;
                    _adapternt!=param.digest_site.end()-param.digest_pos;++_adapternt,++_readnt) {
    		        m+=(*_adapternt!=*_readnt)&&(*_adapternt!='C'||*_readnt!='T');
    	        }
                if(k>=m*5) {
                    _pread->seq.erase(pos); 
                    if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                    return 1;
                }        
                if(param.pairend) { //pair-end
            		m=m0;
                    for(_adapternt=param.digest_site.begin(),_readnt=_pread->seq.begin()+pos-param.digest_site.size()+param.digest_pos;
                        _adapternt!=param.digest_site.end()-param.digest_pos;++_adapternt,++_readnt) {
    	        	    m+=(*_adapternt!=*_readnt)&&(*_adapternt!='G'||*_readnt!='A');
    		        }
                    if(k>=m*5) {
                        _pread->seq.erase(pos); 
                        if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                        return 1;
                    }
                }
            }
        }
    }
    else{
        for(i=0; i<param.n_adapter; i++){
            for(pos=param.seed_size+param.index_interval-1;pos<_pread->seq.size()-4;pos++){
                m0=0;
                for(k=0,_readnt=_pread->seq.begin()+pos,_adapternt=param.adapter[i].begin();k<param.adapter[i].size()&&k<15&&_readnt!=_pread->seq.end();k++,++_readnt,++_adapternt){
                    if((m0+=(*_adapternt!=*_readnt))>4) break;
                }
                if(k>=m0*5&&k>3) {
                    _pread->seq.erase(pos); 
                    if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                    return 1;
                }    
            }
        }    
    }
    return 0;
}


void SingleAlign::ClearHits() {
    bit32_t i;
	for(i=0; i<=MAXSNPS; i++) _cur_n_hit[i]=_cur_n_chit[i]=0;
	for(i=0; i<param.total_ref_seq; i++) {hitset[i].clear(); ghitset[i].clear();} 
	map_readlen=_pread->seq.size();
}

int SingleAlign::RunAlign(RefSeq &ref) {
    bit32_t i, ii;
    //cout <<_pread->name.c_str()<< "\t" << _pread->seq.c_str()<< endl;
	ClearHits();
    seedseg_num=min((int)((map_readlen-param.index_interval+1)/param.seed_size),(int)(read_max_snp_num+1));
	ConvertBinaySeq(); 
    snp_thres=read_max_snp_num;
    cseed_offset=_pread->seq.size()%param.seed_size;
    ReorderSeed(ref);
    for(i=0; i<seedseg_num; i++){
        //GenerateSeeds_1(i);
        SnpAlign(ref,i);  
        if(!param.RRBS_flag) for(ii=0;ii<=i;ii++) if(_cur_n_hit[ii]||_cur_n_chit[ii]) return 1;
	}
    for(i=0; i<=read_max_snp_num; i++) if(_cur_n_hit[i]||_cur_n_chit[i]) return 1;         
    return 0;
}

void SingleAlign::ReorderSeed(RefSeq &ref) {
    bit32_t i,ii,s, total, tt;
	//for(i=0;i<map_readlen-param.seed_size;i++) {disp_bfa(param.map3to4(seed_array[i])); cout<<i<<":"<<ref.index2[seed_array[i]].n[0]<<":"<<ref.index2[seed_array[i]].n[1]<<endl;}
    for(read_chain_index=0;read_chain_index<2;read_chain_index++) {
        if(!xflag_chain[read_chain_index]) continue; 
        if(param.RRBS_flag) xseed_start_offset[read_chain_index]=0;
        else {
            total=0xffffffffu;
            ii=(map_readlen-param.index_interval+1)%param.seed_size;
            for(i=0;i<ii;i++) {
                tt=GetTotalSeedLoc(ref, i);
                if(tt<total) {total=tt; xseed_start_offset[read_chain_index]=i;}
            }
        }
        //cout<<"start="<<xseed_start_offset[0]<<" cstart="<<xseed_start_offset[1]<<endl;
        AdjustSeedStartArray(ref);
        xseedindex[read_chain_index].clear();
        for(i=0; i<seedseg_num; i++){
            GenerateSeeds(i, xseed_start_array[read_chain_index][i]);  
            s=0;
            if(param.RRBS_flag) s+=ref.index[xseeds[read_chain_index][i][0]].n1;
            else for(ii=0; ii!=param.index_interval; ii++) s+=ref.index2[xseeds[read_chain_index][i][ii]].n[0];
            //if(s==0) s=0x7FFFFFFFL;
            xseedindex[read_chain_index].push_back(pair<int,int>(s,i)); 
        }
        //for(vector<pair<int,int> >::iterator sit=xseedindex[read_chain_index].begin();sit!=xseedindex[read_chain_index].end();++sit) cout <<sit->first<<" "<<sit->second<<"   "; cout<<endl;
        sort(xseedindex[read_chain_index].begin(), xseedindex[read_chain_index].end());
        //for(vector<pair<int,int> >::iterator sit=xseedindex[read_chain_index].begin();sit!=xseedindex[read_chain_index].end();++sit) cout <<sit->first<<" "<<sit->second<<"   "; cout<<endl;
    }
}

void SingleAlign::AdjustSeedStartArray(RefSeq &ref) {
    bit32_t i, ii, tt, total;
    bit32_t ptr, start, end, max_offset;
    //cout<<"max offset:"<<max_offset<<" seedseg_num="<<seedseg_num<<" seed_start_offset="<<seed_start_offset<<endl;
    for(i=0;i<seedseg_num;i++) xseed_start_array[read_chain_index][i]=xseed_start_offset[read_chain_index];
    if(param.RRBS_flag) return;
    max_offset=(map_readlen-param.index_interval+1)%param.seed_size; 
    for(i=0;i<seedseg_num;i++) {
        if(i%2==0) ptr=i/2; else ptr=seedseg_num-1-i/2;
        total=0xffffffff; 
        //cout<<"i="<<i<<" ptr="<<ptr<<endl;
        if(ptr==0) start=0; else start=xseed_start_array[read_chain_index][ptr-1];
        if(ptr==seedseg_num-1) end=max_offset; else end=xseed_start_array[read_chain_index][ptr+1];
        //cout<<"start="<<start<<" end="<<end<<endl;
        xseed_start_array[read_chain_index][ptr]=start;
        for(ii=start;ii<=end;ii++) {
            tt=CountSeeds(ref,ptr,ii);
            //cout<<"ii="<<ii<<" tt="<<tt<<endl;
            if(tt<total) {total=tt; xseed_start_array[read_chain_index][ptr]=ii;}
        }
        //for(ii=0;ii<seedseg_num;ii++) cout<<ii<<":"<<xseed_start_array[read_chain_index][ii]<<" "; cout<<endl;
    }
}    

int SingleAlign::CountSeeds(RefSeq &ref, int n, bit32_t start) {
	bit32_t i, total=0; bit32_t s;
    for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) {
        s=xseed_array[read_chain_index][*_pro+start-(int)i];
        total+=ref.index2[s].n[0];
    }
    return total;
}

bit32_t SingleAlign::GetTotalSeedLoc(RefSeq &ref, bit32_t start) {
    bit32_t i, total=0;
    for(i=0; i<seedseg_num; i++) total+=CountSeeds(ref, i, start);
    return total;
}

int SingleAlign::FilterReads() {
    //if(_pread->bam_flag&0x200) return 1;
    if(param.max_snp_num<100) read_max_snp_num=param.max_snp_num;
    else read_max_snp_num=(bit32_t)((param.max_snp_num-100)/100.0*_pread->seq.size()+0.5);
    if(read_max_snp_num>MAXSNPS) read_max_snp_num=MAXSNPS;
    TrimAdapter();
    if(TrimLowQual()!=0) return 1;
	if(_pread->seq.size()<param.min_read_size) return 1;
    if(CountNs()>param.max_ns) return 1;
    read_max_snp_num=(read_max_snp_num+1)*(_pread->seq.size()-1)/raw_readlen;
	return 0;
}

void SingleAlign::Do_Batch(RefSeq &ref) {
	_str_align.clear(); 
	bit32_t tt;
	//alignment for normal sequences
	for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
        //cout<<_pread->index<<" "<<_pread->name.c_str()<<endl<< _pread->seq.c_str()<<"\t"<<_pread->qual.c_str() << endl;
        total_reads++;
		if(FilterReads()) {
			if(param.report_repeat_hits) s_OutHit(0, -1, 0, hits[0], 0, ref, _str_align);
        }
        else {
    		RunAlign(ref);
        	StringAlign(ref, _str_align);
        }
	}
}

void SingleAlign::Process_3nt() {
    int ii, jj, mm;
    for(ii=read_max_snp_num; ii>=0;--ii) {
        //cout<<ii<<" "<<_cur_n_hit[ii]<<" "<<_cur_n_chit[ii]<<endl;
        for(jj=_cur_n_hit[ii]-1;jj>=0;--jj) {
            _ghit=hits[ii][jj]; mm=_ghit.chr>>11; _ghit.chr&=0x07FFU;
            //cout<<" "<<jj<<" "<<mm<<" "<<(int)_ghit.chr<<" "<<_ghit.loc<<endl;
            if(mm<ii) {cerr<<mm<<" "<<ii<<endl; exit(0);}
            if(mm==ii) {hits[ii][jj]=_ghit;}
            else if(mm>ii) {
                hits[ii][jj]=hits[ii][--_cur_n_hit[ii]];
                if(mm<=read_max_snp_num&&_cur_n_hit[mm]+_cur_n_chit[mm]<param.max_num_hits) 
                    hits[mm][_cur_n_hit[mm]++]=_ghit;
            }
        }
        for(jj=_cur_n_chit[ii]-1;jj>=0;--jj) {
            _ghit=chits[ii][jj]; mm=_ghit.chr>>11; _ghit.chr&=0x07FFU;
            //cout<<jj<<" "<<mm<<" "<<(int)_ghit.chr<<endl;
            if(mm<ii) {cerr<<mm<<" "<<ii<<endl; exit(0);}
            if(mm==ii) {chits[ii][jj]=_ghit;}
            else if(mm>ii) {
                chits[ii][jj]=chits[ii][--_cur_n_chit[ii]];
                if(mm<=read_max_snp_num&&_cur_n_hit[mm]+_cur_n_chit[mm]<param.max_num_hits) 
                    chits[mm][_cur_n_chit[mm]++]=_ghit;
            }                
        }
    }
}

//output align hits
void SingleAlign::StringAlign(RefSeq &ref, string &os) {
	bit32_t ii, sum=0, j=0;
	//snp align:
    //if(param.nt3) Process_3nt();
    for(ii=0; ii<=read_max_snp_num; ii++) {
        if((sum=_cur_n_hit[ii]+_cur_n_chit[ii])>0) break;
    }
    
    if(sum==0) s_OutHit(0, sum, ii, hits[0], 0, ref, os);
    else {
        j=myrand(_pread->index,&rand_rSeed)%sum;
   		if(j<_cur_n_hit[ii]) s_OutHit(0, sum, ii, &hits[ii][j], 0, ref, os);
   		else s_OutHit(1, sum, ii, &chits[ii][j-_cur_n_hit[ii]], 0, ref, os);
   	}
}

//write output according to types of hits
/* n: # of hits; chain: 0+/1-; flag: class of read; sig: 1, detect snp sites, 0, not */
void SingleAlign::s_OutHit(int chain, int n, bit8_t nsnps, gHit *hit, int insert_size, RefSeq &ref, string &os) {
    bit32_t ii, rev_seq, mm;
    mm=hit->chr>>11; hit->chr&=0x07FFU;
    rev_seq=chain^(hit->chr%2);
    if(param.out_sam){ //output in .sam format
        int flag; 
        flag=0x40*_pread->readset;
        if(n<0){
            if(!param.out_unmap) return;
            flag|=0x204; //QC
          	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);
        }
        else if(n==0||mm>read_max_snp_num) {
            if(!param.out_unmap) return;
            flag|=0x4; //NM
           	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);
        }
        else if(n>1&&param.report_repeat_hits==0) {
            if(!param.out_unmap) return;
            flag|=0x104; //NM
           	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);        
        }
        else{
            n_aligned++;
            if(n==1) flag|=0x0; //UM
            else flag|=0x100; //MA & OF
            if(rev_seq&&n) flag|=0x010; //reverse read seq
            //cout<<(int)hit->chr<<" "<<hit->loc<<" "<<(int)nsnps<<endl;
            if(hit->gap_size==0) sprintf(cigar, "%uM", map_readlen);
            else if(hit->gap_size>0) sprintf(cigar, "%dM%dD%dM", (int)hit->gap_pos, (int)hit->gap_size, (int)map_readlen-(int)hit->gap_pos);
            else sprintf(cigar, "%dM%dI%dM", (int)hit->gap_pos, -(int)hit->gap_size, (int)map_readlen-(int)hit->gap_pos+(int)hit->gap_size);
            sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d",_pread->name.c_str(),flag,ref.title[hit->chr].name.c_str(),hit->loc+1,cigar,_outseq[rev_seq].c_str(),_outqual[rev_seq].c_str(),nsnps);             
            os.append(_ch);
            if(param.out_ref) {
                int ptr=0;
                for(ii=2;ii>0;ii--) {
                    if(hit->loc<ii) continue;
                    _mapseq[ptr++]=param.useful_nt[ref.bfa[hit->chr&0xfffeU].s[(hit->loc-ii)/SEGLEN]>>(SEGLEN*2-2-((hit->loc-ii)%SEGLEN)*2)&0x3]+32;
                }
                for(ii=0;ii<map_readlen+2;ii++) {
                    _mapseq[ptr++]=param.useful_nt[ref.bfa[hit->chr&0xfffeU].s[(hit->loc+ii)/SEGLEN]>>(SEGLEN*2-2-((hit->loc+ii)%SEGLEN)*2)&0x3];
                }
                _mapseq[ptr]=0; _mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;
                sprintf(_ch, "\tXR:Z:%s",_mapseq);
                os.append(_ch);
            }
            
            if(param.RRBS_flag){
                seg_info=ref.CCGG_seglen(hit->chr, hit->loc, map_readlen);
                sprintf(_ch,"\tZP:i:%d\tZL:i:%d",seg_info.first,seg_info.second); 
                os.append(_ch);
            }

            sprintf(_ch,"\tZS:Z:%c%c\n",chain_flag[hit->chr%2], chain_flag[chain]); 
            os.append(_ch);
        }
    }

    else{ //output in .bsp format    
        if(!param.out_unmap&&(n<=0||(n>1&&param.report_repeat_hits==0))) return;
        sprintf(_ch, "%s\t", _pread->name.c_str());
    	os.append(_ch);
    	
    	sprintf(_ch, "%s\t%s\t", _outseq[rev_seq].c_str(), _outqual[rev_seq].c_str());
    	os.append(_ch);
    	
    	if (n<0) os.append("QC");
    	else if(n==0||mm>read_max_snp_num) os.append("NM"); 
    	else if(n==1) os.append("UM");
    	else if(n>=(int)param.max_num_hits) os.append("OF");
    	else os.append("MA");
        
        if(mm<=read_max_snp_num)
        if((n>0&&param.report_repeat_hits==1)||(n==1&&param.report_repeat_hits==0)){
            n_aligned++;
            int ptr=0;
            for(ii=2;ii>0;ii--) {
                if(hit->loc<ii) continue;
                _mapseq[ptr++]=param.useful_nt[ref.bfa[hit->chr&0xfffeU].s[(hit->loc-ii)/SEGLEN]>>(SEGLEN*2-2-((hit->loc-ii)%SEGLEN)*2)&0x3]+32;
            }
            for(ii=0;ii<map_readlen+2;ii++) {
                _mapseq[ptr++]=param.useful_nt[ref.bfa[hit->chr&0xfffeU].s[(hit->loc+ii)/SEGLEN]>>(SEGLEN*2-2-((hit->loc+ii)%SEGLEN)*2)&0x3];
            }
            _mapseq[ptr]=0; _mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;

        	sprintf(_ch, "\t%s\t%u\t%c%c\t%d\t%s\t", ref.title[hit->chr].name.c_str(), hit->loc+1, chain_flag[hit->chr%2],  chain_flag[chain], insert_size, _mapseq);
            os.append(_ch); 

            if(hit->gap_size) sprintf(_ch, "%d:%d:%d\t", nsnps, (int)hit->gap_size, (int)hit->gap_pos); 
            else sprintf(_ch, "%d\t", nsnps);
            os.append(_ch);
        
            for(ii=0; ii<read_max_snp_num; ii++){
                sprintf(_ch, "%d:", _cur_n_hit[ii]+_cur_n_chit[ii]);
                os.append(_ch);
            }
            sprintf(_ch, "%d", _cur_n_hit[ii]+_cur_n_chit[ii]);
            os.append(_ch);
        }
    	os.append("\n");
    }	

    //1	QNAME	Query (pair) NAME
    //2	FLAG	bitwise FLAG
    //3	RNAME	Reference sequence NAME
    //4	POS	1-based leftmost POSition/coordinate of clipped sequence
    //5	MAPQ	MAPping Quality (Phred-scaled)
    //6	CIAGR	extended CIGAR string
    //7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
    //8	MPOS	1-based Mate POSistion
    //9	ISIZE	Inferred insert SIZE
    //10	SEQ	query SEQuence on the same strand as the reference
    //11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
    //12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

}


void SingleAlign::Fix_Unpaired_Short_Fragment(RefSeq &ref) {
    bit32_t ii, j, k;
    if(map_readlen>=param.min_insert) return;
    for(ii=0; ii<=read_max_snp_num; ii++) {
        for(j=0;j<_cur_n_hit[ii];j++) {
            seg_info=ref.CCGG_seglen(hits[ii][j].chr, hits[ii][j].loc, map_readlen);
            if(seg_info.second<param.min_insert||seg_info.second>param.max_insert) {
                _cur_n_hit[ii]--; 
                for(k=j;k<_cur_n_hit[ii];k++) hits[ii][k]=hits[ii][k+1];
                j--;
            }
        }
        for(j=0;j<_cur_n_chit[ii];j++) {
            seg_info=ref.CCGG_seglen(chits[ii][j].chr, chits[ii][j].loc, map_readlen);
            if(seg_info.second<param.min_insert||seg_info.second>param.max_insert) {
                _cur_n_chit[ii]--; 
                for(k=j;k<_cur_n_chit[ii];k++) chits[ii][k]=chits[ii][k+1];
                j--;
            }
        }
        if(_cur_n_hit[ii]+_cur_n_chit[ii]>0) break;
    }
}
