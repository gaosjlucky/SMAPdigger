import sys, time, os, array, optparse
usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output file name. (required)", default="")
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
parser.add_option("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
parser.add_option("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
parser.add_option("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
parser.add_option("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios.", default=False)
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
parser.add_option("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
parser.add_option("-t", "--trim-fillin", dest="trim_fillin", type="int", metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 2]", default=2)
parser.add_option("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
parser.add_option("-m", "--min-depth", dest="min_depth", type="int", metavar='FOLD', help="report loci with sequencing depth>=FOLD. [default: 1]", default=1)
parser.add_option("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
parser.add_option("-i", "--ct-snp", dest="CT_SNP", help='how to handle CT SNP ("no-action", "correct", "skip"), default: "correct".', default="correct")

options, infiles = parser.parse_args()

if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(options.outfile) == 0: parser.error("Missing output file name, use -o or --out option.")
if len(infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if any(options.chroms): options.chroms = options.chroms.split(',')
CT_SNP_val = {"no-action": 0, "correct": 1, "skip": 2}
try: options.CT_SNP = CT_SNP_val[options.CT_SNP.lower()]
except: parser.error('Invalid -i value, select "no-action", "correct" or "skip"')
if options.min_depth <= 0: parser.error('Invalid -m value, must >= 1')
if options.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')

if any(options.sam_path): 
    if options.sam_path[-1] != '/': options.sam_path += '/'

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, ''.join(['\t' for i in xrange(nt)]+['@ ',time.asctime(),': ',txt])

def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        flag = col[1] 
        if 'u' in flag: return []
        if options.unique and 's' in flag: return []
        if options.pair and 'P' not in flag: return []
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
        if cr not in options.chroms: return []
        for aux in col[11:]:
            if aux[:5] == 'ZS:Z:': 
                strand = aux[5:7]
                break
        assert strand, 'missing strand information "ZS:Z:xx"'
        gap_pos, gap_size = 0, 0
        while 'I' in cigar or 'D' in cigar:
            for sep in 'MID':
                try: gap_size = int(cigar.split(sep, 1)[0])
                except ValueError: continue
                break
            if sep == 'M': gap_pos += gap_size
            elif sep == 'I': seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            elif sep == 'D': 
                seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]                        
                gap_pos += gap_size
            cigar = cigar[cigar.index(sep)+1:]
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if options.unique and flag != 'UM': return []
        if options.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in options.chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if options.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if options.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-': seq = seq[:-options.trim_fillin]
        elif strand == '--': seq, pos = seq[options.trim_fillin:], pos+options.trim_fillin
        elif insert != 0 and len(seq) > abs(insert) - options.trim_fillin:
            trim_nt = len(seq) - (abs(insert) - options.trim_fillin)
            if strand == '++': seq = seq[:-trim_nt]
            elif strand == '-+': seq, pos =seq[trim_nt:], pos+trim_nt
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (seq, strand[0], cr, pos)

ref, cr, seq = {}, '', ''
disp('reading reference %s ...' % options.reffile)
for line in open(options.reffile):
    if line[0] == '>': 
        if any(cr): 
            if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
del seq

meth, depth, coverage, meth1, depth1 = {}, {}, {}, {}, {}
for cr in ref:
    meth[cr] = array.array('H', [0]) * len(ref[cr])
    depth[cr] = array.array('H', [0]) * len(ref[cr])
    if options.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])
    if options.CT_SNP > 0:
        meth1[cr] = array.array('H', [0]) * len(ref[cr])
        depth1[cr] = array.array('H', [0]) * len(ref[cr])

options.chroms = set(ref.keys())

BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
nmap = 0
for infile in infiles:
    nline = 0
    disp('reading %s ...' % infile)
    if infile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (options.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (options.sam_path, infile))
    else: sam_format, fin = False, open(infile)
    for line in fin:
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info
        depthcr = depth[cr]
        if pos + len(seq) > len(depthcr): continue
        nmap += 1
        methcr = meth[cr]
        refseq = ref[cr][pos:pos+len(seq)]
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        index = refseq.find(match)
        while index >= 0:
            if depthcr[pos+index] < 65535:
                if seq[index] == convert: depthcr[pos+index] += 1
                elif seq[index] == match: 
                    methcr[pos+index] += 1
                    depthcr[pos+index] += 1
            index = refseq.find(match, index+1)
        if options.CT_SNP == 0: continue
        methcr1 = meth1[cr]
        depthcr1 = depth1[cr]
        index = refseq.find(rc_match)
        while index >= 0:
            if depthcr1[pos+index] < 65535:
                if seq[index] == rc_convert: depthcr1[pos+index] += 1
                if seq[index] == rc_match: 
                    methcr1[pos+index] += 1
                    depthcr1[pos+index] += 1
            index = refseq.find(rc_match, index+1)
    
    fin.close()

if options.combine_CpG:
    disp('combining CpG methylation from both strands ...')
    for cr in depth:
        depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
        if options.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
        pos = refcr.find('CG')
        while pos >= 0:
            if depthcr[pos] + depthcr[pos+1] <= 65535:
                depthcr[pos] += depthcr[pos+1]
                methcr[pos] += methcr[pos+1]
            else:
                depthcr[pos] = (depthcr[pos] + depthcr[pos+1]) / 2
                methcr[pos] = (methcr[pos] + methcr[pos+1]) / 2
            depthcr[pos+1] = 0
            methcr[pos+1] = 0
            if options.CT_SNP > 0:
                if depthcr1[pos] + depthcr1[pos+1] <= 65535:
                    depthcr1[pos] += depthcr1[pos+1]
                    methcr1[pos] += methcr1[pos+1]
                else:
                    depthcr1[pos] = (depthcr1[pos] + depthcr1[pos+1]) / 2
                    methcr1[pos] = (methcr1[pos] + methcr1[pos+1]) / 2
            pos = refcr.find('CG', pos+2)

disp('writing %s ...' % options.outfile)
ss = {'C': '+', 'G': '-'}
fout = open(options.outfile, 'w')
if not options.no_header: 
    fout.write('chr\tpos\tstrand\tcontext\tratio\teff_CT_count\tC_count\tCT_count\trev_G_count\trev_GA_count\tCI_lower\tCI_upper\n')
z95, z95sq = 1.96, 1.96 * 1.96
nc, nd, dep0 = 0, 0, options.min_depth
for cr in sorted(depth.keys()):
    depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
    if options.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
    for i, dd in enumerate(depthcr):
        if dd < dep0: continue
        if options.CT_SNP > 0: 
            m1, d1 = methcr1[i], depthcr1[i]
            if m1 != d1:
                if options.CT_SNP == 2: continue
                d = float(dd) * m1 / d1
            else: d = float(dd)
        else: d = float(dd)
        nc += 1
        nd += d
        m = methcr[i]
        if m == 0 and not options.meth0: continue
        seq = refcr[i-2:i+3]
        try: strand = ss[refcr[i]]
        except KeyError: continue
        try: ratio = float(min(m,d)) / d
        except ZeroDivisionError:
            if options.CT_SNP:
                fout.write('%s\t%d\t%c\t%s\tNA\t%.2f\t%d\t%d\t%d\t%d\tNA\tNA\n' % (cr, i+1, strand, seq, d, m, dd, m1, d1))
            else:
                fout.write('%s\t%d\t%c\t%s\tNA\t%.2f\t%d\t%d\tNA\tNA\tNA\tNA\n' % (cr, i+1, strand, seq, d, m, dd))
            continue     
        pmid = ratio + z95sq / (2 * d)
        sd = z95 * ((ratio*(1-ratio)/d + z95sq/(4*d*d)) ** 0.5)
        norminator = 1 + z95sq / d
        CIl, CIu = (pmid - sd) / norminator, (pmid + sd) / norminator
        if options.CT_SNP:
            fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' 
                % (cr, i+1, strand, seq, ratio, d, m, dd, m1, d1, CIl, CIu))
        else:
            fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\tNA\tNA\t%.3f\t%.3f\n'
                % (cr, i+1, strand, seq, ratio, d, m, dd, CIl, CIu))
                            
fout.close()
disp('done.')
if nc > 0: print 'total %d valid mappings, %d covered cytosines, average coverage: %.2f fold.' % (nmap, nc, float(nd)/nc)
