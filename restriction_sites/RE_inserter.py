
from fuzzysearch import find_near_matches
from typing import List, Dict

# Human
freqs = {
    'TTT':  17.1,  'TAT':  12.1,  'TCT':  16.9,  'TGT':  10.4,  'ATT':  16.4,  'AAT':  18.4,  'ACT':  14.2,  'AGT':  14.0,
    'TTC':  17.4,  'TAC':  13.4,  'TCC':  17.3,  'TGC':  10.8,  'ATC':  18.6,  'AAC':  18.3,  'ACC':  17.8,  'AGC':  19.6,
    'TTA':   8.7,  'TAA':   0.4,  'TCA':  14.1,  'TGA':   0.7,  'ATA':   8.0,  'AAA':  27.4,  'ACA':  16.5,  'AGA':  13.2,
    'TTG': 13.44,  'TAG':   0.3,  'TCG':  4.03,  'TGG':  11.6,  'ATG': 21.53,  'AAG':  31.7,  'ACG':  5.59,  'AGG':  12.1,
    'CTT':  14.0,  'CAT':  11.8,  'CCT':  19.3,  'CGT':   4.5,  'GTT':  11.7,  'GAT':  24.0,  'GCT':  18.9,  'GGT':  10.8,
    'CTC':  17.8,  'CAC':  14.6,  'CCC':  19.1,  'CGC':   8.7,  'GTC':  13.4,  'GAC':  24.2,  'GCC':  25.8,  'GGC':  19.7,
    'CTA':   7.4,  'CAA':  14.0,  'CCA':  18.9,  'CGA':   6.4,  'GTA':   7.6,  'GAA':  33.6,  'GCA':  17.0,  'GGA':  17.1,
    'CTG': 36.10,  'CAG':  35.5,  'CCG':  6.22,  'CGG':  10.7,  'GTG':  25.8,  'GAG':  39.6,  'GCG':   5.9,  'GGG':  15.3,
}

code = {
    'TTT':   'F',  'TAT':   'Y',  'TCT':   'S',  'TGT':   'C',  'ATT':   'I',  'AAT':   'N',  'ACT':   'T',  'AGT':   'S',
    'TTC':   'F',  'TAC':   'Y',  'TCC':   'S',  'TGC':   'C',  'ATC':   'I',  'AAC':   'N',  'ACC':   'T',  'AGC':   'S',
    'TTA':   'L',  'TAA':   '*',  'TCA':   'S',  'TGA':   '*',  'ATA':   'I',  'AAA':   'K',  'ACA':   'T',  'AGA':   'R',
    'TTG':   'L',  'TAG':   '*',  'TCG':   'S',  'TGG':   'W',  'ATG':   'M',  'AAG':   'K',  'ACG':   'T',  'AGG':   'R',
    'CTT':   'L',  'CAT':   'H',  'CCT':   'P',  'CGT':   'R',  'GTT':   'V',  'GAT':   'D',  'GCT':   'A',  'GGT':   'G',
    'CTC':   'L',  'CAC':   'H',  'CCC':   'P',  'CGC':   'R',  'GTC':   'V',  'GAC':   'D',  'GCC':   'A',  'GGC':   'G',
    'CTA':   'L',  'CAA':   'Q',  'CCA':   'P',  'CGA':   'R',  'GTA':   'V',  'GAA':   'E',  'GCA':   'A',  'GGA':   'G',
    'CTG':   'L',  'CAG':   'Q',  'CCG':   'P',  'CGG':   'R',  'GTG':   'V',  'GAG':   'E',  'GCG':   'A',  'GGG':   'G',
}

aas = set(code.values())
aa_codons = {aa: [(k, v) for k, v in freqs.items() if code[k] == aa] for aa in aas}
sorted_codons = {
    aa: [c for c, f in sorted(codons, key=lambda c: c[1], reverse=True)] for aa, codons in aa_codons.items()
}

def protein_to_DNA(seq: str):
    return ''.join([sorted_codons[aa][0] for aa in seq])

RE_sites = {
    'EcoRI': 'GAATTC',
    'NheI':  'GCTAGC',
    'BamHI': 'GGATCC',
    'BstBI': 'TTCGAA',
    'ClaI':  'ATCGAT',
    'BsiWI': 'CGTACG',
    'AatII': 'GACGTC',
    'MscI':  'TGGCCA',
    'SwaI':  'ATTTAAAT',
    'BsaI':  'GGTCTC',
    'XhoI':  'CTCGAG',
    'XbaI':  'TCTAGA',
}

def rev_comp(seq):
    nts = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    tr = str.maketrans(nts)
    return seq[::-1].translate(tr)

def try_insert_site(site, seq):
    candidates = []
    rr = find_near_matches(site, seq, max_substitutions=len(site), max_insertions=0, max_deletions=0)
    for m in rr:
        s = m.start - m.start % 3
        e = s + 6 if m.start % 3 == 0 else s + 9
        old_aas = [code[seq[i:i+3]] for i in range(s, e, 3)]
        new_seq = seq[:m.start] + site + seq[m.end:]
        new_aas = [code[new_seq[i:i+3]] for i in range(s, e, 3)]
        if (all(aa1 == aa2 for aa1, aa2 in zip(old_aas, new_aas))):
            candidates.append(new_seq)
    return(candidates)

def try_insert_sites(seq: str, sites: Dict[str, str] = RE_sites, protein_seq: bool = False) -> Dict[str, List[str]]:
    '''
    Tries to modify the coding sequence of a protein to include a restriction enzyme site

    Args
        seq: Coding sequence of a protein, must be in +1 frame. Interpreted as a protein sequence instead if protein_seq is True
        sites: A dict mapping a restriction enzyme name to its recognition site, e.g. {'EcoRI': 'GAATTC'}
        protein_seq: If True, seq is interpreted as a protein sequence
    
    Returns
        A dict mapping restriction enzyme names to a list of modified coding sequences
    '''

    if protein_seq:
        seq = protein_to_DNA(seq)
    res = {}
    for name, site in sites.items():
        candidates = try_insert_site(site, seq)
        # checking complementary strand for non-palindromic recognition sites
        if (rc := rev_comp(site)) != site:
            candidates.extend(try_insert_site(rc, seq))
        if len(candidates) > 0:
            res[name] = candidates
    return res

if __name__ == '__main__':
    # An example for T2A peptide
    from pprint import PrettyPrinter
    T2A = 'EGRGSLLTCGDVEENPGP'
    PrettyPrinter().pprint(try_insert_sites(T2A, protein_seq=True))
    # An example of a long protein
    GCN1 = 'MAADTQVSETLKRFAGKVTTASVKERREILSELGKCVAGKDLPEGAVKGLCKLFCLTLHRYRDAASRRALQAAIQQLAEAQPEATAKNLLHSLQSSGIGSKAGVPSKSSGSAALLALTWTCLLVRIVFPSRAKRQGDIWNKLVEVQCLLLLEVLGGSHKHAVDGAVKKLTKLWKENPGLVEQYLSAILSLEPNQNYAGMLGLLVQFCTSHKEMDVVSQHKSALLDFYMKNILMSKVKPPKYLLDSCAPLLRYLSHSEFKDLILPTIQKSLLRSPENVIETISSLLASVTLDLSQYAMDIVKGLAGHLKSNSPRLMDEAVLALRNLARQCSDSSAMESLTKHLFAILGGSEGKLTVVAQKMSVLSGIGSVSHHVVSGPSSQVLNGIVAELFIPFLQQEVHEGTLVHAVSVLALWCNRFTMEVPKKLTEWFKKAFSLKTSTSAVRHAYLQCMLASYRGDTLLQALDLLPLLIQTVEKAASQSTQVPTITEGVAAALLLLKLSVADSQAEAKLSSFWQLIVDEKKQVFTSEKFLVMASEDALCTVLHLTERLFLDHPHRLTGNKVQQYHRALVAVLLSRTWHVRRQAQQTVRKLLSSLGGFKLAHGLLEELKTVLSSHKVLPLEALVTDAGEVTEAGKAYVPPRVLQEALCVISGVPGLKGDVTDTEQLAQEMLIISHHPSLVAVQSGLWPALLARMKIDPEAFITRHLDQIIPRMTTQSPLNQSSMNAMGSLSVLSPDRVLPQLISTITASVQNPALRLVTREEFAIMQTPAGELYDKSIIQSAQQDSIKKANMKRENKAYSFKEQIIELELKEEIKKKKGIKEEVQLTSKQKEMLQAQLDREAQVRRRLQELDGELEAALGLLDIILAKNPSGLTQYIPVLVDSFLPLLKSPLAAPRIKNPFLSLAACVMPSRLKALGTLVSHVTLRLLKPECVLDKSWCQEELSVAVKRAVMLLHTHTITSRVGKGEPGAAPLSAPAFSLVFPFLKMVLTEMPHHSEEEEEWMAQILQILTVQAQLRASPNTPPGRVDENGPELLPRVAMLRLLTWVIGTGSPRLQVLASDTLTTLCASSSGDDGCAFAEQEEVDVLLCALQSPCASVRETVLRGLMELHMVLPAPDTDEKNGLNLLRRLWVVKFDKEEEIRKLAERLWSMMGLDLQPDLCSLLIDDVIYHEAAVRQAGAEALSQAVARYQRQAAEVMGRLMEIYQEKLYRPPPVLDALGRVISESPPDQWEARCGLALALNKLSQYLDSSQVKPLFQFFVPDALNDRHPDVRKCMLDAALATLNTHGKENVNSLLPVFEEFLKNAPNDASYDAVRQSVVVLMGSLAKHLDKSDPKVKPIVAKLIAALSTPSQQVQESVASCLPPLVPAIKEDAGGMIQRLMQQLLESDKYAERKGAAYGLAGLVKGLGILSLKQQEMMAALTDAIQDKKNFRRREGALFAFEMLCTMLGKLFEPYVVHVLPHLLLCFGDGNQYVREAADDCAKAVMSNLSAHGVKLVLPSLLAALEEESWRTKAGSVELLGAMAYCAPKQLSSCLPNIVPKLTEVLTDSHVKVQKAGQQALRQIGSVIRNPEILAIAPVLLDALTDPSRKTQKCLQTLLDTKFVHFIDAPSLALIMPIVQRAFQDRSTDTRKMAAQIIGNMYSLTDQKDLAPYLPSVTPGLKASLLDPVPEVRTVSAKALGAMVKGMGESCFEDLLPWLMETLTYEQSSVDRSGAAQGLAEVMAGLGVEKLEKLMPEIVATASKVDIAPHVRDGYIMMFNYLPITFGDKFTPYVGPIIPCILKALADENEFVRDTALRAGQRVISMYAETAIALLLPQLEQGLFDDLWRIRFSSVQLLGDLLFHISGVTGKMTTETASEDDNFGTAQSNKAIITALGVERRNRVLAGLYMGRSDTQLVVRQASLHVWKIVVSNTPRTLREILPTLFGLLLGFLASTCADKRTIAARTLGDLVRKLGEKILPEIIPILEEGLRSQKSDERQGVCIGLSEIMKSTSRDAVLYFSESLVPTARKALCDPLEEVREAAAKTFEQLHSTIGHQALEDILPFLLKQLDDEEVSEFALDGLKQVMAIKSRVVLPYLVPKLTTPPVNTRVLAFLSSVAGDALTRHLGVILPAVMLALKEKLGTPDEQLEMANCQAVILSVEDDTGHRIIIEDLLEATRSPEVGMRQAAAIILNIYCSRSKADYTSHLRSLVSGLIRLFNDSSPVVLEESWDALNAITKKLDAGNQLALIEELHKEIRLIGNESKGEHVPGFCLPKKGVTSILPVLREGVLTGSPEQKEEAAKALGLVIRLTSADALRPSVVSITGPLIRILGDRFSWNVKAALLETLSLLLAKVGIALKPFLPQLQTTFTKALQDSNRGVRLKAADALGKLISIHIKVDPLFTELLNGIRAMEDPGVRDTMLQALRFVIQGAGAKVDAVIRKNIVSLLLSMLGHDEDNTRISSAGCLGELCAFLTEEELSAVLQQCLLADVSGIDWMVRHGRSLALSVAVNVAPGRLCAGRYSSDVQEMILSSATADRIPIAVSGVRGMGFLMRHHIETGGGQLPAKLSSLFVKCLQNPSSDIRLVAEKMIWWANKDPLPPLDPQAIKPILKALLDNTKDKNTVVRAYSDQAIVNLLKMRQGEEVFQSLSKILDVASLEVLNEVNRRSLKKLASQADSTEQVDDTILT'
    GCN1_RE = try_insert_sites(GCN1, protein_seq=True)
    PrettyPrinter().pprint({k: len(v) for k, v in GCN1_RE.items()})