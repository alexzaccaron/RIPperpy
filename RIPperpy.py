

from Bio import SeqIO
from Bio.SeqUtils import GC
from dinucl_freq import *


def main():
    wsize = 1000       # size of sliding window
    ssize = 500        # step size of sliding window


    fileWindows = open("windows_out.txt", "w")   # file to write indices values for sliding windows


    # header for BED file
    fileWindows.write("seq\tstart\tend\tproduct\tsubstrate\tcomposite\tgc\tresult\n")


    # for each fasta sequence
    for record in SeqIO.parse("example/NC_063025.1.fasta", "fasta"):
    
        w_start = 0
        w_end   = w_start+wsize
    
    
        rip_islant_status = 'RIP_notaffected'
        rip_island_start = 0
        rip_island_end   = 0
    
        while w_start < len(record.seq):
        
             
            wseq = record.seq[w_start:w_end]    # get sequence from window
            
            
            wseq_subs = round( calculate_substrate(wseq.upper()), ndigits=2)
            wseq_prod = round( calculate_product(wseq.upper()), ndigits=2)
            wseq_comp = round( calculate_composite(wseq.upper()), ndigits=2)
            wseq_gc   = round( GC(wseq), ndigits=2)
            rip_res   = classify_rip_window(wseq_subs, wseq_prod, wseq_comp)
            
            
            fileWindows.write(record.id + "\t" + str(w_start) + "\t" + str(w_end) + "\t" + str(wseq_prod) + "\t" + str(wseq_subs) + "\t" + str(wseq_comp) + "\t" + str(wseq_gc) + "\t" + rip_res + "\n")
                

            w_start = w_start+ssize
            w_end   = min((w_start + wsize), len(record.seq) )
    
    fileWindows.close()

if __name__ == "__main__":
    main()



