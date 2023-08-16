

from Bio import SeqIO
from Bio.SeqUtils import GC
from dinucl_freq import *


def main():
    wsize = 1000       # size of sliding window
    ssize = 500        # step size of sliding window


    fileWindows = open("windows_out.txt", "w")   # file to write indices values for sliding windows


    # header for BED file
    fileWindows.write("seq\tstart\tend\tproduct\tsubstrate\tcomposite\tgc\tresult\n")


    for record in SeqIO.parse("example/NC_063025.1.fasta", "fasta"):    # for each fasta sequence
    
        w_start = 0                  # start of first window
        w_end   = w_start+wsize      # end of first window
    

        while w_start < len(record.seq):        # until the start of window reaches the sequence's end

            wseq = record.seq[w_start:w_end]    # get sequence from window
            
            wseq_subs = round( calculate_substrate(wseq.upper()), ndigits=2)    # calculate substrate index for window
            wseq_prod = round( calculate_product(wseq.upper()), ndigits=2)      # calculate product index for window
            wseq_comp = round( calculate_composite(wseq.upper()), ndigits=2)    # calculate composite index for window
            wseq_gc   = round( GC(wseq), ndigits=2)                             # get GC% for window
            rip_res   = classify_rip_window(wseq_subs, wseq_prod, wseq_comp)    # find out if window is affected or not by RIP
            
            # write to indices values to file
            fileWindows.write(record.id + "\t" + str(w_start) + "\t" + str(w_end) + "\t" + str(wseq_prod) + "\t" + str(wseq_subs) + "\t" + str(wseq_comp) + "\t" + str(wseq_gc) + "\t" + rip_res + "\n")
            

            if w_end == len(record.seq):    # Break here if the window's end coordinate reaches the sequence's end
                break

            w_start = w_start+ssize                               # start coordinate of next window
            w_end   = min((w_start + wsize), len(record.seq) )    # end coordinate of next window


    
    fileWindows.close()

if __name__ == "__main__":
    main()



