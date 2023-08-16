

from Bio import SeqIO
from Bio.SeqUtils import GC


# function that calculates the substrate index of a DNA sequence
def calculate_substrate(seq):
    CpA = seq.count('CA')
    TpG = seq.count('TG')
    ApC = seq.count('AC')
    GpT = seq.count('GT')
    
    if ApC + GpT > 0:
        return((CpA + TpG)/(ApC + GpT))
    else:
        return(1)
    
# function that calculates the product index of a DNA sequence
def calculate_product(seq):
    TpA = seq.count('TA')
    ApT = seq.count('AT')
       
    if ApT > 0:
        return(TpA/ApT)
    else:
        return(-1)


# function that calculates the composite index of a DNA sequence
def calculate_composite(seq):
    TpA = seq.count('TA')
    ApT = seq.count('AT')
    CpA = seq.count('CA')
    TpG = seq.count('TG')
    ApC = seq.count('AC')
    GpT = seq.count('GT')

    if ApT > 0 and ApC + GpT > 0:
        return( (TpA/ApT) - ((CpA + TpG) / (ApC + GpT)) )
    else:
        return(-1)

# function that decides whether the given values of substrate, product, and composite indices imply evidence of RIP
def classify_rip_window(w_substrate, w_product, w_composite):
    
    # thresholds
    max_substrate = 0.75    # values less than 0.75 are RIP affected
    min_product   = 1.1     # values greater than 1.1 are RIP affected
    min_composite = 0.01    # a positive RIP composite is RIP affected
    
    if w_substrate <= max_substrate and w_product >= min_product and w_composite >= min_composite:
        return("RIP_affected")
    else:
       return("RIP_notaffected")



min_LRAR = 4000    # minimum size of LRARs to consider
LRAR_count = 1     # starting LRAR counter

wsize = 1000       # size of sliding window
ssize = 500        # step size of sliding window


fileLRARs   = open("LRARs.bed", "w")         # file to write LRAR locations
fileWindows = open("windows_out.txt", "w")   # file to write indices values for sliding windows


# header for BED file
fileWindows.write("seq\tstart\tend\tproduct\tsubstrate\tcomposite\tgc\tresult\n")


# for each fasta sequence
for record in SeqIO.parse("PB655_02_Race4.bp.p_ctg_sorted_edited.fasta", "fasta"):
    
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
        
        
        if rip_islant_status == 'RIP_notaffected' and rip_res == "RIP_affected":
            rip_island_start = w_start
            rip_islant_status = "RIP_affected"
            
        if rip_islant_status == "RIP_affected" and rip_res == "RIP_affected":
            rip_island_end = w_end
        
        if rip_islant_status == "RIP_affected" and rip_res == "RIP_notaffected":
            rip_islant_status = "RIP_notaffected"
            if rip_island_end - rip_island_start >= min_LRAR:
                fileLRARs.write(record.id + "\t" + str(rip_island_start) + "\t" + str(rip_island_end) + "\tLRAR_" + str(LRAR_count) + "\n")
                LRAR_count+=1
            

        w_start = w_start+ssize
        w_end   = min((w_start + wsize), len(record.seq) )
        
        
        

    if rip_islant_status == "RIP_affected":
            rip_island_end = w_start
            rip_islant_status = "RIP_notaffected"
            if rip_island_end - rip_island_start >= min_LRAR:
                fileLRARs.write(record.id + "\t" + str(rip_island_start) + "\t" + str(rip_island_end) + "\tLRAR_" + str(LRAR_count) + "\n")
                LRAR_count+=1


fileLRARs.close()
fileWindows.close()