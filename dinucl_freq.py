

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