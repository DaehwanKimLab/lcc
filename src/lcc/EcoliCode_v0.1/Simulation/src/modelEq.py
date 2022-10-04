
def saturation(conc, KM):
    """ Calculate the relative saturation of substrate
    Output will be within range (0,1)
    """
    return (conc / (KM + conc))
    #return (1 + (KM/conc)) ** -1

def alloAct(conc, Ka):
    """ Calculate the activation effect of allosteric modifier"""
    return(1 + (conc/Ka)) 

def alloInhib(conc,Ki):
    """ Calculate the inhibition effect of allosteric modifier"""
    return (alloAct(conc, Ki) ** -1)