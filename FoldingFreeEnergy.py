#FoldingFreeEnergy

'''Given a mRNA sequence, RNAFold can calculate the folding free energy of the mRNA secondary structure.
   Please donwload the library on your system from 'https://www.tbi.univie.ac.at/RNA/'. '''

import RNA

mrna='ATGC'  #input sequence.
mfe = RNA.fold(str(mrna))
print(float('%6.2f'%mfe[1]))   #Contain to two-place decimal.
