#Absolute codon counts.

'''The following script evaluates codon occurrences and represents them as a bar graph. 
   Alternatively, you can save the images on your system instead of displaying them instantly.'''


import csv
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC,generic_dna
import matplotlib.pyplot as plt


table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'O', 'TAG':'O', 
        'TGC':'C', 'TGT':'C', 'TGA':'O', 'TGG':'W' }     #Initially provide the codon-amino acid table.


AA_0 = dict.fromkeys(table.values(),0 )   #dict of amino acid count

'''The following is used only to read and write to a file.'''
genes=['Arca','CRP','FnR','Fur','hns','rob','soxS']
h=5  #Hamming Distance as required
f_w=csv.writer(open('/Users/sanatmishra27/Desktop/'+'xylR'+'/'+'Ham'+str(h)+'-'+'xylR'+'CodonFreq.csv','w+'))
f_w.writerow(["codonfrequency"])
f_open=open('/Users/sanatmishra27/Desktop/SUMMER PROJECT \'19/Plots/SynonymousSeq/Local/'+'xylR'+str('/Ham')+str(h)+'-'+'xylR'+str('.csv'),"r")
r=csv.reader(f_open)
next(r)

'''The actual implementation begins here. The code runs for each sequence provided in the file.'''
for row in r:
    C_0=dict.fromkeys(table.keys(),0)     #dict of codon counts.
  
    G=row[0].replace('U','T')  #reads the sequence off the file.
    a1=[G[i:i+3] for i in range(0,len(G),3)]   #a1 contains the codons in the sequence.
    for i in a1:
        i=Seq(i)
        if(i.translate() in AA_0):
            AA_0[i.translate()]+=1   #Increase amino acid count.

        if(i in C_0):
            C_0[i]+=1                #Increase codon count.
            
                 
    '''Now that we have the amino acid counts, next step is preparing the data so as to plot it.'''

    x1=[]
    y1=[]
    amacid_codes=[]
    for item in AA_0:
        for z in C_0:
            z=Seq(z)
            if(z.translate()==item):
                amacid_codes.append(item)
                x1.append(str(z))   #contains the codon names.
                y1.append(C_0[z])   #contains number of codon occurrences.
                
  
    fig = plt.figure()
    
    '''Setting the dimensions for the plot.'''
    fig.set_figheight(15)
    fig.set_figwidth(25)

    plt.xticks(rotation=90,size=35) #Specify size of points on x axis.
    plt.yticks(size=35)             #Specify size of points on y axis. 

    for i in range(len(x1)):
        plt.text(x=i-0.4,y=y1[i]+0.5, s = amacid_codes[i], size = 35,color='red')   #This labels each bar with the amino acid code it represents.
    
    plt.xlabel('Codons',size=35)
    plt.ylabel('No. of occurrences.',size=35)
    plt.bar(x1, y1)  #Plot bar graph.
    #plt.savefig('Insert path of file where image will be stored'.png) Use this in case you want to save the image directly.
    plt.show();
