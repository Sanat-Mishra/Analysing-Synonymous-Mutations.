#Genetic Algo

'''A genetic algorithm (GA) mimics the process of natural selection to generate 
solutions given a set of constraints. Given a target solution, the GA produces 
solutions (sequences, in my case) trying to maximise a fitness function until it reaches the target.
Sub-optimal solution are removed, justas through selection.
I have written functions specifically for each component-
1) Generating the target sequence.
2) Introduce random mutations.
3) Initialise parent sequences.
4) Crossing-over.

I have used Hamming Distance as the fitness function, it is used as a fitness measure to generate the target sequence.
'''

from random import random,randint,sample,choice,shuffle,randrange
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import itertools
import csv

'''The following snippet evaluates the Hamming Distance between the sequnces, 
i.e. the number of positions between two strings where they differ. '''

def fitnesseval(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

'''Next, I wanted to find out a sequence containing only synonymous mutations 
(with respect to an input sequence) with the maximum possible Hamming Distance. 
This will become the target sequence for the Genetic Algorithm. '''

def maxHam(inp):
    ogseq=inp.replace('T','U')   #Replace 'T' of the coding sequence with 'U to convert it to mRNA'.
    
    messenger_rna= Seq(ogseq,IUPAC.unambiguous_rna)  #Used BioPython module for translation. Here messenger_rna converts it to a Sequence Object BioPython can recognize. 
    
    amacid=messenger_rna.translate() # amacid now contains the translated form of the input sequnce.
    
    messenger_rna=str(ogseq) #now messenger_rna contains the original sequence
    
    l=messenger_rna[0]+messenger_rna[2]  #l contains the first and third base of the codon. Since by keeping the second base same, almost all synonymous variants for a codon can be generated.
    
    codon=Seq(messenger_rna,IUPAC.unambiguous_rna)
    amacid=str(codon.translate())
    lis=(list(itertools.permutations('AUGC',2)))  #lis is a list containing all 2 base permutation-tuples.
    
    lis.extend([('A','A'),('U','U'),('G','G'),('C','C')]) #now lis includes all possible permutations- with repitition.
    
    shuffle(lis) #lis is now shuffled to ensure while picking element from a list, each time new element is picked.  
   
    for ele in range(0,len(lis)):
        
        '''If amacid, i.e. amino acid does not have synonymous variants, it will be returned as is.'''
        
        if amacid=='M': #singled out methionine since it is 1-fold degenerate (no syn. variant)
            check=ogseq
            break
        if amacid=='W': #singled out tryptophan since it is 1-fold degenerate (no syn. variant)
            check=ogseq
            break
        gt=lis[ele] 
        nt=''.join(gt) #To join the contents of gt (since it was a tuple)
        
        codonlist=list(messenger_rna)  #codonlist is a list of bases in te original sequence.
        
        '''In case the amino acid is Serine, I introduced manual changes to ensure a synoymous variant is produced '''
        
        if amacid=='S' and l!=nt:  # l!=nt ensures that l(original two bases) and nt(randomly chosen two bases) are not the same.
            codonlist[0]=nt[0] #assign the first element of codonlist as the first one in nt.
            if codonlist[1]=='G':
                codonlist[1]='C' 
            elif codonlist[1]=='C':
                codonlist[1]='G' 
            codonlist[len(codonlist)-1]=nt[1]
            newcodon=''.join(codonlist)
            newcodon= Seq(newcodon,IUPAC.unambiguous_rna)
            check=str(newcodon[0:len(newcodon)])
            
            if str(newcodon.translate())==amacid and fitnesseval(check,messenger_rna)==3:
                
                break
            
        codonlist[0]=nt[0]
        codonlist[len(codonlist)-1]=nt[1]
        newcodon=''.join(codonlist)
        newcodon= Seq(newcodon,IUPAC.unambiguous_rna)
        check=str(newcodon[0:len(newcodon)])  #check contains a new codon.
        
        '''Now we match the codon generated (check) with the original one (amacid). '''
        if str(newcodon.translate())==amacid and amacid!='S':
            '''Again, for Arginine and Leucine, things had to be dealt in a case specific manner.'''
            
            if amacid=='R' and l!=nt:
                if fitnesseval(check,messenger_rna)==2:
                    break
                else:
                    for cod in range(ele+1,len(lis)):
                        gt=lis[cod]
                        nt=''.join(gt)
                        codonlist[0]=nt[0]
                        codonlist[len(codonlist)-1]=nt[1]
                        newcodon=''.join(codonlist)
                        newcodon= Seq(newcodon,IUPAC.unambiguous_rna)
                        check=str(newcodon[0:len(newcodon)])
                        if fitnesseval(check,messenger_rna)==2 and newcodon.translate()==amacid:
                            break
                    break
            
            elif amacid=='L' :
                if fitnesseval(check,messenger_rna)==2:
                    break
                else:
                    for cod in range(ele+1,len(lis)):
                        gt=lis[cod]
                        nt=''.join(gt)
                        codonlist[0]=nt[0]
                        codonlist[len(codonlist)-1]=nt[1]
                        newcodon=''.join(codonlist)
                        newcodon= Seq(newcodon,IUPAC.unambiguous_rna)
                        check=str(newcodon[0:len(newcodon)])
                        if fitnesseval(check,messenger_rna)==2 and newcodon.translate()==amacid:
                            break
                    break
            elif l!=nt:
                    
                break 
                
    return check   #check is the synonymous variant with maximum Hamming distance .
    

      
'''Random crossing over, as happens in meiosis.'''
def crossover(p1,p2):  
    x=randrange(0,len(p1)-3,3)
  
    child=p1[:x] + p2[x:]
    return child


'''The following lines take an input DNA sequence, and convert it to a sequence with maximum Hamming Distance.'''
inp='ATGTCCCATCAGAAAATTATTCAGGATCTTATCGCATGGATTGACGAGCATATTGACCAGCCGCTTAACATTGATGTAGTCGCAAAAAAATCAGGCTATTCAAAGTGGTACTTGCAACGAATGTTCCGCACGGTGACGCATCAGACGCTTGGCGATTACATTCGCCAACGCCGCCTGTTACTGGCCGCCGTTGAGTTGCGCACCACCGAGCGTCCGATTTTTGATATCGCAATGGACCTGGGTTATGTCTCGCAGCAGACCTTCTCCCGCGTTTTCCGTCGGCAGTTTGATCGCACTCCCAGCGATTATCGCCACCGCCTGTAA'
inp=inp.replace('T','U')

valseq=''
for i in range(0,len(inp),3):
        j=min(i+3,len(inp))
        valseq+=maxHam(inp[i:j])

        
'''Introducing random mutations. '''

def mutation(inpv):
    l=[]
    done=[]
    r=0
    maxmut=randint(1,fitnesseval(inp,valseq))  #The maximum number of mutations is the Hamming Distance between the input sequence (inp) and one with generated by function maxHam (valseq).
    for i in range(0,len(inpv),3):
        l.append(inpv[i:i+3])  #l contains the codons in the sequence (inpv)
        
    
    while r<maxmut:
        pos=randint(0,len(l)-1) #The position at which the mutation will be introduced.
       
        l[pos]=maxHam(l[pos]) #Now l contains the synonymous codons, and as per the maxHam procedure.
        done.append(pos)
        r+=1

    return ''.join(l)

'''Initialising initial population, i.e. the parents.'''

def initialpop(inp):
    entseq=Seq(inp,IUPAC.unambiguous_rna)
    amacid=entseq.translate()
    parents=[]
    i=0
    while len(parents) <2:
        nt=['A','U','G','C']
        for i in range(len(nt)):
            newnt=nt[i]
            place=randint(0,len(inp)-1)
            if inp[place]!=newnt:
                lis=list(inp)
                lis[place]=newnt
                nlis=''.join(lis)
                nlis=str(nlis)
                newseq=Seq(nlis,IUPAC.unambiguous_rna)
                if newseq.translate()==amacid:
                    parents.append(nlis)
    parents.append(valseq)
    return parents       

'''Now we prepare the csv files onto which the data will be written. '''

f_w=csv.writer(open('/Users/sanatmishra27/Desktop/soxS-MG1655.csv','w'))
f_w.writerow(['Generation no.','Hamming Dist','Sequence'])
f_w.writerow(['Parent',fitnesseval(valseq,inp),inp])

poplist=initialpop(inp) #poplist contains the parents.
comvalues=[(fitnesseval(par,inp),par) for par in poplist] #Now we create a list comvalues, containing the format (fitness function value (Hamming Dist.), sequence)

i=0  

while i<100:
    print(i)  #Just a counter to see the progress of the code- helpful when we have to generate a large number of sequences.

    comvalues=sorted(comvalues)
    comvalues.insert(randint(0,len(comvalues)-1),(1,poplist[randint(0,len(poplist)-1)])) #We randomly add an individual to poplist to include immigration.
    val=comvalues[randint(0,len(comvalues)-1)][1]  #val picks up one sequence randomly, designated to be one parent.
    oval=comvalues[randint(0,len(comvalues)-1)][1] #oval picks up one sequence randomly, designated to be another parent
    p1=mutation(val)
    p2=mutation(oval) 
    
    n1=(fitnesseval(crossover(p1,p2),inp),crossover(p1,p2)) #n1 contains the crossover sequences in the format of comvalues.
    f_w.writerow([i,n1[0],n1[1]])
   
    n2=(fitnesseval(crossover(p1,p2),inp),crossover(p1,p2))  #n1 contains the crossover sequences in the format of comvalues.
    f_w.writerow([i,n2[0],n2[1]])
    
    comvalues.extend((n1,n2)) #comvalues now contains n1 and n2.
    
    '''Then we remove the first three inidividuals to ensure that the starting population size is maintained after every reproductive event.'''
    comvalues.remove(comvalues[0])
    comvalues.remove(comvalues[0])
    comvalues.remove(comvalues[0])
  
    i+=1
    



    
