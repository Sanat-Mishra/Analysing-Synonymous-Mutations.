#Obtaining Translation Initiation rate (TIR) from the web.

'''The TIR values were available from an online service at "https://sbi.postech.ac.kr/utr_designer/". 
On inputting the UTR sequence, coding sequence and 16s rRNA sequence (is common for an organism), 
the TIR values were retrieved. Since we had to perform this task for hundreds of sequences, 
automating it was necessary.
The TIR values were then written onto a csv file. 

NOTE- 1) While running the code, it is best to leave the system undisturbed, since not only will it 
         interrupt whatever you plan on doing, there are chances you might accidentaly close the browser 
         this uses and terminate the code.
      2) Right click on the respective search bar/text to access the page source and see
         the HTML code for it, this is done multiple times.
      3) Install Selenium on your system.  
      4) The Firefox driver can be installed from 'https://github.com/mozilla/geckodriver/releases' .
'''

#Ensure Selenium is installed and the correct driver has been downloaded.

from selenium import webdriver   #Required module, selenium easily automates form-submission and data retrieval.
from selenium.webdriver.common.keys import Keys
import csv
import math


'''The following few lines are not necessary and results maybe retrieved in any manner you decide.'''

h=[1,2,5,10,50]  #list of Hamming Distances
genes=['Arca','CRP','FnR','Fur','hns','rob','soxS']  #list of genes 
gene_exp=csv.writer(open('/Users/sanatmishra27/Desktop/Plots/SynonymousSeq/Global/'+genes[0]+'/'+'Ham'+str(h[0])+'-'+genes[0]+'GeneExp.csv','w'))   #The path of the file onto which the TIR is written.
gene_exp.writerow(["TIRValues"])
dir='/Users/sanatmishra27/Desktop/Plots/SynonymousSeq/Global/'+genes[0]+str('/Ham')+str(h[0])+'-'+genes[0]+str('.csv')   #The path of the file containing the sequences- which will be read.
i=1

'''From here on the actual implementation begins. '''

with open(dir,'r') as readFile:  #opens the file containing sequences (variable dir) in read mode.
    reader=csv.reader(readFile)  #initialising the reader
    next(reader)  #skipping the first row, since it contains the column heads/titles ('Sequences',etc.)
    for line in reader:

        print(i)  #Counter to see progress of code.
       
        driver = webdriver.Firefox()  #I have used the driver for Firefox, similar ones are available for Chrome/IE.
        driver.implicitly_wait(45)    #Wait for 45 sec, until the page has loaded.
        driver.maximize_window()      #Maximise the browser.
        driver.get("https://sbi.postech.ac.kr/utr_designer/")   #The website we wish to access.

        search_fieldUTR = driver.find_element_by_id("pre")      #The required fields' HTML id was used to fill it up.
        search_fieldPCS = driver.find_element_by_id("gene35")   #Similar to above step.
        search_fieldUTR.clear() #Ensuring it is clear
        search_fieldPCS.clear() #Ensuring it is clear

        a='TTAGTTGGCAATTTAGGTAGCAAAC'          #The UTR for the gene is assigned variable a.
        b=str(line[0]).replace('U','T')[:-3]   #Barring the stop codon, the rest of the sequence is assigned variable 'b'.
     
        search_fieldUTR.send_keys(a)           #Input variable a
        search_fieldPCS.send_keys(b)           #Input variable b
        
        '''Now we submit the variables into the form.'''
        
        search_fieldUTR.submit()                
        search_fieldPCS.submit()
        
        '''Then we access the value we need from the results page.
           The exact location of the result can be viewed in HTML by 
           inspecting source (right click on required text).'''
        
        data=driver.find_element_by_xpath("/html/body/table/tbody/tr/td/table/tbody/tr[2]/td[11]").text
        
        gene_exp.writerow([(2500*math.exp(-0.45*float(data)))])   #This formula was used to convert the value to TIR.
  
        driver.quit()  #Quit the browser.
        i+=1
