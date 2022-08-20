import random
from Bio.SeqRecord import SeqRecord
import streamlit as st
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq 
from Bio import SeqIO
#from collections import Counter
#import neatbio.sequtils as utils
import numpy as np 
#from PIL import Image
import joblib
import pandas as pd
#from bokeh.plotting import figure
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from PIL import Image
from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab


def Prediction_OfDiabetes():
   @st.cache(allow_output_mutation=True)
   def load(scaler_path, model_path):
       sc = joblib.load(scaler_path)
       model = joblib.load(model_path)
       return sc , model

   def inference(row, scaler, model, feat_cols):
       df = pd.DataFrame([row], columns = feat_cols)
       x= scaler.transform(df)
       feature=pd.DataFrame(x,columns=feat_cols)
       if model.predict(feature)==0:
        return "This is a healthy person!"
       else: return "This person has high chances of having diabetics!"
   age =           st.number_input("Age in Years", 1, 150, 25, 1)
   pregnancies =   st.number_input("Number of Pregnancies", 0, 20, 0, 1)
   glucose =       st.slider("Glucose Level", 0, 200, 25, 1)
   skinthickness = st.slider("Skin Thickness", 0, 99, 20, 1)
   bloodpressure = st.slider('Blood Pressure', 0, 122, 69, 1)
   insulin =       st.slider("Insulin", 0, 846, 79, 1)
   bmi =           st.slider("BMI", 0.0, 67.1, 31.4, 0.1)
   dpf =           st.slider("Diabetics Pedigree Function", 0.000, 2.420, 0.471, 0.001)

   row = [pregnancies, glucose, bloodpressure, skinthickness, insulin, bmi, dpf, age]

   if (st.button('Find Health Status')):
       feat_cols = ['Pregnancies', 'Glucose', 'BloodPressure', 'SkinThickness', 'Insulin', 'BMI', 'DiabetesPedigreeFunction', 'Age']
       sc, model = load('models/scaler.joblib', 'models/model.joblib')
       result=inference(row, sc, model, feat_cols)
       st.write(result)
        
def Read_CalculateGc():
    seq_file = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
    
    seq_file2 = st.file_uploader("Upload Pig FASTA File", type = ["fasta"])
   
    str_Pig_dna_record = ""
    str_Human_dna_record=""
    if seq_file is not None:  
        Human_dna_record = seq_file.read()           
        st.write(Human_dna_record)
        str_Human_dna_record = str(Human_dna_record) 
        if seq_file2 is not None:
            Pig_dna_record = seq_file2.read()           
            st.write(Pig_dna_record)
            str_Pig_dna_record = str(Pig_dna_record)
            
            
    if seq_file is  None:          
        if seq_file2 is not None:
            Pig_dna_record = seq_file2.read()           
            st.write(Pig_dna_record)
            str_Pig_dna_record = str(Pig_dna_record)

            
    if st.button("Calc GC_Content"):
            if(len(str_Human_dna_record)>0):
                Human_gc_content = (str_Human_dna_record.count("C") + str_Human_dna_record.count("G")) / len(str_Human_dna_record) * 100
                st.write('Human GC content =', Human_gc_content)
            if len(str_Pig_dna_record)> 0:           
                Pig_gc_content = (str_Pig_dna_record.count("C") + str_Pig_dna_record.count("G")) / len(str_Pig_dna_record) * 100
                st.write( 'Pig GC content =', Pig_gc_content)
                
def Visualization():
    seq_file = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
    
   
    str_Pig_dna_record = ""
    if seq_file is not None:  
        Human_dna_record = seq_file.read()           
        st.write(Human_dna_record)
        str_Human_dna_record = str(Human_dna_record) 
            
        if st.button("Visualize"): 
             Human_dna_record = pd.DataFrame({"c":[str_Human_dna_record.count("C")],

                                              "g":[str_Human_dna_record.count("G")],
                         "t":[str_Human_dna_record.count("T")],
                         "a":[str_Human_dna_record.count("A")],
                         })
             
            
             fig, ax = plt.subplots()
             ax.hist(Human_dna_record, bins=20)
             ax.legend(['A','C','G','T'] ,loc='upper right')
             st.pyplot(fig)



def Alignment():
    seq_file = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
    
    seq_file2 = st.file_uploader("Upload Pig FASTA File", type = ["fasta"])
   
    c=0
    str_Human_dna_record=""
    str_Pig_dna_record=""
    if seq_file is not None and seq_file2 is not None:  
        Human_dna_record = seq_file.read()   
        Pig_dna_record = seq_file2.read()      
        str_Human_dna_record = str(Human_dna_record) 
        str_Pig_dna_record = str(Pig_dna_record)
        
    if st.button("Align"): 
        if(seq_file is None or seq_file2 is None):
            st.write("if you Don't Have the File Please Select the Generate Button from Pannel and Re-Alignment..")   
              #if(st.button("Upload Again")):
                 # Alignment()
        else:
            alignments = pairwise2.align.globalxx(str_Human_dna_record, str_Pig_dna_record)
            for a in alignments:
                c=0
                file1 = open("Alignment.txt", "w") 
                file1.write(format_alignment(*a)[1:])
                file1.close()  
                c=c+1
            
        if(c>0):
            with open("Alignment.txt") as f:
                    st.download_button("Download the Alignment file", data=f, file_name=("Alignment.txt"))

def Generate():
    newfh = open('randomseqs.fasta','w')
    seq_num =st.number_input("Enter Number of Seq", 1, 1000, 500, 1)
    st.text("Select The number of nitrogenous bases")
    From=st.number_input("From", 1, 10000, 4000, 1)
    To =st.number_input("To", 1, 20000, 15000, 1)
    if st.button("Generate"):
        if(From < 80 and To<=200 ):
            st.write("You Must Enter Number range from 4000 , 15000")
        else:
            for i in range(1,seq_num):
                # Creates a random number in the range of 4000-15000
                rsl = random.randint(From,To)
                # Generate the random sequence
                rawseq = new_rnd_seq(rsl)
                # Generate a correlative name
                seqname = 'Sequence_number_' + str(i)
                rec = SeqRecord(Seq(rawseq),id=seqname,description='')
                SeqIO.write([rec],newfh,'fasta')
            with open("randomseqs.fasta") as ran:
                st.download_button("Download the Sequence file", data=ran, file_name=("randomseqs.fasta"))
   
def Translated():
    table = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                  '---': '-',
                  }
    seq_file = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
    seq=""
    if(seq_file  is not None):
       seq_file.readline()
       seq = Seq(seq_file.read())
       seq= seq.upper()
       seq = seq.replace("\n", "")
       seq = seq.replace("\r", "")
       seq = seq.replace("N", "")
       seq = seq.replace(">", "")
       
    protein =""
    if (len(seq)%3 == 0):
        for i in range(0, len(seq), 3):
                  codon = seq[i:i + 3]
                  protein+= table[codon]
    elif (len(seq)%3 != 0):
      for i in range(0, len(seq)-2, 3):
              codon = seq[i:i + 3]
              if(codon not in table):
                  pass
              else:
                  protein+= table[codon]

    if st.button("Click Here"):
        chunk_seq("Insulin",protein,20,20)


def filter():
  handle = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
  #handle = open(seq_file, "rU")
  if(handle is not None):
      handle.readline()
      seq = Seq(handle.read())
      if seq.count('N') > 0 :
          seq=seq.replace('N', "") 
  if st.button("Remove"):
    output_handle = open("N_removed.fasta", "w")
    output_handle.write(str(seq))
    output_handle.close()
    with open("N_removed.fasta") as ran:
      st.download_button("Download the Sequence file", data=ran, file_name=("N_removed.fasta"))



def chunk_seq(seq_name, sequence, chunk_size, chunk_increment):
        newfh= open("output.fasta",'w')
        index = 0
        chunks = []
        seqname=""
        while (index < len(sequence)):
            chunk_seq = sequence[index:index+chunk_size]
            if len(chunk_seq) != chunk_size:
                break
            chunks.append(('_'.join([seq_name, str(index), str(index+chunk_size)]), chunk_seq))
            index += chunk_increment
        newfh.write(str(chunks))
        newfh.close()
        with open("output.fasta") as ran:
            st.download_button("Download the Sequence file", data=ran, file_name=("output.fasta"))

                
# Generate Random Sequence
def new_rnd_seq(sl):
    """ Generate a random DNA sequence with a sequence length of "sl" (int). """
    s = ''
    for x in range(sl):
        s += random.choice('ATCG')
        # s += random.sample(’ATCG’,1)[0] is not so fast.
    return s

    
def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    #xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    #yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()
    st.set_option('deprecation.showPyplotGlobalUse', False)

def complement():
    comp=""
    f = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
    
    if(f is not None):
        f.readline()
        seq = Seq(f.read())
        comp=seq.complement()   
        comp2=seq.reverse_complement()

   
    if st.button("Complement"):
        newfh= open("outputComplent.fasta",'w')
        newfh.write(str(comp))
        newfh.close()
        with open("outputComplent.fasta") as ran:
          st.download_button("Download the Sequence file", data=ran, file_name=("outputComplent.fasta"))
  
    if st.button("Reverse"):
        newfh= open("outputReverse.fasta",'w')
        newfh.write(str(comp2))
        newfh.close()
        with open("outputReverse.fasta") as ran:
          st.download_button("Download the Sequence file", data=ran, file_name=("outputReverse.fasta"))


def Txt_to_Fasta():
    seq=""
    f = st.file_uploader("Upload Human FASTA File", type = ["txt"])
    if(f is not None):
        seq = (f.read())
        
    if st.button("Convert"):
      newfh= open("outputFasta.fasta",'w')
      newfh.write(str(seq)[1:])
      newfh.close()
      with open("outputFasta.fasta") as ran:
        st.download_button("Download the Sequence file", data=ran, file_name=("outputFasta.fasta"))


def About():
    st.header("Team : ")
    
    st.subheader("Mohamed Mamdouh Allam Ahmed")
    st.text("ID = 1620195237 ")
    st.text("Section (1), Group (4)")
    
    st.subheader("Mahmoud Sayed Youssef Kotb")
    st.text("ID = 1620195241 ")
    st.text("Section (1), Group (4)")
    
    st.subheader("Mena Nashat Fayez")
    st.text("ID = 1620195277 ")
    st.text("Section (3), Group (4)")
    
    st.subheader("Doaa Sayed Ibrahim Morsy")
    st.text("ID = 1620195096 ")
    st.text("Section (1), Group (2)")
    
    st.subheader("Rahma Yasser Mahmoud")
    st.text("ID = 1620195110 ")
    st.text("Section (2), Group (2)")
    
    st.subheader("Yasmeen Hossam ELdin")
    st.text("ID = 1620195237 ")
    st.text("Section (4), Group (4)")


def main():
    st.title("Simple Bioinformatics App")
    
    menu = ["Inrto", "Sequence Generator" ,"Translated Protein","Complement","Convert to Fasta","Filter","Calculate GC","Alignment","Visualization Plot", "DotPlot", "Diabetes Prediction", "About"]
    choice = st.sidebar.selectbox("Select Activity", menu)
    
    if choice == "Inrto":
        st.subheader("Intro")
        st.write("In this We have tried to predict the probability of a person having diabetes based on some data fields...And do some operations on the Fasta files")
        image = Image.open('C:\\Users\\Electronica Care\\Desktop\\Python prject 2\\static\\Images\\14.jpg')
        st.image(image, use_column_width=True)
    elif choice == "Calculate GC":
        st.subheader("GC Content")
        Read_CalculateGc()
    elif choice == "Filter":
         st.subheader("Filter DNA File")
         filter()
    elif choice == "Convert to Fasta":
       st.subheader("Convert from txt file to FASTA file")
       Txt_to_Fasta()
    elif choice == "Visualization Plot":
        st.subheader("Visualization of Human and Pig Insulin")
        Visualization()
    elif choice == "Alignment":
        st.subheader("Alignment between Human and Pig Insulin")
        Alignment()
    elif choice == "Complement":
        st.subheader("Complement of DNA")
        complement()
    elif choice == "Translated Protein":
        st.subheader("Translated Protein of DNA Sequence ")
        Translated()

    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file = st.file_uploader("Upload Human FASTA File", type = ["fasta"])
        seq_file2 = st.file_uploader("Upload Pig FASTA File", type = ["fasta"])
        
        str_Pig_dna_record = ""
        if seq_file and seq_file2 is not None:  
            Human_dna_record = seq_file.read()     
            Pig_dna_record = seq_file2.read() 
            
            str_Human_dna_record = str(Human_dna_record) 
            str_Pig_dna_record = str(Pig_dna_record)
            
            cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                dotplotx(str_Human_dna_record[0:cus_limit],str_Pig_dna_record[0:cus_limit])  
                st.pyplot()
            

    elif choice == "About":
        st.subheader("About Team")
        image = Image.open('C:\\Users\\Electronica Care\\Desktop\\Python prject 2\\static\\Images\\8.png')
        st.image(image,width=250, use_column_width= 'auto' )
        About()
    elif choice == "Diabetes Prediction":
         st.subheader("Prediction of Diabetes in Human")
         image = Image.open('C:\\Users\\Electronica Care\\Desktop\\Python prject 2\\static\\Images\\1_INSggrGiQ1lCgU8YTsfEVw.png')
         st.image(image,width=500, use_column_width= 'auto' )
         Prediction_OfDiabetes()
    elif choice == "Sequence Generator":
        st.subheader("You Can Generate Sequence Here.")
        Generate()
   
        
if __name__ == '__main__':
    main()