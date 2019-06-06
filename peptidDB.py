#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:04:21 2017

@author: hilla
"""
from Bio import SeqIO
import re
#import difflib
#from difflib import SequenceMatcher
import numpy as np
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


## read gene list
File_GeneList='/home/hilla/Aparna/GeneList.txt'
Out_File='/home/hilla/Aparna/Out.fasta'
Input_File='/home/hilla/Aparna/uniprot-human.fasta'

def ReadGeneUniqueList(File_GeneList):
    return [Gene.rstrip('\n') for Gene in open(File_GeneList)]

def FindIsoforms(GeneList,C4):
    return [[k for k in range(len(C4)) if Gene in C4[k]] for Gene in GeneList]

def FindUniquIsoforms (IsoformList,GeneList):
    return [[IsoformList[k],GeneList[k]] for k in range(len(IsoformList)) if len(IsoformList[k])==1]

def WriteUniqueIsoform (Out_File,UniqIso,FastaSeq,Information,UniqGenes):
    fileID=open(Out_File,"w+")
    for k in xrange(len(UniqGenes)):
        Gene=UniqGenes[k]
        ind=UniqIso[k][0]
        s = Information[ind]
        Protein = re.search('sp\\|(.*)\\|', s).group(1)
        Uniport=str(re.search('%s(.*)%s' % (Protein, '_HUMAN'), s).group(1)+'_HUMAN')[1:]
        IsofurmNum='1/1'
        s=s.split(Uniport,1)[1] 
        s=s.replace('OS=Homo sapiens','')
        Description=s.replace(Gene,'')
        fileID.write('>%s\t %s\t %s\t %s\t %s\t\n' % (Gene,Protein,Uniport,IsofurmNum,Description))
        fileID.write('%s\n' % (FastaSeq[ind]))
    fileID.close()

def FindMultipleIsoforms (IsoformList,GeneList):
    return [[IsoformList[k],GeneList[k]] for k in range(len(IsoformList)) if len(IsoformList[k])>1]

#def ConcentrateSeq(MultiIso,MultiGenes,FastaSeq):
#    VSeq= [None] * len(MultiIso)
#    for k in range(len(MultiIso)):
#        Seq=[FastaSeq[iso] for iso in MultiIso[k]]
#        a=Seq[0]
#        for l in range(len(Seq)-1):
#           b=Seq[l+1]
#           f = SequenceMatcher(None, a, b)
#           Point=f.find_longest_match(0, len(a), 0, len(b))
#           f.get_matching_blocks()
#           a=a[Point[0]:Point[2]]
#        VSeq[k]=a

def ConcentrateSeq(MultiIso,MultiGenes,FastaSeq,Information,Out_File):
#    VSeq= [None] * len(MultiIso)
#    for k in range(len(MultiIso)):
#        Seq=[FastaSeq[iso] for iso in MultiIso[k]]
#        Length=min([len(l) for l in Seq])
#        
#        for l in Length:
#            for i in len(Seq):
#                Seq[l,:]
#                
#        for l in range(len(Seq)-1):
#           b=Seq[l+1]
#           f = SequenceMatcher(None, a, b)
#           Point=f.find_longest_match(0, len(a), 0, len(b))
#           f.get_matching_blocks()
#           a=a[Point[0]:Point[2]]
#        VSeq[k]=a
#
#VSeq= [None] * len(MultiIso)
    VSum=np.zeros(len(MultiIso))
    for k in range(len(MultiIso)):
        Seq=[FastaSeq[iso] for iso in MultiIso[k]]
        Length=min([len(l) for l in Seq]) 
        #iterator=[Se[l] for Se in Seq]
        C=[checkEqual2([Se[l] for Se in Seq]) for l in xrange(Length)]
        C=np.array([0>1]+C+[0>1])
        Vind=np.squeeze(np.where(np.logical_not(C)))
        #Vind=np.array(np.where(np.logical_not(C)))
        VV=np.diff(Vind)
        ind=np.argmax(VV)
        Loc=Vind[ind]
        val=VV[ind]
        ind=np.argmax(np.diff(np.where(np.logical_not(C))))
        P1=ind
        P2=val-1+ind
        A=Se[P1:P2]
        VSe=[Se[:P1]+Se[P2:] for Se in Seq]
        #Se[0].find(A, beg=0, end=len(Se[0]))  
        #VInitial=Se.find(A)
        #VInitial=[Se.find(A) for Se in Seq]
        
        VInf=[Information[n] for n in MultiIso[k]] 
        WriteMultipleIsoforms (A,Out_File,MultiIso[k],VInf,VSe,MultiGenes[k],P1) 
#VSum[k]=np.argmax(np.diff(np.where(np.logical_not(C))))

#newstr = oldstr[:4] + oldst[5:]
#VSum.astype(int)
#VV=np.sort(VSum)  
#C=[[checkEqual2[Se[l] for Se in Seq]] for l in xrange(10)]
def WriteMultipleIsoforms (CommonSeq,Out_File,MIso,VInf,VSe,Gene,Point):
    fileID=open(Out_File,"a")
    Num=len(VInf)
    IsoformNum='IsoformNum '+str(Num)
    fileID.write('>%s\t %s\t %s\t\n' % (Gene,IsoformNum,Point))
    fileID.write('%s\n' % (CommonSeq))
    for l in xrange(len(VInf)):
        Isoform=MIso[l]
        s = VInf[l]
        Protein = re.search('sp\\|(.*)\\|', s).group(1)
        Uniport=str(re.search('%s(.*)%s' % (Protein, '_HUMAN'), s).group(1)+'_HUMAN')[1:]
        #str(k)+'/'+str(Num)
        IsofurmNum=str(l+1)+'/'+str(Num)
        s=s.split(Uniport,1)[1] 
        s=s.replace('OS=Homo sapiens','')
        Description=s.replace(Gene,'')
        fileID.write('>%s\t %s\t %s\t %s\t %s\t %s\t\n' % (Gene,Protein,Uniport,IsofurmNum,Point,Description))
        fileID.write('%s\n' % (VSe[l]))
    fileID.close()

def checkEqual2(iterator):
   return len(set(iterator)) <= 1

#####  Run
FastaSeq,Information=ReadInput(Input_File)
GeneList=ReadGeneUniqueList(File_GeneList)
IsoformList=FindIsoforms(GeneList,C4)
UniqIso,UniqGenes=zip(*FindUniquIsoforms(IsoformList,GeneList))
WriteUniqueIsoform (Out_File,UniqIso,FastaSeq,Information,UniqGenes)
MultiIso,MultiGenes=zip(*FindMultipleIsoforms (IsoformList,GeneList))
ConcentrateSeq(MultiIso,MultiGenes,FastaSeq,Information,Out_File)
#####
from difflib import SequenceMatcher
a = "Now is the time for all good men to come to the aid of their party."
b = "All good men are created equal."
s = SequenceMatcher(None, a, b)
s.find_longest_match(0, len(a), 0, len(b))


ll=np.array([len(seq) for seq in VSeq]) 

a0=Seq[0]
a1=Seq[1]
a2=Seq[2]
a3=Seq[3]
a4=Seq[4]
a5=Seq[5]


l=0
a=Seq[0]
b=Seq[l+1]
f = SequenceMatcher(lambda x: x == " ",a,b)

f = SequenceMatcher(None, a, b)
Point=f.find_longest_match(0, len(a), 0, len(b))
f.get_matching_blocks()
a=a[0:380]
VSeq[k]=a
ll=np.array([len(k) for k in MultiIso])       
a=Seq[0]
b=Seq[1]
f = SequenceMatcher(None, a, b)
#f = SequenceMatcher(None, Seq[0], Seq[1])
f.find_longest_match(0, len(a), 0, len(b))
Match(a=0, b=4, size=5)
f.get_matching_blocks()

LL=np.array([len(l) for l in MultiIso])
s = Information[0]
start = 'GN=YWHAB' 
end = 'SV=3'

start = 'sp' 
end = '|1433B'
s=C1[0]
start = 'sp\\|' 
end = '\\|'


s = Information[0]
result = re.search('%s(.*)%s' % (start, end), s).group(1)
result2=re.search('%s(.*)%s' % (result, '_HUMAN'), s).group(1)
"isoform"
'1/1'
ss=s.replace('OS=Homo sapiens','')
ss=s.replace(UniqGenes[k],'')
print(result)
print(result).group(1)



import re
s=Information[0]
s=C1[0]

result = re.search('sp|(.*)|', s)
s = 'asdf=5;iwantthis123jasd'
result = re.search('asdf=5;(.*)123jasd', s)
print result.group(1)


fileID = open('/home/hilla/Galagan/SendB2/Motifs.txt','w')
fileID.write("Logo picture number\t Number of transcription factors\t Number of query sequences\t  List of e-values\t List of transcription factors\t List of query sequences\n" )
for i in range(len(NumTF)):
    fileID.write("%s\t %s\t %s\t" % (VTF[i][0], NumTF[i], NumQuery[i]))
    fileID.write("%s\t"  % (','.join(map(str,Eval[i]))))
    fileID.write("%s\t"  % (','.join(map(str,VTF[i]))))
    fileID.write("%s\n"  % (','.join(map(str,VQuery[i]))))
fileID.close()
zip(*G)
L1,L2=zip(*[[IsoformList[k],GeneList[k]] for k in range(len(IsoformList)) if len(IsoformList[k])==1])

L=[[k for k in range(len(C4)) if Gene in C4[k]] for Gene in GeneList]

L=[[k for k in range(10)] for i in range(10) if i<5]
for Gene in GeneList:
    for k in range(len(C4)):
        if Gene in C4[k]:
            print k
GeneList=[Gene.rstrip('\n') for Gene in open(File_GeneList)]
for Gene in GeneList:
    print(Gene)
    
def ReadInput(Input_File):
    #C1=[seq_record.id for seq_record in SeqIO.parse(input_file, "fasta")]
    FastaSeq=[seq_record.seq.tostring() for seq_record in SeqIO.parse(input_file, "fasta")]
    Information=[seq_record.description for seq_record in SeqIO.parse(input_file, "fasta")]
    return [FastaSeq,Information]

a,b=ReadInput(Input_File)
Input_File='/home/hilla/Aparna/uniprot-human.fasta'
input_file='/home/hilla/Aparna/uniprot-human.fasta'
G = [[int(item) for item in line.split()] for line in open(file1)]
C1=[seq_record.id for seq_record in SeqIO.parse(input_file, "fasta")]
FastaSeq=[seq_record.seq.tostring() for seq_record in SeqIO.parse(input_file, "fasta")]
#C2=[repr(seq_record.seq) for seq_record in SeqIO.parse(input_file, "fasta")]
Information=[seq_record.description for seq_record in SeqIO.parse(input_file, "fasta")]
/
seq_record.description.split()
C4=[seq_record.description.split() for seq_record in SeqIO.parse(input_file, "fasta")]
for 
for seq_record in SeqIO.parse(input_file, "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

from Bio import SeqIO
input_file='/home/hilla/Aparna/uniprot-human.fasta'
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
fasta_sequences.seq.tostring()
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)