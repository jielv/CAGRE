#!/usr/bin/python

from collections import defaultdict
import sys
from optparse import OptionParser
import numpy as np
import scipy.stats
from time import time
import random
import os

parser = OptionParser()

parser.add_option("-x", "--xls", type="string", dest="xls", metavar="FILE",help="danpos output showing locations of selected peaks (center, column 4 ) and its related transcripts (column 9) in txt format" ) 
parser.add_option("-g", "--gene", type="string", dest="gene", metavar="FILE",help="gene annotation file in xls format from ucsc, containint mapping between transcipt id  (col 1) to gene symbol (col 12)" ) 
parser.add_option("-s", "--symbol", type="string", dest="symbol", metavar="FILE",help="mapping between gene names from different sources: e.g. ucsc gene symbol(col 1) to official gene symbol (last col), acounting for the situation that xls and gene list have different name system.") 
parser.add_option("-u", "--up", type="string", dest="up_symbol", metavar="FILE",help="list of first group of genes, e.g., upregulated genes after some treatment") 
parser.add_option("-U", "--down", type="string", dest="down_symbol", metavar="FILE",help="list of second group of genes") 
parser.add_option("-c", "--control", type="string", dest="control_symbol", metavar="FILE",help="list of control genes") 
parser.add_option("-S", "--size", type="int", dest="size", default=None, help="the number of control genes selected from all control genes, default is None and equals the number of first group of genes ")
parser.add_option("-n", "--number", type="int", dest="break_number", default=40, help="number of breaks points for distance between peak and TSS" )
parser.add_option("-w", "--step", type="int", dest="break_step", default=500, help="step size of break points for distance between peak and TSS, the total range examined would be [0,break_number*break_step]" )


(options, args) = parser.parse_args()


def control_genes(file,size,symbol2dis):
    """
    select number of control genes for P-value conparisons.
    """
    if file:
        f=open(file)
    list=[]
    while True:
        line=f.readline()
        if not line:
            break
        symbol=line.strip()
        list.append(symbol)
    random.seed(1)
    controls=random.sample(list,size)
    list=[]
    n=0
    for symbol in controls:
        if symbol in symbol2dis:
            n+=1
            list.append(symbol2dis[symbol])
        else:
            list.append([])
    print >> sys.stderr," ".join((file, 'number of input symbols: '+str(len(list))+";", 'number of symbols with binding peaks: '+str(n)+";",'number of tss-peak distance for the set of symbols: '+str(len([item for sublist in list for item in sublist]))+";"))
    #print >> sys.stderr, file, 'number of control symbols:',size,len(list), 'number of symbols with binding peaks', n,'number of tss-peak distance for the set of symbols:', len([item for sublist in list for item in sublist])
    return size, list
            
def load_tx2gene(file):
    """
    load ucsc xls file and store mapping between transcript id to gene id
    """
    start_time=time()
    if file:
        f=open(file)
    tx2gene=defaultdict(str)
    while True:
        line=f.readline()
        if line.startswith(('#','track','browser')):continue
        if not line:
            break
        cols=line.strip().split("\t")
        if cols[1]=="start":
            continue
        tx=cols[0]
        gene=cols[-1]
        tx2gene[tx]=gene
    #print >> sys.stderr, 'time cost for load_tx2gene:', time()-start_time
    return tx2gene

def load_gene2symbol(file):
    """
    if given gene list if not official gene symbol, load the mapping information between the given gene names to official gene symbol
    """
    start_time=time()
    if file:
        f=open(file)
    gene2symbol=defaultdict(str)
    while True:
        line=f.readline()
        if line.startswith(('#','Input')):continue
        if not line:
            break
        cols=line.strip().split("\t")
        if cols[-1]=="-":
            continue
        gene=cols[0]
        symbol=cols[-1].split(";")[0]
        gene2symbol[gene]=symbol
    #print >> sys.stderr, 'time cost for load_gene2symbol:', time()-start_time
    return gene2symbol
def load_tx2symbol(tx2gene,gene2symbol):
    """
    map transcript id to the gene names if not official gene symbol
    """
    tx2symbol=defaultdict(str)
    for tx in tx2gene:
        gene=tx2gene[tx]
        if gene not in gene2symbol:
            continue
        else:
            tx2symbol[tx]=gene2symbol[gene]
    return tx2symbol
def loadxls(file):
    """
    load mapping information between binidng peak to TSS
    """
    start_time=time()
    if file:
        f=open(file)
    tx2peak_dis=defaultdict(list)
    tx2dis=defaultdict(list)
    while True:
        line=f.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        if cols[-1]=="relatedGenes":
            continue
        gene_set=[]
        chr=cols[0]
        start=int(cols[1])
        end=int(cols[2])
        center=start+int((end-start)/2)
        genes=cols[-1].split(':')[1].split(',')
        for gene in genes:
            id,TSS,TTS=gene.split('/')
            TSS,TTS = map (int, (TSS,TTS))
            dis=center-TSS
            gene_set.append([id,TSS,TTS,dis,abs(dis)])
        gene_set.sort(key=lambda x:x[4])
        #cutoff=max(int(gene_set[0][4]*range),gene_set[0][4]+distance)
        #for each peak, only map it to the closest tx
        gene=gene_set[0][0]
        dis=gene_set[0][4]
        tx2peak_dis[gene].append(dis)
    for tx in tx2peak_dis:
        tx2peak_dis[tx].sort()
        tx2dis[tx]=tx2peak_dis[tx]
    return tx2dis    

def setStep(list,break_number,step):
    breaks=[]
    N=max(list)
    for i in range(break_number+1):
        if step*i > N:
            print >> sys.stderr, "wrong, range of TSS distance is too large", step*i, N
            print >> sys.stderr, "n*w must equal or less than:", N
            sys.exit()
        breaks.append(step*i)
    return breaks
        
def loadList(file,symbol2dis):
    if file:
        f=open(file)
    list=[]
    n=0
    m=0
    while True:
        line=f.readline()
        if not line:
            break
        m+=1
        symbol=line.strip()
        if symbol in symbol2dis:
            n+=1
            list.append(symbol2dis[symbol])
        else:
            list.append([])
    print >> sys.stderr," ".join((file, 'number of input symbols: '+str(len(list))+";", 'number of symbols with binding peaks: '+str(n)+";",'number of tss-peak distance for the set of symbols: '+str(len([item for sublist in list for item in sublist]))+";"))
    return m,list
def listFlat(list):
    list_flated=[item for sublist in list for item in sublist]
    return list_flated
def mannWhitenyU(l1,l2):
    try:
        s,p=scipy.stats.mannwhitneyu(l1,l2,alternative='greater')
    except:
        p=1
    return p
def cumCount(list_all,N_all,list_control,N_control,list_up,N_up,list_down,N_down,break_number,step,output):
    start_time=time()
    of=open(output,"w")
    list_up_flat=listFlat(list_up)
    list_down_flat=listFlat(list_down)
    list_all_flat=listFlat(list_all)
    list_control_flat=listFlat(list_control)
    breaks=setStep(list_all_flat,break_number,step)
    for point in breaks:
        
        up_vec=[sum(i<=point for i in l) for l in list_up]
        down_vec=[sum(i<=point for i in l) for l in list_down]
        all_vec=[sum(i<=point for i in l) for l in list_all]
        control_vec=[sum(i<=point for i in l) for l in list_control]
        
        P_up_all=mannWhitenyU(up_vec,all_vec)
        P_down_all=mannWhitenyU(down_vec,all_vec)
        P_control_all=mannWhitenyU(control_vec,all_vec)
        P_up_down=mannWhitenyU(up_vec,down_vec)
        
        count_up=sum(i<=point for i in list_up_flat)/float(N_up)*100
        count_down=sum(i<=point for i in list_down_flat)/float(N_down)*100
        count_all=sum(i<=point for i in list_all_flat)/float(N_all)*100
        count_control=sum(i<=point for i in list_control_flat)/float(N_control)*100
        
        of.write("\t".join(map(str,[point,count_up,count_down,count_control,count_all,P_up_all,P_down_all,P_control_all,P_up_down]))+"\n")
    of.close()

if __name__ == "__main__":
    start_time=time()
    tx2gene=load_tx2gene(options.gene)
    if options.symbol:
        gene2symbol=load_gene2symbol(options.symbol)
        tx2symbol=load_tx2symbol(tx2gene,gene2symbol)
    else:
        gene2symbol={}
        tx2symbol=tx2gene


    tx2dis=loadxls(options.xls)

    print >> sys.stderr, " ".join(("total transcript: "+str(len(tx2gene))+";", "number of genes from source 1 mappable to that from source 2: "+str(len(gene2symbol))+";","number of transcripts mappable to symbols: "+str(len(tx2symbol))+";", "number of transcripts with binding peaks: "+str(len(tx2dis))+";"))

    symbol2dis=defaultdict(list)

    gene2dis=defaultdict(list)
    for tx in tx2dis:
        if tx not in tx2symbol:
            continue
        else:
        # for a particular gene symbol, put all TSS-peak distances into a pool, 
            symbol2dis[tx2symbol[tx]].extend(tx2dis[tx])


    #for a group of gene symbols, put all the tss-peak distance in a single pool
    list_all=[]
    if options.symbol:
        N_all=len(gene2symbol)
    else:
        N_all=len(set(tx2symbol.values()))
    for symbol in symbol2dis:
        list_all+=symbol2dis[symbol]
    print >> sys.stderr," ".join (('number of total symbols: '+str(N_all)+";", 'number of total symbols with peaks: '+str(len(symbol2dis))+";", 'number of tss-peak distance for all symbols: '+str(len(list_all))+";"))

    N_up,list_up=loadList(options.up_symbol, symbol2dis)
    N_down,list_down=loadList(options.down_symbol, symbol2dis)
    N_control,list_control=loadList(options.control_symbol, symbol2dis)
    if not options.size:
        control_size=N_up
    else:
        control_size=options.size

    N_control_selected,list_control_selected=control_genes(options.control_symbol,control_size, symbol2dis)

    file_head=os.path.basename(options.up_symbol)
    file_head2=os.path.basename(options.xls)
    output=file_head+"."+file_head2+"."+str(options.break_number)+"."+str(options.break_step)+"."+str(control_size)+".cumulative.txt"
    cumCount(list_control, N_control,list_control_selected, N_control_selected,list_up,N_up,list_down,N_down,options.break_number,options.break_step,output)
    print >> sys.stderr, 'time cost for the whole process:', time()-start_time
