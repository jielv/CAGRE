# GRE
Genomic region enrichment tool 

=================
=================

Jie Lv  
jlu@houstonmethodist.org  
Chen Lab (PI: Dr. Kaifu Chen)  
Center for Bioinformatics and Computational Biology Department of Cardiovascular Sciences  
The Methodist Hospital Research Institute Department of Cardiothoracic Surgery  
Weill Cornell Medical College, Cornell University  

# Overview
Identifying signatures of a group of genes or transcripts associated with certein development and disease process on a genome-wide scale is of great biological interests. Thanks to the fast development of next-generation sequencing assays, we have a ever-increasing list of  large-scale biological signature reference databases, in the format of a set of genomic regions (e.g., TF binding sites, DNA or histone modification positions, etc.). Here, we present the Genomic Region Enrichment(GRE) tool, a simple, flexible commandline tool that provides biologist the ability to perferm comparative analyses with signatures of their lists of gene/transcript. 


RUNNING
============

Step 1 - Perform genomic regions enrichment analysis for groups of genes
------------------------
gene_enrich_v1.0.py 

DESCRIPTION:  
  read mapping file from genomic regions to genes and gene list of interet, and perform comparative analyses with the enriched genomic regions for different gene groups.

INPUT:  

```
  -h, --help            show this help message and exit  
  -x FILE, --xls=FILE   danpos output showing locations of selected peaks
                        (center, column 4 ) and its related transcripts  
                        (column 9) in txt format  
  -g FILE, --gene=FILE  gene annotation file in xls format from ucsc,
                        containint mapping between transcipt id  (col 1) to
                        gene symbol (col 12)  
  -s FILE, --symbol=FILE 
                        mapping between gene names from different sources:
                        e.g. ucsc gene symbol(col 1) to official gene symbol
                        (last col), acounting for the situation that xls and
                        gene list have different name system.  
  -u FILE, --up=FILE    list of first group of genes, e.g., upregulated genes
                        after some treatment  
  -U FILE, --down=FILE  list of second group of genes  
  -c FILE, --control=FILE  
                        list of control genes  
  -S SIZE, --size=SIZE  the number of control genes selected from all control
                        genes, default is None and equals the number of first
                        group of genes  
  -n BREAK_NUMBER, --number=BREAK_NUMBER  
                        number of breaks points for distance between peak and
                        TSS  
  -w BREAK_STEP, --step=BREAK_STEP  
                        step size of break points for distance between peak
                        and TSS, the total range examined would be
                        [0,break_number*break_step]     
```

OUTPUT:  
  cumulative genomic region enrichment stored in txt file

EXECUTING EXAMPLE: 
```
python gene_enrich_v1.0.py -x regionsToGenes.xls -g  mm9.20150218.knownGene.xls  -u up_genes.txt -U down_genes.txt -c control_genes.txt  -n 100 -w 10
```

Step 2 - Perform visualization
-------------------------
GRE_visual_v1.0.R

DESCRIPTION:  
Given a output txt file of cumulative genomic region enrichment for groups of genes of interest,  print out the cumulative figure in pdf format.

INPUT:
```
    args[1]: txt file for cumulative genomic region enrichment  
    args[2]: working directory
```

OUTPUT:  
  Figures of cumulative enrichment number and P-value across a given range for distance between TSS to genomic regions in pdf format.

EXECUTING EXAMPLE:
```

Rscript GRE_visual_v1.0.R up_genes.txt.regionsToGenes.xls.100.10.1000.cumulative.txt
```

