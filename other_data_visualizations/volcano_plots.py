import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import os

'''
This script's function:
    - To make volcano plots from the repac or DGE data
    - generate a list of the shortend 3' UTR genes and lengthened 3' UTR genes
    for a given sample into a csv file
'''

#This function does two things:
    #1. get the volcano graph based on the inputted arguments
    #2. get a list of the shortend 3' UTR genes and lengthened 3' UTR genes
def get_another_graph(filename, col_names, title, control_str, fc_thres, out_path, in_path, labels):
    
    #ARGUMENT DESCRIPTIONS:
        #filename: name of the repac output file
        #col_names: names of the columns you want to extract from the repac dataframe (gene, cFC, adj p-val)
        #title: What you want the title of the plot to be
        #control_str: Tells you what is the control repac is compared to. Ex: Uninfected, Mock
        #fc_thres: What cFC threshold you want, impacts # of genes that are considered to have changed 3' UTR length
        #out_path: The location in your computer files that you want the plot to save to
        #in_path: the path to the folder that contains the filename
        #labels: a list of your axes labels (fist element: x label, second element: y label)
        
    # reads in the repac output file as a dataframe (like a table in R)
    path = in_path+filename
    df = pd.read_csv(path)
    
    plt.scatter(x=df[col_names[1]],y=df[col_names[2]].apply(lambda x:-np.log10(x)),s=1,label="Not significant",color='black')

    # highlight down- or up- regulated genes
    #gets a filtered df containing only upregulated or downregulated genes
    down = df[(df[col_names[1]]<=-fc_thres)&(df[col_names[2]]<=0.05)]
    up = df[(df[col_names[1]]>=fc_thres)&(df[col_names[2]]<=0.05)]
    
    #plots these upregulated and downregulated dataframes
    # x is cFC and y is -log10(adj p-val)
    plt.scatter(x=down[col_names[1]],y=down[col_names[2]].apply(lambda x:-np.log10(x)),s=3,label="Shortened n="+str(len(down)),color='xkcd:ocean blue')
    plt.scatter(x=up[col_names[1]],y=up[col_names[2]].apply(lambda x:-np.log10(x)),s=3,label="Lengthened n="+str(len(up)),color='xkcd:coral')
    
    #IF YOU NEED THIS...
    #this gets the lengthened and shortened genes, their gene name, adj p-value, and cFC, and whether they are 
    #lengthened or shortened.
    #and puts them in a dataframe. this can then be exported as a .csv file and viewed in excel.
    gene_dict = {}
    for i,r in up.iterrows():
        gene_dict[r[col_names[0]]] = ['lengthened',r[col_names[1]]]
        #plt.text(x=r['InfectedvsUninfected.cFC'],y=-np.log10(r['InfectedvsUninfected.adjusted.P.Value']),s=i)
        
    for i,r in down.iterrows():
        gene_dict[r[col_names[0]]] = ['shortened',r[col_names[1]]]
        #plt.text(x=r['InfectedvsUninfected.cFC'],y=-np.log10(r['InfectedvsUninfected.adjusted.P.Value']),s=i)    
    
    gene_df = pd.DataFrame.from_dict(gene_dict, orient='index', 
        columns=['short_or_long',col_names[1]])
    
    #Only get the file for when there is a cFC threshold. Not interested in the genes
    #with 3' UTR change with no cFC threshold.
    if fc_thres > 0:
        gene_df.to_csv(out_path+filename[:-4]+'.csv')
    
    
    #This plots the volcano plot
    plt.axvline(-fc_thres,color="grey",linestyle="--") #plots line for cFC threshold
    plt.axvline(fc_thres,color="grey",linestyle="--") #plots line for cFC threshold
    plt.grid(color='gray', linestyle='--', linewidth=0.5) #makes grid lines
    plt.xlabel(labels[0]) #x-axis label
    plt.ylabel(labels[1]) #y-axis label
    plt.axhline(-np.log10(0.05),color="grey",linestyle="--") #plots the -log10(adj p-val)
    plt.legend() #plots legend
    plt.title(filename[:-4] + ' vs ' + control_str) #plots title based off file name and control
    plt.rcParams['svg.fonttype'] = 'none' #allows you to save it as .svg
    plt.savefig(out_path+filename[:-4]+'_fc'+str(fc_thres)+'.svg',bbox_inches='tight') #saves figure to out_path
    plt.show() #shows plot in spyder IDE


    return gene_df 

#This function does two things:
    #1. get the volcano graph based on the inputted arguments
    #2. get a list of the downregulated genes and upregulated genes
def DGE_graph(filename, col_names, title, control_str, fc_thres, out_path, in_path, labels):
    
    #reads in the file as a dataframe
    path = in_path+filename
    df = pd.read_csv(path)
    
    #plots all the genes based of Log2FC and -log10(adj p-val)
    plt.scatter(x=df[col_names[1]],y=df[col_names[2]].apply(lambda x:-np.log10(x)),s=1,label="Not significant",color='black')

    # subsets the upregulated and downregulated genes
    down = df[(df[col_names[1]]<=-fc_thres)&(df[col_names[2]]<=0.05)]
    up = df[(df[col_names[1]]>=fc_thres)&(df[col_names[2]]<=0.05)]
    
    #replots the up and downreg genes but with colors
    plt.scatter(x=down[col_names[1]],y=down[col_names[2]].apply(lambda x:-np.log10(x)),s=3,label="Downregulated n="+str(len(down)),color='xkcd:ocean blue')
    plt.scatter(x=up[col_names[1]],y=up[col_names[2]].apply(lambda x:-np.log10(x)),s=3,label="Upregulated n="+str(len(up)),color='xkcd:coral')
    

    plt.axvline(-fc_thres,color="grey",linestyle="--") #pos fold change line
    plt.axvline(fc_thres,color="grey",linestyle="--") #neg fold change line
    plt.grid(color='gray', linestyle='--', linewidth=0.5) #add grid lines
    plt.xlabel(labels[0]) #x-axis label
    plt.ylabel(labels[1]) #y-axis label
    plt.xlim(-4,4) #this just makes the volcano plot smaller but it loses some datapoints
    plt.ylim(0,150) #so, I would recommend not adding this unless really necessary
    

    plt.axhline(-np.log10(0.05),color="grey",linestyle="--") #-log10(0.05) line
    plt.legend() #plots legend
    plt.title(filename[:-4] + ' vs ' + control_str) #plot title
    plt.rcParams['svg.fonttype'] = 'none' #allows it to be saved as an .svg
    plt.savefig(out_path+filename[:-4]+'_fc'+str(fc_thres)+'.svg',bbox_inches='tight') #saves plot in out_path
    plt.show() #shows plot in IDE



# these are the paths to the folders containing the repac or dge output
# and the paths to where you want to save the pictures!
# You will need to change these to your paths
pm_out = '/Users/nstegman/Desktop/hultquist_rotation/final_plots/pm_repac/'
macro_base = '/Users/nstegman/Desktop/hultquist_rotation/primary_macro_repac_out/'#pm_out = '/Users/nstegman/Desktop/hultquist_rotation/final_plots/PM_repac/'
ht1080_base = '/Users/nstegman/Desktop/hultquist_rotation/repac_out_ht1080_new/'
ht_out = '/Users/nstegman/Desktop/hultquist_rotation/final_plots/HT1080_repac/'
dge_ht1080_base = '/Users/nstegman/Desktop/hultquist_rotation/new_dge_out/'
dge_out = '/Users/nstegman/Desktop/hultquist_rotation/final_plots/HT1080_10_dge/'


#since the dge and repac output have different column names for the gene name, fold change, and adj p-val
#I need to specify what the column names are. 
repac_colnames = ['InfectedvsUninfected.gene_name', 'InfectedvsUninfected.cFC', 'InfectedvsUninfected.adjusted.P.Value']
#The very first column in the dge output is empty, so pandas labels it as 'Unnamed: 0'
dge_colnames = ['Unnamed: 0', 'log2FoldChange', 'padj']

#These are what I want my axes labels to be.
dge_labels = ['logFC',"-log10 (p-value)"]
repac_labels = ["Change in 3' UTR Length (cFC)","-log10 (p-value)"]

#makes a list of the file names contained in the folder
ht1080_files = os.listdir(ht1080_base)
#iterates through the files in this list. gets the volcano plot both with and without 
#the cFC threshold. 
for i in ht1080_files:
    get_another_graph(i, repac_colnames,i[:-4],'Uninfected',0, ht_out, ht1080_base,repac_labels)
    get_another_graph(i, repac_colnames,i[:-4],'Uninfected',0.25, ht_out, ht1080_base,repac_labels)


#makes a list of the file names contained in this folder
dge_files = os.listdir(dge_ht1080_base)

#iterates through these folders
for i in dge_files:
    #if this string is in the filename, then run these functions. this makes sure 
    #that you arent running hidden files sometimes that can cause an issue.
    if '_out' in i:
        DGE_graph(i,dge_colnames ,i[:-4],'_Uninfected_24h',0, dge_out, dge_ht1080_base, dge_labels)
        DGE_graph(i,dge_colnames,i[:-4],'_Uninfected_24h',1, dge_out, dge_ht1080_base,dge_labels)

#get the macro volcano plots but with no log2 FC threshold...instead of 0.25 or -0.25 make it 0.   
#Testing volcano plots from concatenated bams

#makes a list of the file names contained in this folder
macro_files = os.listdir(macro_base)

#iterates through these file names
for i in macro_files:
    #I only want the output of the repac from the concatenated bam files, so
    #I made sure I only want to run files with 'concat' in it. You can get rid of this if you 
    #do not need.
    if 'concat' in i:
        #gets the volcano plot with and without fold change thresholds
        get_another_graph(i, repac_colnames,i[:-4],'Uninfected',0, pm_out, macro_base, repac_labels)
        get_another_graph(i, repac_colnames,i[:-4],'Uninfected',0.25, pm_out, macro_base,repac_labels)
