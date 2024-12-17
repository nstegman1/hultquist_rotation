import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

'''
This script's functions:
- gets a csv of all of the shared genes of various conditions, including:
    - REPAC vs DGE (the get_correlation_DGE() function):
        - all shared genes 
        - shortened 3' UTR and downregulated genes
        - lengthend 3' UTR and upregulated genes
    - REPAC v REPAC (the get_correlation_fil() function):
        - all shared genes
        - shortened 3' UTR of one condition and shortened 3' UTR of a different condition
        
- This script also produces correlation plots but these do not include the gene count labels
    for each quadrant.
    - I still left it in just in case you wanted correlation plots without gene counts
'''

def get_correlation_DGE(fname1, fname2, out_p, plot_title, labels):    
    
    #reads in files as dataframes
    f1_df = pd.read_csv(fname1)
    f2_df = pd.read_csv(fname2)
    
    #gets just the file name from fname1 or 2. Not the full path.
    split1 = fname1.split('/')
    f1 = split1[-1]
    split2 = fname2.split('/')
    f2 = split2[-1]
    
    #This function is only for repac v DGE plots, so this gets column names from the 
    #file. REMEMBER to always input the repac file as fname1! DGE is always fname2!
    cFC_column_name1 = f1_df.filter(regex='cFC$').columns[0]
    gn_column_name1 = f1_df.filter(regex='gene_name$').columns[0]
    p_column_name1 = f1_df.filter(regex='adjusted.P.Value$').columns[0]
    
    #This function is only for repac v DGE plots, so column name for DGE are always these
    cFC_column_name2 = 'log2FoldChange'
    gn_column_name2 = 'Unnamed: 0'
    p_column_name2 = 'padj'
    

    #extracts the threes columns we want: gene name, fold change, and adj p-value
    f1_df_trunc = f1_df[[gn_column_name1,cFC_column_name1,p_column_name1]]
    f2_df_trunc = f2_df[[gn_column_name2,cFC_column_name2,p_column_name2]]
    
    # Rename columns by setting df.columns to the new names
    new_column_names = ['gene_name','cFC','adj_p_value']
    f1_df_trunc.columns = new_column_names
    f2_df_trunc.columns = new_column_names
    
    # Merge on the 'gene' column to get shared genes
    # ONLY includes the gene names that exist in both dfs
    shared_genes = pd.merge(f1_df_trunc, f2_df_trunc, on='gene_name', suffixes=('_'+f1[:-4], '_'+f2[:-4]))
    
    #removes duplicate genes (this happens when there is a gene name, 
    #but there are different fold changes and p vals for the same gene name)
    #we loose some information, but it does not really matter because we are only looking
    #at shared genes and not accounting for p-vals or fold changes.
    background_genes = shared_genes['gene_name'].unique()
    
    #gets the genes that are have shortened 3' UTR and downregulated
    shared_down_down = shared_genes[
    (shared_genes['cFC_' + f1[:-4]] <= -0.25) &
    (shared_genes['cFC_' + f2[:-4]] <= 0) &
    (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
    (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
    
    #removes duplicates of the downregulated and shortened 3' UTR gene list 
    down_down_list = shared_down_down['gene_name'].unique()
    
    #gets the genes that are have shortened 3' UTR and upregulated
    shared_down_up = shared_genes[
    (shared_genes['cFC_' + f1[:-4]] >= 0.25) &
    (shared_genes['cFC_' + f2[:-4]] >= 0) &
    (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
    (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
    
    #removes duplicates of the upregulated and shortened 3' UTR gene list 
    down_up_list = shared_down_up['gene_name'].unique()
    
    # Calculates correlation between fold changes if you want
    correlation = shared_genes['cFC_'+f1[:-4]].corr(shared_genes['cFC_'+f2[:-4]])
    print('Correlation: '+str(correlation))
    
    #initiates plot
    fig1, ax1 = plt.subplots()
    
    #makes the scatter plot using the shared genes dataframe (no p val threshold)
    sns.scatterplot(x='cFC_'+f1[:-4], y='cFC_'+f2[:-4], data=shared_genes)
    
    #plots a grid
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
   
    #adds axes labels and plot title
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.title(plot_title)
    
    #saves plot with no p val threshold
    plt.savefig(out_p+f1[:-4]+'_vs_'+f2[:-4]+'_unfil.png',format='png', dpi=300, bbox_inches='tight')
    plt.show()
    
    #filters the shared gene list with p = 0.05
    shared_genes_fil = shared_genes[
    (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
    (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
    
    #creates the scatter plot of repac vs DGE fold change
    sns.scatterplot(x='cFC_'+f1[:-4], y='cFC_'+f2[:-4], data=shared_genes_fil)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    
    #plots title and axes labels
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.title(plot_title)
    
    #saves plot 
    plt.savefig(out_p+f1[:-4]+'_vs_'+f2[:-4]+'_fil.png',format='png', dpi=300, bbox_inches='tight')
    plt.show()
      
    #this returns the background genes, the shortened 3' UTR and downregulated genes, and the shortened 3' UTR and upregulated genes
    #in a dictionary
    final_dict = {f1[:-4]+ ' and '+f2[:-4]+'_RvDGE_background_genes':background_genes, f1[:-4]+ ' and '+f2[:-4]+'_RvDGE_down_down_genes':down_down_list,f1[:-4]+' and '+f2[:-4]+'_RvDGE_down_up_genes':down_up_list}
    return final_dict
    
    

def get_correlation_fil(fname1, fname2, need_pval, out_path):
    '''
    This gets the correlation plots for just RvR comparisons
    
    fname1: path to the file of first REPAC output
    fname2: path to the file of second REPAC output
    need_pval: if True, then filter the dataset you want to plot with p-val.
        If False, then do not filter dataset based on p-value.
    out_path: The folder where you want the plots to be exported to
    '''
    
    #reads in the two files as dataframes
    f1_df = pd.read_csv(fname1)
    f2_df = pd.read_csv(fname2)
    
    #gets just the name of the file 
    split1 = fname1.split('/')
    f1 = split1[-1]
    
    split2 = fname2.split('/')
    f2 = split2[-1]
    
    cFC_column_name1 = f1_df.filter(regex='cFC$').columns[0]
    gn_column_name1 = f1_df.filter(regex='gene_name$').columns[0]
    p_column_name1 = f1_df.filter(regex='adjusted.P.Value$').columns[0]
    
    cFC_column_name2 = f2_df.filter(regex='cFC$').columns[0]
    gn_column_name2 = f2_df.filter(regex='gene_name$').columns[0]
    p_column_name2 = f2_df.filter(regex='adjusted.P.Value$').columns[0]
            
    f1_df_trunc = f1_df[[gn_column_name1,cFC_column_name1,p_column_name1]]
    f2_df_trunc = f2_df[[gn_column_name2,cFC_column_name2,p_column_name2]]
    
    new_column_names = ['gene_name','cFC','adj_p_value']
    
    # Rename columns by setting df.columns to the new names
    f1_df_trunc.columns = new_column_names
    f2_df_trunc.columns = new_column_names
    
    #get the shared genes
    shared_genes = pd.merge(f1_df_trunc, f2_df_trunc, on='gene_name', suffixes=('_'+f1[:-4], '_'+f2[:-4]))
    
    #removes duplicates from the gene names of the shared_genes df.
    #this gets the total genes shared by both conditions
    back_genes = shared_genes['gene_name'].unique()
 
    #If we want to filter by a p-value, then do this
    if need_pval == True:
        
        #This created a list of genes that have shortened 3' UTRs
        #and are down regulated filtered by p-value
        shared_down_genes = shared_genes[
        (shared_genes['cFC_' + f1[:-4]] <= -0.25) &
        (shared_genes['cFC_' + f2[:-4]] <= -0.25) &
        (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
        (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
        
        #this gets the unique gene names from this dataframe
        down_down_genes = shared_down_genes['gene_name'].unique()
        
        #shared_down_genes.to_csv(shared_path+f1[:-4]+'_vs_'+f2[:-4]+'_fil_llqV2.csv')
    
    else: 
        #This created a list of genes that have shortened 3' UTRs
        #and are down regulated not filtered by p-value
        shared_down_genes = shared_genes[
        (shared_genes['cFC_' + f1[:-4]] <= -0.25) &
        (shared_genes['cFC_' + f2[:-4]] <= -0.25)]
        
        #this gets the unique gene names from this dataframe
        down_down_genes = shared_down_genes['gene_name'].unique()
        
        #shared_down_genes.to_csv(shared_path+f1[:-4]+'_vs_'+f2[:-4]+'_unfil_llqV2.csv')
    
    f1_cfc = shared_genes.iloc[:,2]  # Selects the first and third columns
    f2_cfc = shared_genes.iloc[:,3]
    
    
    corr, _ = pearsonr(f1_cfc, f2_cfc)
    print('Pearson correlation of '+f1+' vs '+f2+': '+str(corr))

    #This makes correlation plots
    if need_pval == True:
        correlation_in = shared_genes[
        (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
        (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
        # Plot
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='cFC_'+f1[:-4], y='cFC_'+f2[:-4], data=correlation_in)
        plt.grid(color='gray', linestyle='--', linewidth=0.5)
        plt.axvline(-0.25,color="grey",linestyle="--")
        plt.axvline(0.25,color="grey",linestyle="--")
        plt.axhline(0,color="grey",linestyle="--")
        plt.xlabel(f1[:-4]+' Expression')
        plt.ylabel(f2[:-4]+' Expression')
        plt.title('Gene Expression Correlation between '+f1[:-4]+ ' and '+f2[:-4])
        #plt.show() #shows the plot in your IDE
        plt.savefig(out_path+f1[:-4]+'_vs_'+f2[:-4]+'_rvr_fil.png',format='png', dpi=300, bbox_inches='tight')
    else:
        correlation_in = shared_genes
        # Plot
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='cFC_'+f1[:-4], y='cFC_'+f2[:-4], data=correlation_in)
        plt.grid(color='gray', linestyle='--', linewidth=0.5)  
        plt.axvline(-0.25,color="grey",linestyle="--")
        plt.axvline(0.25,color="grey",linestyle="--")
        plt.axhline(0,color="grey",linestyle="--")
        plt.xlabel(f1[:-4]+' Expression')
        plt.ylabel(f2[:-4]+' Expression')
        plt.title('Gene Expression Correlation between '+f1[:-4]+ ' and '+f2[:-4])
        #plt.show() #shows the plot in your IDE
        plt.savefig(out_path+f1[:-4]+'_vs_'+f2[:-4]+'_rvr_UNfil.png',format='png', dpi=300, bbox_inches='tight')
        

    final_dict = {f1[:-4]+ ' and '+f2[:-4]+'RvR_background_genes':back_genes, f1[:-4]+ ' and '+f2[:-4]+'RvR_down_down_genes':down_down_genes }

    
    return final_dict


     
#REPAC OUT
base = '/Users/nstegman/Desktop/hultquist_rotation/gene_analysis/repac_out/'
#out_path = '/Users/nstegman/Desktop/hultquist_rotation/FINAL_DATA/cor_plots_test/'
out_path = '/Users/nstegman/Desktop/hultquist_rotation/FINAL_DATA/shared_gene_data/corr/'
base_KO = '/Users/nstegman/Desktop/hultquist_rotation/targetscan/CPSF6_KO_repac_DGE_data/'
ht1080_nvp_wt_fname = base+'HT1080_WT_NVP_10.csv'
ht1080_wt_fname = base+'HT1080_WT_DMSO_10.csv'
ht1080_n74d_fname = base+'HT1080_N74D_DMSO_10.csv'
T_cell_wt_fname = base+'Primary_CD4_WT_vs_Uninfected_REPAC.csv'
T_cell_KO_fname = base_KO+'CPSF6_KO_Primary_Cell_REPAC.csv'
T_cell_N74D_fname = base+'Primary_CD4_N74D.csv'


#DGE OUT
base_dge = '/Users/nstegman/Desktop/hultquist_rotation/gene_analysis/dge_out/'
ht1080_WT_NVP_DGE = base_dge+'ht1080_WT_10_NVP_diffreg.csv'
ht1080_WT_DGE = base_dge+'WT_10_DMSO_diffreg.csv'
ht1080_N74D_DGE = base_dge+'N74D_10_DMSO_diffreg.csv'
T_KO_dge = base_KO+'all3donors_CPSF6vNT_includesallsamples.csv'

#MAKE SURE TO SET YOUR FUNCTIONS EQUAL TO A VARIABLE! we will be 
#doing something with this later 

# index 0: background, 1: down down of RvR
ht1080_nvp_v_t_wt_f = get_correlation_fil(ht1080_nvp_wt_fname,T_cell_wt_fname, True, out_path)
ht1080_nvp_v_t_KO_f = get_correlation_fil(ht1080_nvp_wt_fname,T_cell_KO_fname, True, out_path) # misnamed?
ht1080_nvp_v_t_n74d_f = get_correlation_fil(ht1080_nvp_wt_fname,T_cell_N74D_fname, True, out_path)


#ALWAYS PUT REPAC CSV as first argument!!! 
# index 0: background, 1: down down, 2: down up
al1 = ['HT1080 WT NVP REPAC (cFC)', 'HT1080 WT NVP DGE (Log2FC)']
ft1 = 'Gene Expression Correlation HT1080 WT NVP REPAC vs. HT1080 WT NVP DGE'
ht1080_DGE_v_repac_NVP_WT = get_correlation_DGE(ht1080_nvp_wt_fname, ht1080_WT_NVP_DGE, out_path, ft1,al1)

al2 = ['HT1080 WT REPAC (cFC)', 'HT1080 WT DGE (Log2FC)']
ft2 = 'Gene Expression Correlation HT1080 WT REPAC vs. HT1080 WT DGE'
ht1080_DGE_v_repac_WT = get_correlation_DGE(ht1080_wt_fname, ht1080_WT_DGE, out_path,ft2,al2)

al3 = ['HT1080 N74D REPAC (cFC)', 'HT1080 N74D DGE (Log2FC)']
ft3 = 'Gene Expression Correlation HT1080 N74D REPAC vs. HT1080 N74D DGE'
ht1080_DGE_v_repac_N74D = get_correlation_DGE(ht1080_n74d_fname,ht1080_N74D_DGE, out_path,ft3,al3)


final_title1 = 'Gene Expression Correlation CPSF6 KO T Cell REPAC vs. CPSF6 KO T Cell DGE'
axes_labels1 = ['CPSF6 KO T Cell REPAC (cFC)', 'CPSF6 KO T Cell DGE (Log2FC)']
# this is for repac vs dge CPSF6 KO shared gene list
T_KO_DGE_v_repac = get_correlation_DGE(T_cell_KO_fname,T_KO_dge,out_path, final_title1, axes_labels1)


'''
The variables that you set your functions equal to hold the shared gene lists 
in a data structure called a dictionary. In order to organize all this output into 
a csv file, we have to "unwind" all these dictionaries and combine them into 
a single dataframe
'''

#makes a list of the dictionaries from the functions you called above
lol = [ht1080_nvp_v_t_wt_f,ht1080_nvp_v_t_KO_f,ht1080_nvp_v_t_n74d_f,ht1080_DGE_v_repac_NVP_WT,ht1080_DGE_v_repac_WT,ht1080_DGE_v_repac_N74D,T_KO_DGE_v_repac]
#initializes an empty dictionary
unwound = dict()

#this for loop adds the items from each dictionary into the unwound dictionary
for i in lol:
    for key,value in i.items():
        print(key,value)
        unwound[key] = value

#this converts the unwound dictionary into a dataframe.
#but because each value (gene list) in unwound is a different length,
#you have to fill the gene lists with NA's so that they are all the same length
df = pd.DataFrame(dict([(col_name,pd.Series(values)) for col_name,values in unwound.items()]))

#export the final df to where you want on your computer
df.to_csv('/Users/nstegman/Desktop/hultquist_rotation/FINAL_DATA/shared_genes.csv')

