import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
This script's function:
    - Make correlation plots between repact vs repac data or
    repac vs dge data.
    
    - These plots also count the genes in each quadrant and labels that
    on the plot
'''

#This function gets the counts of each quadrant. For now, the x and y intercept are 0,
#but you can change that if needed
def get_gene_counts(df, x_intercept, y_intercept, f1, f2):
    # Define quadrant boundaries
    x_intercept = 0
    y_intercept = 0
    
    # Count points in each quadrant
    q1 = np.sum((df['cFC_'+f1[:-4]] > x_intercept) & (df['cFC_'+f2[:-4]] > y_intercept)) #top right
    q2 = np.sum((df['cFC_'+f1[:-4]] < x_intercept) & (df['cFC_'+f2[:-4]] > y_intercept)) #top left
    q3 = np.sum((df['cFC_'+f1[:-4]] < x_intercept) & (df['cFC_'+f2[:-4]] < y_intercept)) #bottom left
    q4 = np.sum((df['cFC_'+f1[:-4]] > x_intercept) & (df['cFC_'+f2[:-4]] < y_intercept)) #bottom right
    
    print(q1,q2,q3,q4)
    
    return q1,q2,q3,q4

#This plots the shared genes between two conditions
def get_plot (data,x_threshold, y_threshold, f1, f2, labels, plot_title, out_p, distance, distance2, is_fil) :  
    fig, ax = plt.subplots()
    ax.scatter(x=data['cFC_'+f1[:-4]], y=data['cFC_'+f2[:-4]],s=15,edgecolor='white', linewidth=.5)
    
    points = get_gene_counts(data, x_threshold, y_threshold, f1, f2)
    
    # Mark the section of interest
    ax.axvline(x_threshold, color='gray', linestyle='--')
    ax.axhline(y_threshold, color='gray', linestyle='--')
    
    # Find an empty space for the label
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()


    y_min -= 0.5
    y_max -= 0.5
    
    # Default placement for annotation
    label_x_ll = (x_min-1-distance-distance2 + x_threshold) / 2
    label_y_ll = (y_min-1-distance-distance2 + y_threshold) / 2
    
    # Annotate the count on the plot
    ax.text(
        label_x_ll, label_y_ll, 'n= '+str(points[2]), 
        fontsize=12, color="black", 
        ha='center', va='center',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
    
    
    label_x_lr = (x_max+.7+distance-distance2 + x_threshold) / 2
    label_y_lr = (y_min-1-distance-distance2 + y_threshold) / 2
    ax.text(
        label_x_lr, label_y_lr, 'n= '+str(points[3]), 
        fontsize=12, color="black", 
        ha='center', va='center',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))


    label_x_upl = (x_min-1-distance + x_threshold) / 2
    label_y_upl = (y_max+1+distance + y_threshold) / 2
    ax.text(
        label_x_upl, label_y_upl, 'n= '+str(points[1]), 
        fontsize=12, color="black", 
        ha='center', va='center',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
    
    label_x_upr = (x_max+.7+distance + x_threshold) / 2
    label_y_upr = (y_max+1+distance + y_threshold) / 2

    ax.text(
        label_x_upr, label_y_upr, 'n= '+str(points[0]), 
        fontsize=12, color="black", 
        ha='center', va='center',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.title(plot_title)
    
    if is_fil == True:
        plt.savefig(out_p+f1[:-4]+'_vs_'+f2[:-4]+'_fil.png',format='png', dpi=300, bbox_inches='tight')
    else:
        plt.savefig(out_p+f1[:-4]+'_vs_'+f2[:-4]+'_unfil.png',format='png', dpi=300, bbox_inches='tight')
    plt.show()
    
#this makes the shared gene lists for the correlation plots and further analyses (GSEA, TargetScan)
def get_correlation(fname1, fname2, out_p, plot_title, labels, distance,distance2,is_both):
    
    #change this path to where you want your shared gene lists (the csv files) to export to 
    shared_path = '/Users/nstegman/Desktop/hultquist_rotation/gene_analysis/shared_out_repac_v_dge/'
    
    #reads in the genes we are comparing against each other as dataframes
    f1_df = pd.read_csv(fname1)
    f2_df = pd.read_csv(fname2)
    
    #this makes sure we just get the file name of each file alone and not the whole path string
    split1 = fname1.split('/')
    f1 = split1[-1]
    split2 = fname2.split('/')
    f2 = split2[-1]
    
    #If you are doing a REPAC vs. DGE comparison, set the is_both argument to False, because 
    #not both of the files are REPAC output. This is done because the DGE output has 
    # different column names than the REPAC output.
    if is_both == False:
        
        #This extracts the column names of the input DGE file
        cFC_column_name2 = 'log2FoldChange'
        gn_column_name2 = 'Unnamed: 0'
        p_column_name2 = 'padj'
        
        #This extracts the column names of the input REPAC file that we want to extract.
        cFC_column_name1 = f1_df.filter(regex='cFC$').columns[0] #col name of cFC column
        gn_column_name1 = f1_df.filter(regex='gene_name$').columns[0] #col name of gene name column
        p_column_name1 = f1_df.filter(regex='adjusted.P.Value$').columns[0] #col name of adj p-pal column
        
   
    else:
        #if is_both is set to True, then that means both input files are REPAC output,
        #which means the column names we want are in the same format.
        
        #This extracts the column names of the first REPAC file that we want to extract.
        cFC_column_name1 = f1_df.filter(regex='cFC$').columns[0]
        gn_column_name1 = f1_df.filter(regex='gene_name$').columns[0]
        p_column_name1 = f1_df.filter(regex='adjusted.P.Value$').columns[0]
        
        #This extracts the column names of the second input REPAC file that we want to extract.
        cFC_column_name2 = f2_df.filter(regex='cFC$').columns[0]
        gn_column_name2 = f2_df.filter(regex='gene_name$').columns[0]
        p_column_name2 = f2_df.filter(regex='adjusted.P.Value$').columns[0]
    
    #this makes new dataframes that only have the three columns we need from each input file
    f1_df_trunc = f1_df[[gn_column_name1,cFC_column_name1,p_column_name1]]
    f2_df_trunc = f2_df[[gn_column_name2,cFC_column_name2,p_column_name2]]
     
    #This renames the df's column names so that we can  merge them in the next step
    new_column_names = ['gene_name','cFC','adj_p_value']
    # Rename columns by setting df.columns to the new names
    f1_df_trunc.columns = new_column_names
    f2_df_trunc.columns = new_column_names
    
    #pd.merge merges the two dataframes based on the 'gene_name' column
    #this means that in order for this merged data frame to include a row, both dataframes
    #MUST have that gene name. Therefore, this df ONLY includes the gene names that
    #exist in both dataframes
    shared_genes = pd.merge(f1_df_trunc, f2_df_trunc, on='gene_name', suffixes=('_'+f1[:-4], '_'+f2[:-4]))

    #time to plot fil and unfil
    #def get_DGE_plot (data,x_threshold, y_threshold, f1, f2, labels, plot_title, out_p, distance) :  
    shared_genes_fil = shared_genes[
    (shared_genes['adj_p_value_' + f1[:-4]] < 0.05) &  # Add p-value filter for f1
    (shared_genes['adj_p_value_' + f2[:-4]] < 0.05)]    # Add p-value filter for f2
    
    
    get_plot(shared_genes_fil,0,0,f1,f2,labels, plot_title, out_p, distance,distance2,True)
    get_plot(shared_genes,0,0,f1,f2,labels, plot_title+' no p-value', out_p, distance, distance2, False)
    
    # Calculate correlation (if you want)
    correlation = shared_genes['cFC_'+f1[:-4]].corr(shared_genes['cFC_'+f2[:-4]])
    print(f'Correlation: {correlation}')


#REPAC OUT
base = '/Users/nstegman/Desktop/hultquist_rotation/gene_analysis/repac_out/'
out_path = '/Users/nstegman/Desktop/hultquist_rotation/ppt_data/correlation_plots/'
base_KO = '/Users/nstegman/Desktop/hultquist_rotation/targetscan/CPSF6_KO_repac_DGE_data/'
ht1080_nvp_wt_fname = base+'HT1080_WT_NVP_10.csv'
ht1080_wt_fname = base+'HT1080_WT_DMSO_10.csv'
ht1080_n74d_fname = base+'HT1080_N74D_DMSO_10.csv'
T_cell_wt_fname = base+'Primary_CD4_WT_vs_Uninfected_REPAC.csv'
T_cell_KO_fname = base_KO+'CPSF6_KO_Primary_Cell_REPAC.csv'
T_cell_N74D_fname = base+'Primary_CD4_N74D.csv'


#DGE OUT
base_dge = '/Users/nstegman/Desktop/hultquist_rotation/new_dge_out/'
ht1080_WT_NVP_DGE = base_dge+'count_matrix_WT_10_NVP_out.csv'
ht1080_WT_DGE = base_dge+'count_matrix_WT_10_DMSO_out.csv'
ht1080_N74D_DGE = base_dge+'count_matrix_N74D_10_DMSO_out.csv'
T_KO_dge = base_KO+'all3donors_CPSF6vNT_includesallsamples.csv'


final_title1 = 'Gene Expression Correlation CPSF6 KO T Cell REPAC vs. CPSF6 KO T Cell DGE'
axes_labels1 = ['CPSF6 KO T Cell REPAC (cFC)', 'CPSF6 KO T Cell DGE (Log2FC)']
# this is for repac vs dge CPSF6 KO shared gene list
T_KO_DGE_v_repac = get_correlation(T_cell_KO_fname,T_KO_dge,out_path, final_title1, axes_labels1,0,0, False)

final_title2 = 'Gene Expression Correlation WT NVP HT1080 Cell REPAC vs. WT NVP HT1080 Cell DGE'
axes_labels2 = ['WT NVP HT1080 REPAC (cFC)', 'WT NVP HT1080 DGE (Log2FC)']
distance_fil = 0.1
distance_unfil = 1
# this is for repac vs dge HT1080 nvp shared gene list
ht1080nvp_DGE_v_repac = get_correlation(ht1080_nvp_wt_fname,ht1080_WT_NVP_DGE,out_path, final_title2, axes_labels2,-0.5,0, False)

# this is for HT1080 NVP REPAC vs. CPSF6 KO T cell REPAC
final_title3 = 'Gene Expression Correlation WT NVP HT1080 Cell REPAC vs. CPSF6 KO T Cell REPAC'
axes_labels3 = ['WT NVP HT1080 REPAC (cFC)', 'CPSF6 KO T Cell REPAC (cFC)']
distance_fil = 0.1
distance_unfil = 1
ht1080_nvp_v_t_KO_f = get_correlation(ht1080_nvp_wt_fname,T_cell_KO_fname, out_path,final_title3, axes_labels3,-.5,0, True)
