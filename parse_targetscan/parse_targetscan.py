import os
import pandas as pd

'''
parse targetscan
The funtions of this script:
    - once you have ran through your genes of interest through Target scan, and
        saved the table as a csv, make sure that these tables are all in the same folder.
    - rename each file to gene name of interest so you know who is who!
    - Then this script will see what miRNAs are shared among these genes of interest.
    -I am only looking at the columns "conserved sites total" and only in mammals. not vertebrates!!
'''

def get_counts(f_list, in_path, num_shared, out_path, out_fname):
    
    '''
    f_list: list of file names of the TargetScan output tables
    in_path: the path to the folder that holds TargetScan output tables
    '''
    
    #Initialize a dictionary 
    big_dict = {}
    
    #for each file in f_list
    for i in f_list:    
        #My Mac has this hidden file in my folder so I just have to make sure
        #that it doesn't not get inputted or else it will error
        if i != '.DS_Store': 
            #make the full path to the file
            f_path = in_path+i
            #read file in as a dataframe
            df = pd.read_csv(f_path,sep='\t')
            
            # I only want: miRNA family is conserved in: (meh chop this out for now tbh)
            # and: miRNA family
            # and: Conserved sites total
            #
            #just gets the miRNAs conserved in mammals
            mammal_df = df.loc[df['miRNA family is conserved in:'] == 'Mammals']
            
            #Makes a list of the miRNA families in this df
            miRNA_list = mammal_df['miRNA family'].tolist()
            
            #Makes a list of the number of total conserved sites for a given miRNA family
            gene_counts_list = mammal_df['Conserved sites total'].tolist()
            
            #This creates a dictionary of just these two lists of one file
            little_dict = {'miRNA_list':miRNA_list, 'conserved_sites_total':gene_counts_list}
            
            #this adds the little dict from the file you are on to the bid_dict
            big_dict[i] = little_dict
    
    #Once I have the big dict, I get the total list of genes among all the input files
    list_mirnas = []
    for key, value in big_dict.items():
        for j in value['miRNA_list']:
            list_mirnas.append(j)
    
    #I then sort this list and remove duplicates
    all_miRNAs = sorted(set(list_mirnas))
    
    #This initializes the dataframe so that the miRNAs from the list with no duplicates
    #is used as the index, and the columns are the keys of the big_dict, which are the 
    #gene names
    df = pd.DataFrame(0, index=all_miRNAs, columns=big_dict.keys())

    #This adds the gene name, and conserved sites total number in the correct place in the df
    for gene, gene_data in big_dict.items():
        for miRNA, count in zip(gene_data['miRNA_list'], gene_data['conserved_sites_total']):
            df.at[miRNA, gene] = count
            
    # I transpose this df so that the gene names are the rows and the
    # miRNA families are the columns
    final_df = df.transpose()
    
    #This add another row that gets the total count of conserved sites 
    final_df.loc['Total'] = final_df.sum(axis=0)
    
    final_df.loc['Total_unique'] = final_df.astype(bool).sum()
    
    # Now I want to make a df with just columns with total unique > certain threshold
    fil_final_df = final_df.loc[:, final_df.loc['Total_unique'] >= num_shared]
    
    #exports the two dataframes into csv files
    final_df.to_csv(out_path+out_fname)
    fil_final_df.to_csv(out_path+out_fname[:-4]+'_FIL.csv')


# this is for down+down T KO REPAC vs. HT1080 NVP REPAC
KO_dd_path = '/Users/nstegman/Desktop/hultquist_rotation/targetscan/T_KO_dge_v_repac_targetscanout_dd/'
KO_dd_out = '/Users/nstegman/Desktop/hultquist_rotation/targetscan/T_KO_dge_v_repac_Tscan_out_dd/'
KO_dgeR_dd_files = os.listdir(KO_dd_path)

get_counts(KO_dgeR_dd_files, KO_dd_path,10,KO_dd_out,'R_CPSF6_KO_shared_dgeR_dd_miRNAs.csv')
