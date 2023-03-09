#!/usr/bin/env python3
#
# PlotSpecScript
#
# Analyze data with PlotSpec functions.
# Mgmt mice data.
#
#  B. Fedeles, 2022.

import PlotSpec as ps
import ClustPlot as cp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

print(ps.pu_muts)
print(ps.py_muts)

allmuts = ps.pu_muts+ps.py_muts

print(allmuts)


def generate_spectra(group, name, title, kmerfile, colorscheme='COSMIC2'):

    # Function to plot average mutational spectra.
    # Group is a list of msp files
    # Name is Figure/Save filename
    # Title is a string to be added to the plot
    # kmerfile is a file with trinucleotide proportions (for normalization)


    #### Define here paths to save files:
    datafilepath="Datafiles3/"
    picfilepath="Pics3/"
    ####

    kmer = ps.import_kmer_counts(kmerfile)

    
    # Printing with purine labels
    pu_dict=ps.init_spec_dict('purine')
    xlab = list(zip(*pu_dict.keys()))[1]
    
    val_list=[]
    std_list=[]
    xlab_list=[]

    print(group)

    new = ps.combine_csv_files(group, op='avg')
    ps.save_csv_file(datafilepath+name[:-5]+'.csv', new)

    #new2=new
    new2 = ps.normalize_spec(new, kmer)
    ps.save_csv_file(datafilepath+name[:-5]+'-norm.csv', new2)

    vals =list(zip(*new2.values()))[0]
    stds =list(zip(*new2.values()))[1]
    #xlab = [list(zip(*new2.keys()))[1]]

    val_list.append(vals)
    std_list.append(stds)
    xlab_list.append(xlab)

    #plotting individual plots
    mpl.rc("savefig", dpi=320)
    ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts, errorbars=[stds],
                   titles=[title], ylabel='Proportion of mutations', colorscheme=colorscheme)
    plt.savefig(picfilepath+name+'.svg')
    plt.savefig(picfilepath+name+'.png')

    return

def main():

    filespath = 'Datafiles/'
    filespath2= 'Datafiles2/'
    picspath = 'Pics/'
    picspath2= 'Pics2/'
    filesuff1 = 'msp.csv'

    #kmerfile = 'bbmap.count.EG10c.txt'
    kmerfile = 'twnstr-mouse-contexts.txt'

    plotspectra=True
    hist= False

    cossim = False
    cossimhist = False

    

    males = ['8212', '8213', '8214']  # males
    females = ['8215', '8216', '8217'] # females
    malectrl = ['8204', '8205', '8206'] # males ctrl
    femalectrl =['8207', '8208', '8209'] # female ctrl

    maletitle= 'Mgmt-/- NDMA-treated mice males n=3 unique normalized'
    femaletitle= 'Mgmt-/- NDMA-treated mice females n=3 unique normalized'
    maletitlectrl = 'Ctrl mice males n=3 unique normalized'
    femaletitlectrl = 'Ctrl mice females n=3 unique normalized'
    allmgmtndma = 'Mmgt-/- NDMA-treated mice n=6 unique normalized'
    allmgmtctrl='Mgmt-/- ctrl mice n=6 unique normalized'

    malegroup = [filespath+f+filesuff1 for f in males]
    femalegroup = [filespath+f+filesuff1 for f in females]
    ctrlmale = [filespath+f+filesuff1 for f in malectrl]
    ctrlfemale = [filespath+f+filesuff1 for f in femalectrl]
    allndma = malegroup+femalegroup
    allctrl = ctrlmale+ctrlfemale

    #generate_spectra(malegroup, 'mgmt-males-avg-norm-no-title', "", kmerfile)
    #generate_spectra(femalegroup, 'mgmt-females-avg-norm-no-title', "", kmerfile)

    
    if plotspectra:
        
        generate_spectra(malegroup, 'mgmt-males-avg-norm', maletitle, kmerfile, 'COSMIC3')
        generate_spectra(femalegroup, 'mgmt-females-avg-norm', femaletitle, kmerfile, 'COSMIC3')
        generate_spectra(ctrlmale, 'mgmt-ctrl-males-avg-norm', maletitlectrl, kmerfile, 'COSMIC3')
        generate_spectra(ctrlfemale, 'mgmt-ctrl-females-avg-norm', femaletitlectrl, kmerfile, 'COSMIC3')

        generate_spectra(allndma, 'all-mgmt-ndma-avg-norm', allmgmtndma, kmerfile, 'COSMIC3')
        generate_spectra(allctrl, 'all-mgmt-ctrl-avg-norm', allmgmtctrl, kmerfile, 'COSMIC3')

    if hist:
    # Clustplot of heatmap'

        print('Making hist')

        datasamples = ['NDMA-WT', 'NDMA-Mgmt','Temozolomide', 'Streptozotocin', 'MNU', 'SBS11']

        datafiles = ['all_NDMA_sum_bgsub', 'all-mgmt-ndma-avg', 'Temozolomide_bgsub', 'Streptozotocin_bgsub',
                     'MNU_bgsub', 'SBS11_32']
        datafiles2= [file+'_norm.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):
            spec_list[sam] = ps.read_csv_file(filespath2+file)


        ext1='png'
        ext2='svg'
        
        for st_col in [2.3]:
##        for st_col in [2.2, 2.4, 2.6, 2.8]:
            plotfile1 = picspath2+'histplot'+str(st_col)+'.'+ext1
            plotfile2 = picspath2+'histplot'+str(st_col)+'.'+ext2

            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, st_col)  
            cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, st_col)  


    if cossim:

        # Cossim between various spectra

        wt_ctrl = ps.read_csv_file(filespath2+'wt-ctrl-all-norm.csv')
        wt_ctrl_vals = list(wt_ctrl.values())

        mgmt_ctrl = ps.read_csv_file(filespath2+'all-mgmt-ctrl-avg-norm.csv')
        mgmt_ctrl_vals = list(mgmt_ctrl.values())
        
        #for sig, file in zip(csigs,cfiles2):

##        sigspec = ps.read_msp_file(file)
##        sigspecvals = list(sigspec.values())

        cossim = cp.cosine_similarity([wt_ctrl_vals], [mgmt_ctrl_vals])
        print('{}  \t  {}'.format('wt vs mgmt ctrl', cossim[0][0])) 


    if cossimhist:

        # cossim between a collection of spectra

        datasamples = ['WT-ctrl', 'WT-males-ctrl', 'WT-fem-ctrl', 'Mgmt-ctrl', 'Mgmt-males-ctrl', 'Mgmt-fem-ctrl']

        datafiles = ['wt-ctrl-all', 'wt-ctrl-males', 'wt-ctrl-females', 'all-mgmt-ctrl-avg', 'mgmt-ctrl-males-avg',
                     'mgmt-ctrl-females-avg']

        datafiles2 = [file+'-norm.csv' for file in datafiles]

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]

        val_list=[]
        xlab_list=[xlab]
        spec_list={}

        for sam, file in zip(datasamples, datafiles2):
            spec_list[sam] = ps.read_csv_file(filespath2+file)


        ext1='png'
        ext2='svg'

        plotfile1 = picspath2+'ctrl-histplot'+str(2.3)+'.'+ext1
        plotfile2 = picspath2+'ctrl-histplot'+str(2.3)+'.'+ext2

        cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile1, ext1, 2.3)  
        cp.plot_uhc_heatmap(spec_list, datasamples, True, plotfile2, ext2, 2.3)  



if __name__=="__main__":
    main()

