#!/usr/bin/env python3
#
# PlotSpecScript4
#
# Analyze data with PlotSpec functions.
#
# Modified for the lung data
# B. Fedeles, Nov 2022

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
    # Group is a list of msp files (with no headers)
    # Name is Figure/Save filename
    # Title is a string to be added to the plot
    # kmerfile is a file with trinucleotide proportions (for normalization)


    #### Define here paths to save files:
    datafilepath="Datafiles/"
    picfilepath="Pics/"
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

    plotspectra=False

    bsub = False
    hist= False
    hist2 = True

    cossim = False
    cossimhist = False

    

    males = ['9675', '9676']  # males NDMA-treated
    females = ['9677', '9678'] # females NDMA-treated
    malectrl = ['9679', '9680'] # males ctrl
    femalectrl =['9681', '9682'] # female ctrl

    maletitle= 'Lungs of NDMA-treated mice males n=2 unique normalized'
    femaletitle= 'Lungs of NDMA-treated mice females n=2 unique normalized'
    maletitlectrl = 'Lungs Ctrl mice males n=2 unique normalized'
    femaletitlectrl = 'Lungs Ctrl mice females n=2 unique normalized'
    allmgmtndma = 'Lungs of NDMA-treated mice n=4 unique normalized'
    allmgmtctrl='Lungs of ctrl mice n=4 unique normalized'

    malegroup = [filespath+f+filesuff1 for f in males]
    femalegroup = [filespath+f+filesuff1 for f in females]
    ctrlmale = [filespath+f+filesuff1 for f in malectrl]
    ctrlfemale = [filespath+f+filesuff1 for f in femalectrl]
    allndma = malegroup+femalegroup
    allctrl = ctrlmale+ctrlfemale

    #generate_spectra(malegroup, 'mgmt-males-avg-norm-no-title', "", kmerfile)
    #generate_spectra(femalegroup, 'mgmt-females-avg-norm-no-title', "", kmerfile)

    
    if plotspectra:
        
        generate_spectra(malegroup, 'lung-males-avg-norm', maletitle, kmerfile, 'COSMIC3')
        generate_spectra(femalegroup, 'lung-females-avg-norm', femaletitle, kmerfile, 'COSMIC3')
        generate_spectra(ctrlmale, 'lung-ctrl-males-avg-norm', maletitlectrl, kmerfile, 'COSMIC3')
        generate_spectra(ctrlfemale, 'lung-ctrl-females-avg-norm', femaletitlectrl, kmerfile, 'COSMIC3')

        generate_spectra(allndma, 'all-lung-ndma-avg-norm', allmgmtndma, kmerfile, 'COSMIC3')
        generate_spectra(allctrl, 'all-lung-ctrl-avg-norm', allmgmtctrl, kmerfile, 'COSMIC3')

    
    if bsub:

        # average together the 4 ctrl and 4 NDMA-treated lung spectra
        # subtract the ctrl background from the NDMA-treated
        
        # allctrl has all the ctrl files
        # allndma has all the ndma files

        kmer = ps.import_kmer_counts(kmerfile)

        sumfiles = ['all_ctr_sum', 'all_NDMA_sum']
        sumspecs = {}

        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]
        
        val_list=[]
        std_list=[]
        xlab_list=[]


        
        datafiles = [allctrl, allndma]
        
        for file, group in zip(sumfiles, datafiles):
            sumspecs[file] = ps.combine_csv_files(group, op='sum')
            ps.save_csv_file(filespath+file+'.csv', sumspecs[file])

        subspec = ps.subtract_background(sumspecs[sumfiles[1]], sumspecs[sumfiles[0]])
        ps.save_csv_file(filespath+'all_NDMA_sum_bgsub.csv', subspec)

        savename = 'all_NDMA_sum_bgsub_norm'

        subspec2 = ps.normalize_spec(subspec, kmer)
        ps.save_csv_file(filespath+savename+'.csv', subspec2)

        #vals =list(zip(*subspec2.values()))[0]
        vals = list(subspec2.values())
        val_list.append(vals)
        xlab_list.append(xlab)
    
        mpl.rc("savefig", dpi=320)
        ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts,
                           titles=['All NDMA bkgr sub n=4'], ylabel='Proportion of mutations')
        plt.savefig(picspath+savename+'.svg')
        plt.savefig(picspath+savename+'.png')

        




    if hist:
    # Clustplot of heatmap'

        print('Making hist')

        datasamples = ['WT-NDMA-liver', 'WT-NDMA-lung', 'Mgmt-NDMA-liver','Temozolomide', 'Streptozotocin', 'MNU', 'SBS11']

        datafiles = ['WT_liver_NDMA_sum_bgsub', 'WT_lung_NDMA_sum_bgsub', 'Mgmt_liver_NDMA_avg', 'Temozolomide_bgsub', 'Streptozotocin_bgsub',
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



    if hist2:

    # Clustplot of heatmap including controls

        print('Making hist')

        datasamples = ['WT-NDMA-liver', 'WT-NDMA-lung', 'Mgmt-NDMA-liver','Temozolomide', 'Streptozotocin', 'MNU', 'SBS11', 'WT-ctrl-liver', 'WT-ctrl-lung', 'Mmgt-ctrl-liver', 'WT-ctrl-MEFs']

        datafiles = ['WT_liver_NDMA_sum_bgsub', 'WT_lung_NDMA_sum_bgsub', 'Mgmt_liver_NDMA_avg', 'Temozolomide_bgsub', 'Streptozotocin_bgsub',
                     'MNU_bgsub', 'SBS11_32', 'WT_liver_ctrl_avg', 'WT_lung_ctrl_avg', 'Mgmt_liver_ctrl_avg', 'WT_MEFs_ctrl_avg']
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
            plotfile1 = picspath2+'histplot-all2'+str(st_col)+'.'+ext1
            plotfile2 = picspath2+'histplot-all2'+str(st_col)+'.'+ext2

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


















##    filenames = ['males_ctr', 'females_ctr', 'males_NDMA', 'females_NDMA', 'all_ctr', 'all_NDMA']
##
##    males_ctr = [354, 355, 358, 77, 78]
##    females_ctr = [356, 357, 74, 75, 76]
##    males_NDMA = [349, 350, 351, 80, 83]
##    females_NDMA = [352, 353, 79, 81, 82]
##    all_ctr =males_ctr+females_ctr
##    all_NDMA=males_NDMA+females_NDMA
##
##
##    datafiles = [males_ctr, females_ctr, males_NDMA, females_NDMA, all_ctr, all_NDMA]
##    titlenames= ['Untreated males, n=5', 'Untreated females, n=5',
##             'NDMA-treated males, n=5', 'NDMA-treated females, n=5',
##             'Untreated mice, n=10', 'NDMA-treated mice, n=10']
##
##    if process_avg:
##
##        #print(datafiles)
##
##        for group in datafiles:
##            for idx, item in enumerate(group):
##                if item>100:
##                    group[idx]=filespath+str(item)+filesuff1
##                else:
##                    group[idx]=filespath+str(item)+filesuff2
##
##        print(datafiles)
##
##        # Printing with purine labels
##        pu_dict=ps.init_spec_dict('purine')
##        xlab = list(zip(*pu_dict.keys()))[1]
##        
##        val_list=[]
##        std_list=[]
##        xlab_list=[]
##
##        for name,titlename, group in zip(filenames, titlenames, datafiles):
##            new = ps.combine_csv_files(group, op='avg')
##            ps.save_csv_file(filespath+name+'.csv', new)
##
##            new2 = ps.normalize_spec(new, kmer)
##            ps.save_csv_file(filespath+name+'-norm.csv', new2)
##
##            vals =list(zip(*new2.values()))[0]
##            stds =list(zip(*new2.values()))[1]
##            #xlab = [list(zip(*new2.keys()))[1]]
##
##            val_list.append(vals)
##            std_list.append(stds)
##            xlab_list.append(xlab)
##
##            #plotting individual plots
##            mpl.rc("savefig", dpi=320)
##            ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts, errorbars=[stds],
##                           titles=[titlename], ylabel='Proportion of mutations')
##            plt.savefig(picspath+name+'.svg')
##            plt.savefig(picspath+name+'.png')
##
##
##    if subtract_bkgr:
##
##        datafiles2 = datafiles[4:6]
##        titles2 = ['Untreated mice sum, n=10', 'NDMA-treated mice sum, n=10']
##
##        for group in datafiles2:
##            for idx, item in enumerate(group):
##                if item>100:
##                    group[idx]=filespath+str(item)+filesuff1
##                else:
##                    group[idx]=filespath+str(item)+filesuff2
##
##        print(datafiles2)
##
##        # Printing with purine labels
##        pu_dict=ps.init_spec_dict('purine')
##        xlab = list(zip(*pu_dict.keys()))[1]
##        
##        val_list=[]
##        std_list=[]
##        xlab_list=[]
##
##        sumfiles = ['all_ctr_sum', 'all_NDMA_sum']
##        sumspecs = {}
##        
##        for file, group in zip(sumfiles, datafiles2):
##            sumspecs[file] = ps.combine_csv_files(group, op='sum')
##            ps.save_csv_file(filespath+file+'.csv', sumspecs[file])
##
##        subspec = ps.subtract_background(sumspecs[sumfiles[1]], sumspecs[sumfiles[0]])
##        ps.save_csv_file(filespath+'all_NDMA_sum_bgsub.csv', subspec)
##
##        savename = 'all_NDMA_sum_bgsub_norm'
##
##        subspec2 = ps.normalize_spec(subspec, kmer)
##        ps.save_csv_file(filespath+savename+'.csv', subspec2)
##
##        #vals =list(zip(*subspec2.values()))[0]
##        vals = list(subspec2.values())
##        val_list.append(vals)
##        xlab_list.append(xlab)
##    
##        mpl.rc("savefig", dpi=320)
##        ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts,
##                           titles=['All NDMA bkgr sub n=10'], ylabel='Proportion of mutations')
##        plt.savefig(picspath+savename+'.svg')
##        plt.savefig(picspath+savename+'.png')
##
##
##
##        
##    #plotting composite figs
####    mpl.rc("savefig", dpi=600)
####    
####    ps.spec_figure(2, 2, val_list, xlabels=xlab_list, labels=ps.pu_muts, errorbars=std_list,
####                   titles=titlenames, ylabel='Proportion of mutations')
####
####    plt.savefig('fourplots-2x2-newcolor-hr.png')
####    plt.savefig('fourplots-2x2-newcolor-hr.svg')
##
##    #new = ps.read_csv_file(file)
##
##    #print(new)



if __name__=="__main__":
    main()

