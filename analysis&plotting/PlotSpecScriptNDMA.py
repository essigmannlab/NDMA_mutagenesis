#!/usr/bin/env python3
#
# PlotSpecScript
#
# Analyze data with PlotSpec functions.
#
#

import PlotSpec as ps
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

print(ps.pu_muts)

def main():
    #print("Hello there!")

    #control vars

    individual_spec=True
    process_avg=False
    subtract_bkgr=False   
    


    filespath = '../CSVfiles2/'
    picspath = '../Pics3/'
    filesuff1 = '_FINAL_unique.csv'
    filesuff2 = '_FINAL2_unique.csv'

    kmerfile = 'bbmap.count.EG10c.txt'
    kmer = ps.import_kmer_counts(kmerfile)

    filenames = ['males_ctr', 'females_ctr', 'males_NDMA', 'females_NDMA', 'all_ctr', 'all_NDMA']

    males_ctr = [354, 355, 358, 77, 78]
    females_ctr = [356, 357, 74, 75, 76]
    males_NDMA = [349, 350, 351, 80, 83]
    females_NDMA = [352, 353, 79, 81, 82]
    all_ctr =males_ctr+females_ctr
    all_NDMA=males_NDMA+females_NDMA


    datafiles = [males_ctr, females_ctr, males_NDMA, females_NDMA, all_ctr, all_NDMA]
    titlenames= ['Untreated males, n=5', 'Untreated females, n=5',
             'NDMA-treated males, n=5', 'NDMA-treated females, n=5',
             'Untreated mice, n=10', 'NDMA-treated mice, n=10']

    if individual_spec:

        #We need a direct plotting function in PlotSpec; supply csv file with data, and output a plot.

        


    if process_avg:

        #print(datafiles)

        for group in datafiles:
            for idx, item in enumerate(group):
                if item>100:
                    group[idx]=filespath+str(item)+filesuff1
                else:
                    group[idx]=filespath+str(item)+filesuff2

        print(datafiles)

        # Printing with purine labels
        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]
        
        val_list=[]
        std_list=[]
        xlab_list=[]

        for name,titlename, group in zip(filenames, titlenames, datafiles):
            new = ps.combine_csv_files(group, op='avg')
            ps.save_csv_file(filespath+name+'.csv', new)

            new2 = ps.normalize_spec(new, kmer)
            ps.save_csv_file(filespath+name+'-norm.csv', new2)

            vals =list(zip(*new2.values()))[0]
            stds =list(zip(*new2.values()))[1]
            #xlab = [list(zip(*new2.keys()))[1]]

            val_list.append(vals)
            std_list.append(stds)
            xlab_list.append(xlab)

            #plotting individual plots
            mpl.rc("savefig", dpi=320)
            ps.spec_figure(1, 1, [vals], xlabels=[xlab], labels=ps.pu_muts, errorbars=[stds],
                           titles=[titlename], ylabel='Proportion of mutations', colorscheme='COSMIC3')
            plt.savefig(picspath+name+'.svg')
            plt.savefig(picspath+name+'.png')


    if subtract_bkgr:

        datafiles2 = datafiles[4:6]
        titles2 = ['Untreated mice sum, n=10', 'NDMA-treated mice sum, n=10']

        for group in datafiles2:
            for idx, item in enumerate(group):
                if item>100:
                    group[idx]=filespath+str(item)+filesuff1
                else:
                    group[idx]=filespath+str(item)+filesuff2

        print(datafiles2)

        # Printing with purine labels
        pu_dict=ps.init_spec_dict('purine')
        xlab = list(zip(*pu_dict.keys()))[1]
        
        val_list=[]
        std_list=[]
        xlab_list=[]

        sumfiles = ['all_ctr_sum', 'all_NDMA_sum']
        sumspecs = {}
        
        for file, group in zip(sumfiles, datafiles2):
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
                           titles=['All NDMA bkgr sub n=10'], ylabel='Proportion of mutations')
        plt.savefig(picspath+savename+'.svg')
        plt.savefig(picspath+savename+'.png')



        
    #plotting composite figs
##    mpl.rc("savefig", dpi=600)
##    
##    ps.spec_figure(2, 2, val_list, xlabels=xlab_list, labels=ps.pu_muts, errorbars=std_list,
##                   titles=titlenames, ylabel='Proportion of mutations')
##
##    plt.savefig('fourplots-2x2-newcolor-hr.png')
##    plt.savefig('fourplots-2x2-newcolor-hr.svg')

    #new = ps.read_csv_file(file)

    #print(new)



if __name__=="__main__":
    main()

