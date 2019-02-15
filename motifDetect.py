from optparse import OptionParser #modul for 
import inspect
import os
import sys
from pathlib import Path
import re # modul for regular expression
import random
import xml.etree.ElementTree as ET # creating xml output file
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker# for visualization
import matplotlib.lines as mlines
import numpy as np
from itertools import cycle#color  cycle for motifs
#import matplotlib.colors as pltc
from random import sample
from operator import itemgetter#für Darstellung =>Liste umsortieren


#Visualization:
global plt
global motifpos
fig = plt.figure(figsize=(13,6))
#plt.style.use('Solarize_Light2')
hax = fig.add_subplot(1,1,1)
#########################
#Initialization of variabels:
global proteins # dictionary with all sequences of proteins: header = key; values = sequences
global yseq
global all_legend
global seqlengths
global all_colors # list with all colors for the visualization of motifs
all_colors = ['aqua', 'blueviolet', 'darkred', 'darkslateblue', 'deeppink', 'dimgrey', 'dodgerblue', 'fuchsia', 'gold', 'green', 'grey', 'khaki', 'lightcoral', 'lime', 'mediumblue', 'mediumpurple', 'olive', 'plum', 'rebeccapurple', 'red', 'rosybrown', 'royalblue', 'sienna', 'steelblue', 'tan', 'tomato', 'turquoise', 'yellow']

all_legend =[]
yseq =[]
proteins = { "":""}
seqlengths = []
seqlengths.append(0)

tmp = inspect.getfile(inspect.currentframe())
expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
tmp = ""
##############################################################################
parser = OptionParser()
parser.add_option("-r", "--read", dest="file", default="None", help="Please specify your input sequence")#path of multifasta
parser.add_option("-s", "--show", dest="show",action = 'store_true', help="print all sequences of fasta file")
parser.add_option("-n", "--naiv", dest="naiv", default =None, help="Naive pattern matcher, plese give a pattern")
parser.add_option("-e", "--regularEx", dest="regular", default = None, help="Search of motifs via regular expression, please give a filename")
parser.add_option("-m", "--motif_score", dest="score", default=None, help="Finding motifs by rates and score, please give a filename")#filename with rates
parser.add_option("-p", "--pssm_score", dest="pssm", default=None, help="Finding motifs with pssm, please give a filename")#filename with pssm matrice
parser.add_option("-f", "--name", dest="namefile",default="motif_output.txt", help= "Please enter a name for your results")
#input : ein pfad für einen outputordner
(options, args) = parser.parse_args();
#############################################################################


############################################################################
#reading and splitting multi-fasta data + saving all sequences into a dictionary

def readFile(filename):
    fileObject = open(filename)
    fileText = fileObject.read()
    #print(fileText)
    splitMultiFastaSeqs (fileText)
    fileObject.close()
    return

def splitMultiFastaSeqs (fileText):
    #splits a multifasta data into single sequences
    curr_index = 0
    header_indices= []
    while curr_index != -1:# -1 because the method find() returns -1 when there is no more "<"
        curr_index = fileText.find(">")
        fileText = fileText[curr_index::]
        seq_startindex = fileText.find('\n') +1
        if curr_index == 0:
             header = fileText[curr_index+1 : seq_startindex-1]
        else:
            header = fileText[curr_index : seq_startindex-1]
        fileText = fileText[seq_startindex::]
        seq_endindex = fileText.find('\n') - 1
        seq = fileText[0: seq_endindex-1]
        fileText = fileText[seq_endindex + 1::]
        #print('seq', seq)
        saveSeq(seq, header)
    #out = open(output,'w')
    oview = open('overview_tab.txt','a+')
    oview.write('motif')#
    i=0
    for key in proteins.keys():
        oview.write(key+'\t')
        #out.write(key+'\t')
        i= i +1
    oview.write('\n')
    oview.close()
    return


def saveSeq (seq, header):
    #saves sequence and its key(header) in a dictionary
    newprotein = { header : seq}
    proteins.update(newprotein)
    return
######################################################################
# various methods of the program
#######################################################################
def parser(pattern):
    '''
    The following function looks for a strict consensus sequence in protein sequences.
    '''
    #search for a certain pattern   -naive search
    #outputfile = open(output,'a+')
    oview = open('overview_tab.txt','a+')
    oview.write(str('Naiveparser: '+str(pattern)+'\t'))
    #pattern = "DIFT
    txt = '\nNaivepar:'+str(pattern)+'\t'
    #outputfile.write(txt)
    i=0
    for key in proteins.keys():
        numb = 0
        seq =proteins[key]
        for j in range(len(seq)-len(pattern)+1):
            subseq = seq[j:j+ len(pattern)]
            if subseq == pattern:
                numb = numb +1
        oview.write(str(numb)+'\t')
    #outputfile.close()
    return
##########################################################################
#2nd method where the user define their motifs as regular expressions
def regularExpression(txtFile, tree):
    '''
    The following function gets a textfile containing patterns as regular expressions and
    looks for motifs in orthologs sequences.
    '''
    #txtFile is a textfile with regular expressions and their headers --> see readMe for format information
    # tree is 
    global oview
    #global outputfile
    #outputfile = open(output,'a+')
    oview = open('overview_tab.txt','a+')
    with open(txtFile) as f:
        content = f.readlines()
    motifcontent = False
    all_motifs ={"":""}
    for line in content:
        if line[0] == '>':
            motifheader= line[1: -1]
            motifcontent = True
            rates= []
        elif line[0] not in [ '>', '\n'] and motifcontent == True:
            if line[-1:] == '\n':
                line = line.replace('\n', '')# delete '\n' in motif sequence
            d = {motifheader: line}
            all_motifs.update(d)
    del all_motifs[""]#dictionary with all defined motifs from the external file
    cycol = cycle(sample(all_colors, len(all_motifs)))
    mcolor =[]#saves the chosen color for every motif (for plotting)
    motifpos = []
    for l in range(len(proteins)):
        motifpos.append([])
    all_legend = []
    for motif in all_motifs:
        #outputfile.write('\nregEx:'+motif+'\t')
        oview.write('\nregEx:'+motif+'\t')
        currentcolor = next(cycol)
        mcolor.append(currentcolor)
        legend_line = mlines.Line2D([], [], color=currentcolor, marker='.', markersize=15, label=motif)
        all_legend.append(legend_line)
    plt.legend(handles = all_legend,bbox_to_anchor=(1, 1), loc=4, borderaxespad=0.) 
    numberOfSeq = 0
    for key in proteins.keys():
        seq =proteins[key]#seq is single sequence of the multifasta.data
        root = tree.getroot()
        items = ET.SubElement(root, 'protein', id =  str(key), length = str(len(seq)))
        plt.hlines(y=numberOfSeq+1, xmin=0, xmax= len(seq), lw =1)
        yseq.append((numberOfSeq+1, str(key), len(seq)))
        c = 0#index
        for motif in all_motifs:
            #outputfile.write('\nregEx:'+motif+'\t')
            item1 = ET.SubElement(items, 'feature', type = motif )
            motif_num = re.findall(all_motifs[motif], seq)#length of this shows the number of occurence
            motifIndex = [(mo.start(0), mo.end(0)) for mo in re.finditer(all_motifs[motif], seq)]
            #print('Index', motifIndex)
            #https://docs.python.org/2/library/re.html  ----for more informations to build up regular expressions
            if motif_num != []:
                mode = True
                #outputfile.write('YES\t')
                oview.write(str(len(motif_num))+'\t')
                for mtuple in motifIndex:
                    start = ET.SubElement(item1,'start', start = str(mtuple[0]))
                    end = ET.SubElement(item1, 'end', end = str(mtuple[1]))
                    item1.set('instance', str(len(motif_num)))
                    color = mcolor[c]
                    onemotif = (numberOfSeq, mtuple[0] , mtuple[1], color)
                    motifpos[numberOfSeq].append(onemotif)
            else:
                #outputfile.write('-\t')
                oview.write('0\t')
            c = c+1
        oview.write('\n')
        numberOfSeq = numberOfSeq +1
    #outputfile.close() 
    return (tree, motifpos)
###########################################################################
def motif_scoreFinder(txtFile):
    global rates
    #fileObject = open(txtFile)
    motif_headers = []
    with open(txtFile) as f:
        content = f.readlines()
    motifcontent = False
    all_rates =[]
    for line in content:
        if line[0] == '>':
            motifheader= line[1: -1]
            motif_headers.append(motifheader)
            motifcontent = True
            rates= []
        elif line[0] not in [ '>', '\n'] and motifcontent == True:
            tab_index = [m.start() for m in re.finditer('\t', line)]#
            amino_index= [ x+1 for x in tab_index ]
            all_dic = {"":""}
            for amino in amino_index:
                d = {line[amino]: float(line[amino+2:amino+5])}
                all_dic.update(d)
            del all_dic[""]
            rates. append(all_dic)
            pass
        elif line[0] == '\n':
            motifcontent = False
            all_rates.append(rates)# length of this list is the number of motifs
    all_rates.append(rates)
    return [all_rates, motif_headers]


def pssm_scoreFinder(txtFile):
    global rates
    #fileObject = open(txtFile)
    with open(txtFile) as f:
        content = f.readlines()
    #print(content)
    motifcontent = False
    rates =[]
    coSeq =""
    count = 0
    for line in content:
        if count == 0:
            amino = line.split('\t')
            if amino[2] == 'Master':
                mode = 'cdd'
                amino= amino[3:-1]
            else:
                mode= 'normal'
                amino= amino[2:-1]
        else:
            all_dic = {"":""}
            line = line.split('\t')
            coSeq=coSeq+line[1]
            if mode == 'cdd':
                line= line[3:-1]
            elif mode == 'normal':
                line= line[2:-1]
            else:
                print(' Your pssm has an invalid format')
            i =0
            for element in line:
                d = {amino[i]: int(element)}
                all_dic.update(d)
                i =i+ 1
            del all_dic[""]
            rates. append(all_dic)
        count = count +1
    return [rates, coSeq]


def rate(seq, rates):
    i = 0
    score = 0
    for amino in seq:
        currentRates = rates[i]
        if amino in currentRates:
            score = score + currentRates[amino]
        i = i+1
    return score

def cutoffFinder(number_randomSeqs, rates, probability):
    alphabet_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", \
            "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    lengthForSeq = len(rates)
    all_scores =[]
    randomSeq =""
    for i in range(number_randomSeqs):
        for j in range(lengthForSeq):
            randomSeq = random.choice(alphabet_aa)+randomSeq
        #print(randomSeq)
        all_scores.append(rate(randomSeq, rates))
        randomSeq = ""
    all_scores.sort()
    all_scores.reverse()
    #print(all_scores)
    cutoff = all_scores[int(probability/100 * len(all_scores))]
    #liste abnehmend sortieren
    #probability/100 * len(liste) - zeigt wie vielte elemente der cutoff ist
    return cutoff

def seqParser(cutoff, rates, outputtxt, motifname, tree, a, motifpos, color):
    #looks for motif in all sequences
    global oview
    outputfile = open('motifRateresults.txt','a+')
    outputfile.write('\nCurrent motif: '+motifname+'\n')
    oview = open('overview_tab.txt','a+')
    oview.write("\nmotifRate:"+motifname+"\t")
    lengthForMotif = len(rates)
    numberOfSeq = 0
    for key in proteins.keys():
        output = str(numberOfSeq+1)+" protein: "+str(key)# 
        outputfile.write(output)
        currentSeq = proteins[key]
        root = tree.getroot()
        if a == 1:
            items = ET.SubElement(root, 'protein', id =  str(key), length = str(len(currentSeq)))
            plt.hlines(y=numberOfSeq+1, xmin=0, xmax= len(currentSeq), lw =1)
            #plt.plot([0, len(currentSeq)], [numberOfSeq+1,numberOfSeq+1], 'k-|')
            yseq.append((numberOfSeq+1, str(key), len(currentSeq)))
            for l in range(len(proteins)):
                motifpos.append([])
        else:
            for proteinnode in root:#searching for the proteinnode
                proteinID = proteinnode.attrib['id']
                if proteinID == str(key):
                    items = proteinnode
        item1 = ET.SubElement(items, 'feature', type = motifname )# feature node
        output = "  length: "+ str(len(currentSeq))+'\n'
        outputfile.write(output)
        instances = 0
        for j in range(len(currentSeq)-lengthForMotif+1):
            subseq = currentSeq[j:j+ lengthForMotif]
            currentscore = rate(subseq, rates)
            if currentscore >= cutoff:
                value1 = str( j+1)
                value2 = str(j+1+lengthForMotif)
                output = "Residue: "+value1+' - '+value2+'\n'
                outputfile.write(output)
                start = ET.SubElement(item1,'start',  start = value1)
                end = ET.SubElement(item1, 'end', end = value2)
                onemotif = (numberOfSeq, j, j+lengthForMotif, color)#for plotting: saves a tuple with (number of current sequnce, start position, end position)
                #plt.hlines(y=numberOfSeq+1, xmin=value1, xmax= value2,  label ='label', lw =10)
                #plt.hlines(colors= color, label = motifname)
                instances = instances +1#number of occurence
                #newnode.set( 'start' ,value1)
                motifpos[numberOfSeq].append(onemotif)
        item1.set('instance', str(instances))
        oview.write(str(instances)+'\t')
        numberOfSeq = numberOfSeq + 1
    instances = 0
    return tree

def xmlmaker (treeroot):
    '''
    The following function gets a root from an ElementTree-Object and creats a xml-file.
    '''
    #create a new XML file 
    treefile = ET.tostring(treeroot, encoding='utf8', method='xml') 
    myfile = open("outputmotifs.xml", "wb")  
    myfile.write(treefile)
    
############################################################################
# The next part includes the handling of options and plotting functions.
############################################################################
def Main():
    global root
    global tree
    overview = open('overview_tab.txt','w')
    filename = options.file
    readFile(filename)#, options.namefile)
    del proteins[""]
#####Jobs:
    def plotter(positions):
        '''
        The following function arranges all motifs in a plot.
        '''
        # The list "positons" holds information of motifs for every proteinsequence.
        seqno = 0#number of proteinsequence
        print(positions)
        for key in proteins.keys():
            seq =proteins[key]#sequence as string of current protein sequence
            seqlengths.append(len(seq))
            seqinfo= positions[seqno]# information about motifs in current sequence
            seqinfo = sorted(seqinfo,key=itemgetter(1))
            already_occ = [[] for i in range(len(seq))]
            def appender( start, end, occList, element):
                for q in range(start, end):
                    occList[q].append(element)
                return occList
            for p in seqinfo:
                currentrow = 0
                while currentrow in already_occ[p[1]]:
                    currentrow = currentrow +1
                else:
                    plt.hlines(y=p[0]+1+ 0.05 * currentrow, xmin=p[1], xmax=p[2], colors= p[3],  lw =5)
                    already_occ = appender(p[1], p[2], already_occ, currentrow)
            seqno = seqno +1
        return
    
    def sc_jobs(root, tree, all_legend):
        '''
        The following function executes all nessecary function and creates all file
        for the scoreMotif option.
        '''
        r_path = options.score
        results = motif_scoreFinder(r_path)
        allrates = results[0]
        motifs_header = results[1]
        n = 0
        m = 0 # Zählindex für tree
        open('motifRateresults.txt','w')
        motifpos =[]# saves the position of motifs for plotting
        all_legend =[]
        cycol = cycle(sample(all_colors, len(allrates)))
        for rates in allrates:
            cutoff =cutoffFinder(10000, rates, 0.2)
            m= m+1
            mcolor = next(cycol)
            tree = seqParser(cutoff, rates, options.namefile, motifs_header[n], tree, m, motifpos, mcolor)
            legend_line = mlines.Line2D([], [], color=mcolor, marker='.', markersize=15, label= motifs_header[n])
            all_legend.append(legend_line)
            n= n+1
        plt.legend(handles = all_legend, bbox_to_anchor=(1, 1), loc=4, borderaxespad=0.) 
        root = tree.getroot()
        #tree.write("outputmotifs.xml")#create a new XML file with the results (update)
        xmlmaker(root)
        plotter(motifpos)

    def pssm_jobs(root, tree, all_legend):#job scoreMotif
        '''
        The following function executes all nessecary function and creates all file
        for the pssm option.
        '''
        r_path = options.pssm
        results = pssm_scoreFinder(r_path)
        rates = results[0]
        motifs_header = results[1]
        n = 0
        m = 0 # Zählindex für tree
        #open('motifpssmresults.txt','w')
        motifpos =[]# saves the position of motifs for plotting
        all_legend =[]
        cutoff =cutoffFinder(10000, rates, 0.2)
        m= m+1
        mcolor = sample(all_colors, 1)
        tree = seqParser(cutoff, rates, options.namefile, motifs_header, tree, m, motifpos, mcolor[0])
        legend_line = mlines.Line2D([], [], color=mcolor[0], marker='.', markersize=15, label= motifs_header)
        all_legend.append(legend_line)
        plt.legend(handles = all_legend, bbox_to_anchor=(1, 1), loc=4, borderaxespad=3.) 
        root = tree.getroot()
        #tree.write("outputmotifs.xml")#create a new XML file with the results (update)
        xmlmaker(root)
        plotter(motifpos)
        
#######################################        
# Handling options and executing right function #
#######################################
    if options.show:
        print(proteins)
    if options.naiv != None:
        parser(options.naiv)#, options.namefile)
    if options.regular != None or options.score != None:#at least one of the options
        root = ET.Element('root', name = 'custom')# tree
        xmlmaker(root)#create a new XML file
        tree = ET.parse('outputmotifs.xml')
        if options.regular != None and options.score == None:#just regular-option
            ex_path = options.regular
            results = regularExpression(ex_path, tree)#, options.namefile)
            root = tree.getroot()
            xmlmaker(root)
            plotter(results[1])
        elif options.regular == None and options.score != None:#just score option
            sc_jobs(tree.getroot(), tree, all_legend)
        elif options.regular == None and options.score == None:#both options
            ex_path = options.regular
            tree = regularExpression(ex_path, tree)#, options.namefile)
            sc_jobs(tree.getroot(), tree, all_legend)
    if options.pssm != None:
        root = ET.Element('root', name = 'custom')# tree
        xmlmaker(root)#create a new XML file
        tree = ET.parse('outputmotifs.xml')
        pssm_jobs(tree.getroot(), tree, all_legend)
        


Main()
hax.set_xlabel('position')
hax.set_ylabel('sequences')
hax.set_title('Features in sequences')


seqlengths.append(seqlengths[-1]+10)
steps = round((seqlengths[-1]+10) / 20)
x_axis = np.arange(0, seqlengths[-1]+10, steps)
hax.set_xticks(x_axis)
#hax.set_xticks(np.arange(0, 1000+20, 5))
# yseq.append((numberOfSeq+1, str(key), len(currentSeq)))
yheaders =[]
ypos = []
for seqY in yseq:# for labelling sequenques
    hax.text(seqY[2], seqY[0], str(seqY[2]), fontsize=10)
    hax.text(-2.5, seqY[0], '0', fontsize=10)
    yheaders.append(seqY[1])
    ypos.append(seqY[0])
#print(yseq)
hax.set_yticks(ypos)
hax.set_yticklabels(yheaders)
plt.show()
