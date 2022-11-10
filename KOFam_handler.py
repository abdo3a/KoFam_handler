#!/usr/bin/env python
 
import argparse
import subprocess
import os, sys
import pandas as pd
import numpy as np
import csv


parser = argparse.ArgumentParser(description='find the best KO in the KofamScan output and based ko list')

parser.add_argument('-ko', '--ko', help='KO list', required= True)
parser.add_argument('-o', '--outfolder', help='outfolder', required= True)

args = parser.parse_args()
ko_ls = args.ko
output = args.outfolder

if not os.path.exists(output):
    os.mkdir(os.path.join('./', output))
out_path = os.path.join('./', (output +'/'))

#read folder files function
def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


for file_n in files('./'):
    if file_n.endswith(".kofam"):
        infile_f1 = os.path.join('./', file_n.replace('.kofam','_f.tsv'))
        infile_f2 = os.path.join('./', file_n.replace('.kofam','_new.tsv'))
        outfile = os.path.join(out_path, file_n.replace('.kofam','.out'))
        subprocess.run("sed 's/\*/ /g;s/ KO definition/\tKO_definition/g;s/           /\t/g;s/          /\t/g;s/         /\t/g;s/        /\t/g;s/       /\t/g;s/      /\t/g;s/     /\t/g;s/    /\t/g;s/   /\t/g;s/  /\t/g;s/ /\t/;s/ /\t/;s/gene\tname/gene_name/g' %s >%s"% (file_n,infile_f1), shell=True, executable='/bin/bash')
        subprocess.run("paste <(cut -f2-4 %s|sed 's/-/0/g') <(cut -f5 %s) <(cut -f7-8 %s |sed 's/\t/ /g') > %s"% (infile_f1,infile_f1,infile_f1,infile_f2), shell=True, executable='/bin/bash')
        os.remove(infile_f1)

        read_file = pd.read_csv(infile_f2,sep='\t')
        df = read_file[["gene_name", "KO", "thrshld", "score", "KO_definition"]]
        df.drop([0], inplace=True)
        df['ration'] = ((df["score"]/(df['thrshld']/100))).round(2)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.dropna(subset=["ration"], how="all", inplace=True)

        with open(ko_ls) as f:
            lines = [line.rstrip('\n') for line in f]
        f =  open(outfile,"w+")
        for line in lines:
            if line == '':
                f.write(line + '\n')
            else:
                x = df[(df['KO'] == line)]
                x = x['ration'].max(axis=0)
                df2 = df[(df['KO'] == line) & (df['ration'] == x )]
                df2['ration'].values[df2['ration'].values >= 100] = 100
                if df2.empty:
                    f.write(line + ',' + 'NO' + '\n')
                else:
                    f.write(df2["KO"].to_string(index=False, header=False) + ',' + df2["ration"].to_string(index=False, header=False) + ',' + df2["gene_name"].to_string(index=False, header=False) + ',' + df2["KO_definition"].to_string(index=False, header=False) + '\n')
        f.close()
        

with open(os.path.join(out_path, 'script'),"w+") as mat:
    mat.write("sed -i 's/,/\t/g' %s\nfor i in `ls %s|head -1`; do cut -f1 $i;done > %s\npaste <(cut -f1 %s)"% (os.path.join(out_path, '*.out'), os.path.join(out_path, '*.out'), os.path.join(out_path, 'ko.txt'), os.path.join(out_path, 'ko.txt')))
    for file_n in files(out_path):
        if file_n.endswith(".out"):
            mat.write(" <(awk 'NR==1{print FILENAME; next}{print $2}' %s )"%(os.path.join(out_path, file_n)))
    mat.write(" >> %s\n"%(os.path.join(out_path, 'matrix.txt')))     
mat.close()


subprocess.run("chmod a+x %s && ./%s && sed  '3 i \n' %s|sed '/^n/s/n//g;s/.out//g;s/.\///g' >> %s"% (os.path.join(out_path, 'script'), os.path.join(out_path, 'script'), os.path.join(out_path, 'matrix.txt'), os.path.join(out_path, 'matrix.out')), shell=True, executable='/bin/bash')

os.remove(os.path.join(out_path, 'ko.txt'))
os.remove(os.path.join(out_path, 'matrix.txt'))
os.remove(os.path.join(out_path, 'script'))





