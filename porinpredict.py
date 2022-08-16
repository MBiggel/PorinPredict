#!/usr/bin/env python3

"""
PorinPredict

Copyright 2022 Michael Biggel (michael.biggel@uzh.ch)

https://github.com/MBiggel/PorinPredict

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import os
import subprocess
import argparse
import datetime
import pandas as pd
import shutil
import logging
import datetime

__version__ = '1.0.0'

def get_arguments(args):
    parser = argparse.ArgumentParser(usage = "porinpredict -i </path/to/genome.fasta> -o </path/to/output/directory/> [other options]")
    parser.add_argument("-i", "--input", help = "Path to input genome in fasta format", nargs = 1, required = True)
    parser.add_argument("-o", "--outdir", help = "Path to desired outdir directory", nargs = 1, required = True)
    parser.add_argument("--summarize", help = "Summarize results in output folder.", action='store_true')
    parser.add_argument("-t", "--threads", help = "Number of CPU threads", type=int, default = 1)
    parser.add_argument("-v","--version", action="version", version='porinpredict v' +  __version__, help="Print version")

    args = parser.parse_args()
    return args

def main(args=None):
    args = get_arguments(args)
    print()
    
    # Check if dependencies available
    path = shutil.which("diamond")
    if path is None:
        logging.error('WARNING: Diamond could not be found')
        print()
        os._exit(0)
        
    path = shutil.which("blastn")
    if path is None:
        logging.error('WARNING: Blastn could not be found')
        print()
        os._exit(0)
        
    path = shutil.which("Rscript")
    if path is None:
        logging.error('WARNING: Rscript could not be found')
        print()
        os._exit(0)
        
    # Check if input file path is available
    input = args.input[0]
    if os.path.isfile(input) == False:
        logging.error('WARNING: Assembly file not found. Check path and spelling.')
        print()
        os._exit(0)

    # get path to PorinPredict executable and databases
    porinpredict_path = os.path.realpath(__file__)
    porinpredict_path = os.path.dirname(porinpredict_path) + "/"

    database_path_diamond = porinpredict_path + r"/db/PA_OprD"
    database_path_blastn_oprd = porinpredict_path + r"/db/PA_oprD.fasta"
    database_path_blastn_prom = porinpredict_path + r"/db/PA_oprD_promoter_200bp_plus_10.fasta"
    
    # get prefix and output arguments
    prefix = input.split("/")[-1].strip()
    prefix = ".".join(prefix.split(".")[0:-1])
    outdir = os.path.abspath(args.outdir[0])
    threads = args.threads

    # add a / to output directory name, if one is not already supplied
    if not outdir.endswith("/"):
        outdir = outdir.strip() + "/"

    # if the outdir directory doesn't yet exist, make it
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    # create directories for raw output and log files
    if not os.path.isdir(outdir + "temp_dir_diamond"):
        os.makedirs(outdir + "temp_dir_diamond/")
    if not os.path.isdir(outdir + "temp_dir_blastn_oprd"):
        os.makedirs(outdir + "temp_dir_blastn_oprd/")
    if not os.path.isdir(outdir + "temp_dir_blastn_prom"):
        os.makedirs(outdir + "temp_dir_blastn_prom/")
    if not os.path.isdir(outdir + "logs"):
        os.makedirs(outdir + "logs/")    
    
    # path to results directory
    temp_dir_blastn_prom = outdir + "temp_dir_blastn_prom/"
    temp_dir_blastn_oprd = outdir + "temp_dir_blastn_oprd/"
    temp_dir_diamond = outdir + "temp_dir_diamond/"
    
    blast_result_prom = temp_dir_blastn_prom + prefix + '.tsv'
    blast_result_oprd = temp_dir_blastn_oprd + prefix + '.tsv'
    diamond_result = temp_dir_diamond+prefix+'.tsv'
    
    # initialize log file
    logging.basicConfig(level = logging.DEBUG, filename = outdir + "logs/" + prefix + ".log", filemode = "a+", format = "%(message)s")
    logging.getLogger().addHandler(logging.StreamHandler())
         
    # run pipeline 
    now = datetime.datetime.now
    print("############")
    logging.info("PorinPredict " + __version__)
    logging.info(now().strftime("%Y-%m-%d %H:%M"))

    run_diamond(input,temp_dir_diamond,diamond_result,database_path_diamond,prefix,outdir,args.threads)
    run_blastn_prom(input,blast_result_prom,temp_dir_blastn_prom,database_path_blastn_prom,prefix,outdir,args.threads)
    run_blastn_oprd(input,blast_result_oprd,temp_dir_blastn_oprd,database_path_blastn_oprd,prefix,outdir,args.threads)
    run_R(blast_result_prom,blast_result_oprd,diamond_result,prefix,outdir,porinpredict_path)
    
    if args.summarize == True:
        run_summarize(outdir,prefix)
   
    logging.info("Results written to " + outdir + prefix + "_PorinPredict.tsv")
    print()
    
def run_diamond(input,temp_dir_diamond,diamond_result,database_path_diamond,prefix,outdir,threads):
    print()
    logging.info("Analyzing "+ prefix)
    print()
     
    diamond_command = ['diamond','blastx','--query-gencode','11','-d',database_path_diamond,'-q',input,'-o',diamond_result,'-p',str(threads),'--id','95','--subject-cover','60','--outfmt','6','qseqid','sseqid','slen','pident','scovhsp','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qstrand','qseq_translated','full_sseq','-b','6','--max-target-seqs','1','--sensitive','-c1','--masking','0','--quiet']
    
    subprocess.run(diamond_command)
    logging.info("Running DIAMOND: " + str(diamond_command).replace(',',"").replace('\'',""))
    print()
    
    # checking if size of output is 0
    if os.stat(diamond_result).st_size == 0:
        with open(diamond_result, "w") as file:
            file.write("no hit")          
    
    # checking if more than one porin detected; keep hit with lowest Evalue
    df_diamond = pd.read_csv(diamond_result, sep = "\t", header=None)
    if df_diamond.count()[0] > 1:
        logging.warning("WARNING: More than one OprD detected, likely representing more than one P. aeruginosa strain in the sample. PorinPredict analyzes only the first hit.")
        print()
        df_diamond = df_diamond.sort_values(by = [12], ascending = True).head(1)       
   
    #add column with prefix
    df_diamond["ID"] = [prefix]
    
    # write to tsv
    df_diamond.to_csv(diamond_result,sep='\t',index=False,header=False)
        
    return diamond_result
    
def run_blastn_prom(input,blast_result_prom,temp_dir_blastn_prom,database_path_blastn_prom,prefix,outdir,threads):
    
    blastn_prom_command = ['blastn','-db',database_path_blastn_prom,'-query',input,'-out',blast_result_prom,'-perc_identity','95','-outfmt',"6 qseqid sseqid slen pident qcovs length mismatch gaps gapopen qstart qend sstart send evalue bitscore",'-max_target_seqs','10','-evalue','1E-20','-culling_limit','1','-num_threads',str(threads)]
    blastn = subprocess.run(blastn_prom_command)
    
    logging.info("Running blastn promoter sequence: " + str(blastn_prom_command).replace(',',"").replace('\'',""))
    print()
    
    # checking if size of output is 0
    if os.stat(blast_result_prom).st_size == 0:
        with open(blast_result_prom, "w") as file:
            file.write("no hit" + "\t" + prefix + "\n")
    else:        
        # keep hit with lowest Evalue (highest bitscore if same Evalue)
        out_blastn = pd.read_csv(blast_result_prom, sep = "\t", header=None).sort_values(by = [14], ascending = False).sort_values(by = [13], ascending = True).head(1) 
        
        # add column with prefix
        out_blastn["ID"] = [prefix]
    
        # write to tsv
        out_blastn.to_csv(blast_result_prom,sep='\t',index=False,header=False)
   
    return blast_result_prom

def run_blastn_oprd(input,blast_result_oprd,temp_dir_blastn_oprd,database_path_blastn_oprd,prefix,outdir,threads):
    
    blastn_oprd_command = ['blastn','-db',database_path_blastn_oprd,'-query',input,'-out',blast_result_oprd,'-perc_identity','90','-outfmt',"6 qseqid sseqid slen pident qcovs length mismatch gaps gapopen qstart qend sstart send evalue bitscore",'-max_target_seqs','10','-evalue','1E-20','-culling_limit','1','-num_threads',str(threads)]
    blastn = subprocess.run(blastn_oprd_command)
    
    logging.info("Running blastn oprD: " + str(blastn_oprd_command).replace(',',"").replace('\'',""))
    print()
    
    # checking if size of output is 0
    if os.stat(blast_result_oprd).st_size == 0:
        with open(blast_result_oprd, "w") as file:
            file.write("no hit" + "\t" + prefix + "\n")
    else:        
        # keep hit with lowest Evalue (highest bitscore if same Evalue)
        out_blastn = pd.read_csv(blast_result_oprd, sep = "\t", header=None).sort_values(by = [14], ascending = False).sort_values(by = [13], ascending = True).head(1) 
        
        # add column with prefix
        out_blastn["ID"] = [prefix]
    
        # write to tsv
        out_blastn.to_csv(blast_result_oprd,sep='\t',index=False,header=False)
   
    return blast_result_oprd
    
def run_R(blast_result_prom,blast_result_oprd,diamond_result,prefix,outdir,porinpredict_path):
    logging.info("Running porinpredict.R")
    print()
    script_path = porinpredict_path + "porinpredict.R"
    args = [blast_result_prom,blast_result_oprd,diamond_result,prefix,outdir]
    cmd = ["Rscript", script_path] + args
    
    subprocess.run(cmd)      
       
def run_summarize(outdir,prefix):
    # create summary file from first output, if not existing 
    if not os.path.isfile(outdir + "PorinPredict_results_table.tsv"):
            shutil.copy(outdir + prefix + "_PorinPredict.tsv",outdir + "PorinPredict_results_table.tsv")
    # append results from additional output
    else:
        file1 = pd.read_csv(outdir + "PorinPredict_results_table.tsv", sep = "\t", index_col=None)
        file2 = pd.read_csv(outdir + prefix + "_PorinPredict.tsv", sep = "\t", index_col=None)
        cat = pd.concat([file1,file2])
        cat.to_csv(outdir + "PorinPredict_results_table.tsv", sep = "\t", index=False)   

if __name__ == "__main__":
    main()