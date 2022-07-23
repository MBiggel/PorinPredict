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
    parser = argparse.ArgumentParser(usage = "PorinPredict -i </path/to/genome.fasta> -o </path/to/output/directory/> [other options]")
    parser.add_argument("-i", "--input", help = "Path to input genome in fasta format", nargs = 1, required = True)
    parser.add_argument("-o", "--outdir", help = "Path to desired outdir directory", nargs = 1, required = True)
    parser.add_argument("--summarize", help = "Summarize results in output folder.", action='store_true')
    parser.add_argument("-t", "--threads", help = "Number of CPU threads", type=int, default = 1)
    parser.add_argument("-v","--version", action="version", version='PorinPredict v' +  __version__, help="Print version")

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
    else:
        pass
        
    path = shutil.which("blastn")
    if path is None:
        logging.error('WARNING: Blastn could not be found')
        print()
        os._exit(0)
    else:
        pass
        
    path = shutil.which("Rscript")
    if path is None:
        logging.error('WARNING: Rscript could not be found')
        print()
        os._exit(0)
    else:
        pass
        
    # Check if input file path is available
    input = args.input[0]
    if os.path.isfile(input) == False:
        logging.error('WARNING: Assembly file not found. Check path and spelling.')
        print()
        os._exit(0)

    # get path to PorinPredict executable and databases
    PorinPredict_path = os.path.realpath(__file__)
    PorinPredict_path = os.path.dirname(PorinPredict_path) + "/"

    database_path_diamond = PorinPredict_path + r"/db/PA_OprD"
    database_path_blastn = PorinPredict_path + r"/db/PA_oprD_promoter_200bp_plus_10.fasta"
    
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
        
    # create temporary directory for CDS analysis
    if not os.path.isdir(outdir + "temp_dir_diamond"):
        os.makedirs(outdir + "temp_dir_diamond/")
    
    # create temporary directory for promoter analysis
    if not os.path.isdir(outdir + "temp_dir_blastn"):
        os.makedirs(outdir + "temp_dir_blastn/")
        
    # create directory for log files
    if not os.path.isdir(outdir + "logs"):
        os.makedirs(outdir + "logs/")    
    
    # path to results directory
    temp_dir_blastn = outdir + "temp_dir_blastn/"
    temp_dir_diamond = outdir + "temp_dir_diamond/"
    
    blast_result = temp_dir_blastn + prefix + '.tsv'
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
    run_blastn(input,blast_result,temp_dir_blastn,database_path_blastn,prefix,outdir,args.threads)
    run_R(blast_result,diamond_result,prefix,outdir,PorinPredict_path)
    
    if args.summarize == True:
        run_summarize(outdir,prefix)

    else:
        pass
    
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
            file.write("missing")          
    else:
        pass
    
    # checking if more than one porin detected
    df_diamond = pd.read_csv(diamond_result, sep = "\t", header=None)
    if df_diamond.count()[0] > 1:
        logging.warning("WARNING: More than one OprD detected, likely representing more than one P. aeruginosa strain in the sample. PorinPredict analyzes only the first hit.")
        print()
        
        # Drop additional hits by writing first row to file.
        df_diamond = df_diamond.loc[1:2]
    else:
        pass    
    
    #add column with prefix
    df_diamond["ID"] = [prefix]
    
    # write to tsv
    df_diamond.to_csv(diamond_result,sep='\t',index=False,header=False)
        
    return diamond_result
    
def run_blastn(input,blast_result,temp_dir_blastn,database_path_blastn,prefix,outdir,threads):
    
    blastn_command = ['blastn','-db',database_path_blastn,'-query',input,'-out',blast_result,'-perc_identity','95','-outfmt',"6 qseqid sseqid slen pident qcovs length mismatch gaps gapopen qstart qend sstart send evalue bitscore",'-max_target_seqs','10','-evalue','1E-20','-culling_limit','1','-num_threads',str(threads)]
    
    blastn = subprocess.run(blastn_command)
    
    logging.info("Running Blastn: " + str(blastn_command).replace(',',"").replace('\'',""))
    print()
    
    # checking if size of output is 0
    if os.stat(blast_result).st_size == 0:
        with open(blast_result, "w") as file:
            file.write("missing" + "\t" + prefix + "\n")
    else:        
        # keep hit with lowest Evalue
        out_blastn = pd.read_csv(blast_result, sep = "\t", header=None).sort_values(by = [13], ascending = True).head(1) # to be checked. make sure not surting by first digit or similar
        # add column with prefix
        out_blastn["ID"] = [prefix]
    
        # write to tsv
        out_blastn.to_csv(blast_result,sep='\t',index=False,header=False)
    
    return blast_result
    
def run_R(blast_result,diamond_result,prefix,outdir,PorinPredict_path):
    logging.info("Running PorinPredict.R")
    print()
    script_path = PorinPredict_path + "PorinPredict.R"
    args = [blast_result,diamond_result,prefix,outdir]
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