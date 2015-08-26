
from __future__ import print_function
from __future__ import division

import numpy as np
import scipy as sp
import pdb
import h5py
import sys
import argparse
import logging
import time
import gzip
import re
import os
import math
import uuid
from collections import defaultdict

verboseprint=lambda *a, **k: None
__version__ = "1.0"

def main():

    parser=argparse.ArgumentParser(description='Extract c-row_stripe from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('-o', '--output', dest='outfile', type=str, help='interaction matrix output file')
    parser.add_argument('-b','--blocksize', dest='blocksize', type=int, default=32, help='block size of HDF5 file')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    infile=args.infile
    verbose=args.verbose
    outfile=args.outfile
    blocksize=args.blocksize
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("\n",end="")
    
    infile_name=os.path.basename(infile)
    infile_name=re.sub(".matrix", "", infile_name)
    infile_name=re.sub(".gz", "", infile_name)
    
    verboseprint("inputFileName",infile_name,sep="\t")
    
    verboseprint("")
    
    verboseprint("parsing matrix ...")
    headers,dim,assembly=get_matrix_info(infile)
    ncol=dim[0]
    nrow=dim[1]
    verboseprint("\tdim\t",nrow,"x",ncol,sep="")
    
    verboseprint("")
    
    # ensure symmetrical
    if nrow!=ncol:
        sys.exit('error: non-symmetrical matrix found [x,y] '+str(dim))
    n=nrow=ncol
    
    blocksize=min(max(blocksize,1),n)
    verboseprint("blocksize",blocksize,sep="\t")

    # calculate optimal block size
    itx_dtype=np.dtype(np.float64)
    itx_dtype_size=itx_dtype.itemsize
    
    verboseprint("assembly",assembly,sep="\t")
    
    infh=input_wrapper(infile)
    
    if outfile==None:
        outfile=infile_name+'.hdf5'
    elif not outfile.endswith('.hdf5') and not outfile.endswith('.h5'):
        outfile=outfile+'.hdf5'

    outhdf=h5py.File(outfile,'w')
    
    # attrs
    outhdf.attrs['genome']=assembly
    # row_stripe
    outhdf.create_dataset('interactions',shape=(n,n),dtype='float64',compression='gzip',chunks=(blocksize,blocksize))
    # my5C formatted headers
    outhdf.create_dataset('headers',data=headers)
    
    bin_chrs=[]
    bin_starts=[]
    bin_ends=[]
    chrs=[]
    row_stripe=np.zeros((blocksize,n),dtype='float64')
    row_stripe.fill(np.nan)
    
    c=0
    r=0
    
    verboseprint("")
    
    verboseprint("loading matrix ...")
    
    linenum=0
    num_zero=0
    # iterate over input file
    for line in infh:
        if line.startswith("#"): # skip headers
            continue
        if not line.lstrip(): # skip blank lines
            continue
        
        linenum+=1
        tmp=line.rstrip("\n").split("\t")

        if linenum == 1:
            continue
        
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',tmp[0])
        if m==None:
            sys.exit('error: incorrect input format!')
            
        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        
        # replace all "NA" with np.nan
        tmp[1:]=[np.nan if ((x=='NA') or (not is_float(x))) else float(x) for x in tmp[1:]]
        
        # once row_stripe is full, dump to hdf
        if (r != 0) and ((r%blocksize)==0):
            outhdf['interactions'][c:c+blocksize,:]=row_stripe
            
            # clear row_stripe
            row_stripe=np.zeros((min(blocksize,n-(c+blocksize)),n),dtype='float64')
            row_stripe.fill(np.nan)
            
            # increase row index
            c+=blocksize
            
        offset=r%blocksize 
        row_stripe[offset,:]=np.array(tmp[1:],dtype=np.float64)
        
        if len(chrs)==0 or chr_id!=chrs[-1]:
            chrs+=[chr_id]

        bin_chrs.append(len(chrs)-1)
        bin_starts.append(int(bin_start))
        bin_ends.append(int(bin_end))

        pc=(r/(n-1))*100
        verboseprint("\r\t"+str(r)+" / "+str(n-1    )+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
        if verbose: sys.stdout.flush()
        
        r += 1
        
    # handle final set of rows - if nrow not divisible by blocksize
    if (len(row_stripe)>0):
        outhdf['interactions'][c:c+row_stripe.shape[0],:]=row_stripe
    
    verboseprint("")
    verboseprint("")
    
    infh.close()

    chrs=np.array(chrs)
    nchrs=chrs.shape[0]
    
    bin_positions=np.c_[bin_chrs,bin_starts,bin_ends]
    
    chr_bin_range=np.zeros((nchrs,2))
  
    for i in xrange(nchrs):
        chr_bins=np.nonzero(bin_positions[:,0]==i)[0]
        chr_bin_range[i]=np.min(chr_bins),np.max(chr_bins)

    outhdf.create_dataset('chrs',data=chrs)
    outhdf.create_dataset('chr_bin_range',data=chr_bin_range,dtype='int64')
    outhdf.create_dataset('bin_positions',data=bin_positions,dtype='int64')
    
    if(outfile==None):
        outfile=re.sub(".matrix", "", outfile)
        outfile=re.sub(".gz", "", outfile)
    outfile=re.sub(".hdf5", "", infile_name)
    
def get_matrix_info(matrix_file):
  
    x_headers=[]
    y_headers=[]
    assembly=None
    num_zero=0
    num_nan=0
        
    infh=input_wrapper(matrix_file)
    
    linenum=0
    # iterate over input file
    for line in infh:
        if line.startswith("#"): # skip headers
            continue
        if not line.lstrip(): # skip blank lines
            continue
            
        linenum+=1
        tmp=line.rstrip("\n").split("\t")
        
        # first line must contain xHeaders
        if linenum == 1:
            # XDIMxYDIM
            dim = [None,None]
            if(len(tmp[0].split('x')) == 2):
                dim = [int(s) for s in tmp[0].split('x')]
            
            # required header format
            pattern = r'(\S+)\|(\S+)\|(\S+)\:(\d+)-(\d+)'
            
            # build bool list
            x_headers_check = (re.match('(\S+)\|(\S+)\|(\S+)\:(\d+)-(\d+)', x) for x in tmp[1:])
            bad_headers=[tmp[i] for i, x in enumerate(x_headers_check, 1) if not x]
            
            if len(bad_headers) > 0:
                print("Found improperly formatted x-headers!")
                print(bad_headers)
                sys.exit("")
            x_headers=tmp[1:]
            
            for m in x_headers_check:                
                bin_id,genome,chr_id,bin_start,bin_end=m.groups()
                if assembly == None:
                    assembly=genome
                elif genome != assembly:
                    sys.exit("genome ["+genome+" / "+assembly+"] is not constant!")
            continue
        
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',tmp[0])
        if m==None:
            print("Found improperly formatted y-headers!")
            print(tmp[0])
            sys.exit("")
        
        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        if assembly == None:
            assembly=genome
        elif genome != assembly:
            sys.exit("genome ["+genome+" / "+assembly+"] is not constant!")
            
        y_headers.append(tmp[0])
        
        tmp[1:]=[np.nan if ((x=='NA') or (not is_float(x))) else float(x) for x in tmp[1:]]
        num_zero += len([1 for x in tmp[1:] if x == 0])
        num_nan += tmp[1:].count(np.nan)
    
    # if matrix dimension is not embedded in tsv
    if dim[0] == None:
        dim[0] = len(x_headers)
    if dim[1] == None:
        dim[1] = len(y_headers)
            
    # ensure matrix dimensions match
    if dim[0] != len(x_headers):
        sys.exit("x-dimension != #xHeaders ["+str(dim[0])+" / "+str(len(x_headers))+"]")
    if dim[1] != len(y_headers):
        sys.exit("y-dimension != #yHeaders ["+str(dim[1])+" / "+str(len(y_headers))+"]")
    
    if not set(x_headers) & set(y_headers):
        sys.exit("xHeaders != yHeaders")
    
    headers=np.array(x_headers)
    
    verboseprint("\tnum_zero",num_zero)
    verboseprint("\tnum_nan",num_nan)
    
    return(headers,dim,assembly)

def input_wrapper(infile):
    if infile.endswith('.gz'):
        fh=gzip.open(infile,'r')
    else:
        fh=open(infile,'r')
        
    return fh
    
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def getSmallUniqueString():  
    tmp_uniq=str(uuid.uuid4())
    tmp_uniq=tmp_uniq.split('-')[-1]
    return(tmp_uniq)
    
def bin2header(bin,genome,chrs,index=getSmallUniqueString()):
    #name|assembly|chr:start-end
    header=str(index)+'|'+genome+'|'+str(chrs[bin[0]])+':'+str(bin[1])+'-'+str(bin[2])
    return(header)

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def deGroupHeader(header,extractBy="liteChr",index=getSmallUniqueString()):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect input format!')
                
    bin_id,genome,chr_id,bin_start,bin_end=m.groups()
    chr_id=chr_id.split('-')[0]

    header=str(bin_id)+'|'+genome+'|'+str(chr_id)+':'+str(bin_start)+'-'+str(bin_end)
    
    return(header)
    
def split_zoom_coord(z):
    """validate and split zoom coordinate.
    coordinates must be UCSC formatted.
    e.g. chr1:500-1000
    chr(colon)start(hyphen)end where start <= end
    """
    z=z.replace(',','')
    zoom_coord=re.search(r'(\S+):(\d+)-(\d+)',z)
    
    if zoom_coord==None:
        return None
        
    zoom_chr,zoom_start,zoom_end=zoom_coord.groups()
    zoom_start=int(zoom_start)
    zoom_end=int(zoom_end)
    
    if(zoom_start > zoom_end):
        return None
        
    return [zoom_chr,zoom_start,zoom_end]
        
def de_dupe_list(input):
    """de-dupe a list, preserving order.
    """
    
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

def byte_to_megabyte(byte):
    """convert bytes into megabytes.
    """
    byte=float(byte)
    return round(((byte / 1000) / 1000), 4) # megabyte
    #return round(((byte / 1024) / 1024),4) # mebibyte

    
def flip_intervals(a,b):
    """flip intervals, to ensure a < b
    """
    
    return(b,a)
    
def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    
    if a[0] < b[0] and a[1] > b[1]:
        return((b[1]-b[0])+1)
    
    if b[0] < a[0] and b[1] > a[1]:   
        return((a[1]-a[0])+1)
        
    if b[0] < a[0]:
        a,b=flip_intervals(a,b)
           
    return max(0, ( min(a[1],b[1]) - max(a[0],b[0]) ) ) 
    
if __name__=="__main__":
    main()

   