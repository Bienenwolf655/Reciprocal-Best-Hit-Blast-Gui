import sys
import subprocess 
packages = ['numpy', 'pandas']
for i in packages:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', f'{i}'])

import os
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
from tkinter import ttk

def main(files):
    os.makedirs(files['Results Folder'], exist_ok=True)
    outdir = os.path.join(files['Results Folder'])
    s1 = os.path.join(files['S1'])
    s2 = os.path.join(files['S2'])

    fwd_out = os.path.join(outdir, 'fwd_results.tab')
    rev_out = os.path.join(outdir, 'rev_results.tab')


    if files["NumResBlast"] == "all":
        pass
    else:
        num = int(files["NumResBlast"])

    fwd_blastp = NcbiblastpCommandline(query=s1, subject=s2, out=fwd_out, max_target_seqs = num,
                                      outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                    )
    rev_blastp = NcbiblastpCommandline(query=s2, subject=s1, out=rev_out, max_target_seqs = num,
                                      outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                    )

    print("FORWARD: %s" % fwd_blastp)
    print("REVERSE: %s" % rev_blastp)

    fwd_stdout, fwd_stderr = fwd_blastp()
    rev_stdout, rev_stderr = rev_blastp()

    print("FWD STDOUT: %s" % fwd_stdout)
    print("FWD STDERR: %s" % fwd_stderr)
    print("REV STDOUT: %s" % rev_stdout)
    print("REV STDERR: %s" % rev_stderr)

    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    rev_results = pd.read_csv(rev_out, sep="\t", header=None)

    headers = ["query", "subject", "identity", "coverage",
            "qlength", "slength", "alength",
            "bitscore", "E-value"]
    fwd_results.columns = headers
    rev_results.columns = headers

    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    rev_results['qcov'] = rev_results.alength/rev_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    rev_results['scov'] = rev_results.alength/rev_results.slength

    threshold = 1
    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=threshold)
    rev_results['qcov'] = rev_results['qcov'].clip(upper=threshold)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=threshold)
    rev_results['scov'] = rev_results['scov'].clip(upper=threshold)

    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                  left_on='subject', right_on='query',
                  how='outer')

    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()

    rbbh.head()

    rbbh.to_csv(os.path.join(outdir,'rbbh.csv'))
    pb.stop()

from tkinter import *
from tkinter import filedialog

global files
files = {}
def browseFiles(name, label_file_explorer):
    
    if name == "S1" or name == "S2":
        filename = filedialog.askopenfilename(initialdir = "/home/",
                                              title = "Select a File",
                                              filetypes = (("Fasta files",
                                                            ".fasta"),
                                                           ("all files",
                                                            ".")))
    else:
        filename = filedialog.askdirectory()
    label_file_explorer.configure(text=f"File Opened: {filename}")
    files[name]=filename
    
def close():
    window.destroy() 
       

window = Tk()
window.title('Reciprocal Best Hit Blast')
  
window.geometry("800x300")
  
window.config(background = "white")
pb = ttk.Progressbar(window, orient='horizontal', mode='indeterminate', length=280)


button_dict={}
label_file_explorer_1 = Label(window,text = f"File Explorer S1",width = 100, height = 1,fg = "blue")
button_dict_1 =Button(window, text="S1", width=25, command=lambda*args: browseFiles("S1", label_file_explorer_1))
label_file_explorer_2 = Label(window,text = "File Explorer S2",width = 100, height = 1,fg = "blue")
button_dict_2 =Button(window, text="S2", width=25, command=lambda*args: browseFiles("S2", label_file_explorer_2))
label_file_explorer_3 = Label(window,text = "File Explorer Results Folder",width = 100, height = 1,fg = "blue")
button_dict_3 =Button(window, text="Results Folder", width=25, command=lambda*args: browseFiles("Results Folder", label_file_explorer_3))

label_file_explorer_1.pack()
button_dict_1.pack()
label_file_explorer_2.pack()
button_dict_2.pack()
label_file_explorer_3.pack()
button_dict_3.pack()
options = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "all"
]

clicked = StringVar()

clicked.set( "1" )

drop = OptionMenu(window ,clicked , *options )
label = ttk.Label(window,  text='Select the number of results to consider from the blast search:')
files["NumResBlast"] = clicked.get()
label.pack()
drop.pack()
pb.pack()
button_dict["Run"] = Button(window, text = "Run RBBH",command = lambda *args: main(files))
button_dict["Run"].pack()

button_dict["Exit"] = Button(window, text = "Exit",command = close)
button_dict["Exit"].pack()
window.mainloop()
del files
