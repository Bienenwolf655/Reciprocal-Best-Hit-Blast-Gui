import sys
import subprocess 
packages = ['numpy', 'pandas', 'bio', 'tk']
for i in packages:
    subprocess.check_call(['python3', '-m', 'pip', 'install', f'{i}'])
import os
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import math
import os
from pathlib import Path


def main(files):
    files["effector_pred"] = files["effector_pred"].get()
    files["gene_list"] = files["gene_list"].get()
    files["exp"] = int(files["exp"].get())
    files["NumResBlast"] = files["NumResBlast"].get()
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
    threshold = math.exp(files["exp"])
    fwd_results = fwd_results[fwd_results["E-value"]<=threshold]
    rev_results = rev_results[rev_results["E-value"]<=threshold]

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
    rbbh = pd.read_csv(os.path.join(outdir,'rbbh.csv'))
    if files["gene_list"] == 1:
        data = list(rbbh['query_x'])
        with open(os.path.join(outdir,'genelist_query_x.txt'), 'w') as fn:
            fn.write('\n'.join(data))
            fn.close()
        data = list(rbbh['subject_x'])
        with open(os.path.join(outdir,'genelist_subject_x.txt'), 'w') as fn:
            fn.write('\n'.join(data))
            fn.close()
    rbbh.to_csv(os.path.join(outdir,'rbbh.csv'))
    if files["effector_pred"] == 1:
        s1_prot = list(rbbh['query_x'])
        s2_prot = list(rbbh['subject_x'])

        with open(s1,'r') as fn:
            line_s1 = fn.readlines() 
        with open(s2,'r') as fn:
            line_s2 = fn.readlines()
        f = open(f"{outdir}/rbbh.fasta", "a")
        for i in s1_prot:
            f.write(*[f"{line_s1[num]}{line_s1[num+1]}\n" for num,_ in enumerate(line_s1) if i in line_s1[num]])  

        for i in s2_prot:
            f.write(*[f"{line_s2[num]}{line_s2[num+1]}\n" for num,_ in enumerate(line_s2) if i in line_s2[num]])
        f.close()

        p = subprocess.Popen(["python", f"{sys.path[0]}/EffectorP-3.0/EffectorP.py","-f",  f"-i{outdir}/rbbh.fasta", 
        f"-o{outdir}/effektor.tab"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()        
        print(f"STDOUT:{stdout}\nSTDERR:{stderr}")
        df_eff = pd.read_csv(f"{outdir}/effektor.tab", sep='\t')
        eff = df_eff.loc[:, 'Cytoplasmic effector':'Prediction']
        final = pd.concat([rbbh, eff], axis=1, join="inner")
        final.to_csv(os.path.join(outdir,'final.csv'))
        messagebox.showinfo('Finished Job', f'Your joined effector prediction and RBBH output has been saved to {outdir}')
    else:
        messagebox.showinfo('Finished Job', f'Your RBBH output has been saved to {outdir}')
    del files


def app():
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
    
    window.geometry("800x400")
    
    window.config(background = "white")

    button_dict={}
    label_file_explorer_1 = Label(window,text = f"File Explorer S1",width = 100, height = 1,fg = "blue", highlightbackground='#3E4149')
    button_dict_1 =Button(window, text="S1", width=25, command=lambda*args: browseFiles("S1", label_file_explorer_1), highlightbackground='#3E4149')
    label_file_explorer_2 = Label(window,text = "File Explorer S2",width = 100, height = 1,fg = "blue", highlightbackground='#3E4149')
    button_dict_2 =Button(window, text="S2", width=25, command=lambda*args: browseFiles("S2", label_file_explorer_2), highlightbackground='#3E4149')
    label_file_explorer_3 = Label(window,text = "File Explorer Results Folder",width = 100, height = 1,fg = "blue", highlightbackground='#3E4149')
    button_dict_3 =Button(window, text="Results Folder", width=25, command=lambda*args: browseFiles("Results Folder", label_file_explorer_3), highlightbackground='#3E4149')

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

    drop = OptionMenu(window ,clicked, *options)
    label = ttk.Label(window,  text='Select the number of results to consider from the blast search:')
    spin = Spinbox(window, from_=-100, to=0, width=10, highlightbackground='#3E4149')
    label.pack()
    drop.pack()
    ttk.Label(window,  text='Select the exponent of the threshold e^ for the E-Value of the Blast Search:').pack()
    spin.pack()
    effector_pred = IntVar()
    gene_list = IntVar()
    Checkbutton(window, text="EffectorP 3.0 Prediction", variable=effector_pred, highlightbackground='#3E4149').pack()    
    Checkbutton(window, text="Gene List", variable=gene_list, highlightbackground='#3E4149').pack()
    files["gene_list"] = gene_list
    files["effector_pred"] = effector_pred
    files["exp"] = spin
    files["NumResBlast"] = clicked
    button_dict["Run"] = Button(window, text = "Run RBBH",command = lambda *args: main(files))
    button_dict["Run"].pack()
    button_dict["Exit"] = Button(window, text = "Exit",command = close)
    button_dict["Exit"].pack()
    window.mainloop()
    
if __name__ == "__main__":
    global files
    files = {}
    app()
    del files
