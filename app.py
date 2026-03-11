import tkinter as tk
from tkinter import scrolledtext, filedialog
import re

# =========================
# Genetic code
# =========================
genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# =========================
# Reverse complement
# =========================
def reverse_complement(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    return "".join(comp[b] for b in seq[::-1])

# =========================
# Translate DNA to Protein
# =========================
def translate_dna(seq, start_codon="ATG"):
    protein=""
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            if i == 0 and codon == start_codon:
                aa = 'M'
            else:
                aa = genetic_code.get(codon, 'X')
            protein += aa
            if aa == '*':
                break
    return protein

# =========================
# Promoter search
# =========================
def similarity(a,b):
    return sum([1 for x,y in zip(a,b) if x==y])

def search_promoters(seq, box10="TATAAT", box35="TTGACA", threshold=5):
    pos10=[]
    pos35=[]
    for i in range(len(seq)-5):
        frag = seq[i:i+6]
        if similarity(frag, box10) >= threshold:
            pos10.append(i)
        if similarity(frag, box35) >= threshold:
            pos35.append(i)
    return pos10, pos35

# =========================
# Rho-independent terminator
# =========================
def is_palindrome(seq):
    return seq == reverse_complement(seq)

def search_terminator(seq):
    positions=[]
    for i in range(len(seq)-20):
        fragment=seq[i:i+8]
        if is_palindrome(fragment):
            if "TTTT" in seq[i+8:i+20]:
                positions.append(i)
    return positions

# =========================
# Shine-Dalgarno
# =========================
def search_SD(seq):
    positions=[]
    for m in re.finditer("AGGAGG", seq):
        positions.append(m.start())
    return positions

# =========================
# Coloring helper
# =========================
def color_sequence(text_widget, positions, length, color):
    for pos in positions:
        start=f"1.{pos}"
        end=f"1.{pos+length}"
        tag=f"{pos}{color}"
        text_widget.tag_add(tag, start, end)
        text_widget.tag_config(tag, foreground=color)

# =========================
# Load FASTA
# =========================
def load_fasta():
    file = filedialog.askopenfilename()
    if file:
        with open(file) as f:
            lines = f.readlines()
        seq = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
        input_box.delete("1.0", tk.END)
        input_box.insert(tk.END, seq)
        seq_length_var.set(str(len(seq)))

# =========================
# Display codons in frames
# =========================
def display_all_codons(seq, text_widget, strand="Forward"):
    text_widget.delete("1.0", tk.END)
    text_widget.insert(tk.END, f"{strand} Strand Frames:\n\n")
    for frame in range(3):
        codons = [seq[i:i+3] for i in range(frame, len(seq)-2, 3)]
        spaced = " ".join(codons)
        text_widget.insert(tk.END, f"Frame {frame+1}:\n{spaced}\n\n")

# =========================
# ORF detection
# =========================
def find_orfs(seq, start_codon="ATG", stop_codons=None):
    if stop_codons is None:
        stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq)-2:
            codon = seq[i:i+3]
            if codon == start_codon:
                j = i+3
                while j < len(seq)-2:
                    stop = seq[j:j+3]
                    if stop in stop_codons:
                        orfs.append((i, j+3, j+3-i, frame))
                        break
                    j +=3
                i = j
            else:
                i +=3
    return orfs

# =========================
# LONGEST ORF
# =========================
def find_longest_orf(seq, start_codon="ATG", stop_codons=None):
    orfs = find_orfs(seq, start_codon, stop_codons)
    if not orfs:
        return None
    longest = max(orfs, key=lambda x: x[2])
    start, stop, length, frame = longest
    dna_seq = seq[start:stop]
    protein_seq = translate_dna(dna_seq, start_codon)
    codons_spaced = " ".join([dna_seq[i:i+3] for i in range(0, len(dna_seq),3)])
    return {
        'start': start, 'stop': stop, 'length': length,
        'frame': frame, 'dna': dna_seq,
        'codons_spaced': codons_spaced,
        'protein': protein_seq,
        'start_codon': start_codon
    }

# =========================
# Display LONGEST ORF
# =========================
def display_longest_orf(orf_info, text_widget, strand="Forward"):
    text_widget.delete("1.0", tk.END)
    if orf_info is None:
        text_widget.insert(tk.END, f"No ORF found on {strand}")
        return
    text_widget.insert(tk.END, f"=== LONGEST ORF - {strand.upper()} ===\n\n")
    text_widget.insert(tk.END, f"Position: {orf_info['start']} - {orf_info['stop']}\n")
    text_widget.insert(tk.END, f"Frame: {orf_info['frame']+1}\n")
    text_widget.insert(tk.END, f"Length: {orf_info['length']} bp ({len(orf_info['protein'])-1} codons)\n")
    start_codon = orf_info['start_codon']
    stop_codon = orf_info['dna'][-3:]
    text_widget.insert(tk.END, f"Start Codon: {start_codon} (position {orf_info['start']})\n")
    text_widget.insert(tk.END, f"Stop Codon: {stop_codon} (position {orf_info['stop']-3})\n\n")
    text_widget.insert(tk.END, "DNA Sequence:\n")
    text_widget.insert(tk.END, f"{orf_info['dna']}\n\n")
    text_widget.insert(tk.END, "DNA (spaced by codons):\n")
    text_widget.insert(tk.END, f"{orf_info['codons_spaced']}\n\n")
    num_aa = len(orf_info['protein'])-1
    text_widget.insert(tk.END, f"Protein Sequence ({num_aa} amino acids):\n")
    text_widget.insert(tk.END, f"{orf_info['protein']}\n")

# =========================
# Analyze sequence
# =========================
def analyze():
    dna = input_box.get("1.0", tk.END).strip().upper()
    if not dna:
        return
    seq_length_var.set(str(len(dna)))

    # ORFs forward and reverse
    display_all_codons(dna, orf_box_fwd, "Forward")
    rev_seq = reverse_complement(dna)
    display_all_codons(rev_seq, orf_box_rev, "Reverse")

    longest_fwd = find_longest_orf(dna, "ATG", ["TAA","TAG","TGA"])
    display_longest_orf(longest_fwd, longest_orf_fwd_box, "Forward Strand")
    longest_rev = find_longest_orf(rev_seq, "TAC", ["ATT","ATC","ACT"])
    display_longest_orf(longest_rev, longest_orf_rev_box, "Reverse Strand")

    # Promoters
    box10, box35 = search_promoters(dna, "TATAAT", "TTGACA")
    box10r, box35r = search_promoters(rev_seq, "ATTATA", "TGTCAA")  # reverse complements exactes

    # Terminators
    rho = search_terminator(dna)

    # Shine-Dalgarno
    sd_pos = search_SD(dna)

    # ORFs start/stop for coloring
    orfs = find_orfs(dna)

    # Clear bars and insert sequence
    for bar in bars:
        bar.delete("1.0", tk.END)
        bar.insert(tk.END, dna)

    # Color features
    for start, stop, length, frame in orfs:
        color_sequence(bars[0], [start], 3, "green")      # start codon
        color_sequence(bars[0], [stop-3], 3, "red")      # stop codon

    color_sequence(bars[1], box10, 6, "green")
    color_sequence(bars[2], box35, 6, "red")
    color_sequence(bars[3], box10r, 6, "lightblue")
    color_sequence(bars[4], box35r, 6, "plum")
    color_sequence(bars[5], rho, 8, "brown")
    color_sequence(bars[6], sd_pos, 6, "orange")

# =========================
# GUI
# =========================
window = tk.Tk()
window.title("DNA Feature Analyzer")
window.geometry("1000x1000")

canvas = tk.Canvas(window)
scrollbar = tk.Scrollbar(window, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=scrollbar.set)
scrollbar.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)

frame_inside = tk.Frame(canvas)
canvas.create_window((0,0), window=frame_inside, anchor="nw")

def on_frame_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

frame_inside.bind("<Configure>", on_frame_configure)

# Sequence length
length_frame = tk.Frame(frame_inside)
length_frame.pack(fill="x", padx=10, pady=5)
tk.Label(length_frame, text="Sequence Length:").pack(side="left")
seq_length_var = tk.StringVar(value="0")
seq_length_entry = tk.Entry(length_frame, textvariable=seq_length_var, width=10, state="readonly")
seq_length_entry.pack(side="left", padx=5)

# DNA Input
input_frame = tk.LabelFrame(frame_inside, text="DNA Input")
input_frame.pack(fill="x", padx=10, pady=5)
input_box = scrolledtext.ScrolledText(input_frame, height=6, font=("Courier",12))
input_box.pack(fill="x")

# Buttons
btn_frame = tk.Frame(frame_inside)
btn_frame.pack(pady=5)
tk.Button(btn_frame, text="Load FASTA", command=load_fasta).pack(side="left", padx=5)
tk.Button(btn_frame, text="Analyze", command=analyze).pack(side="left", padx=5)

# LONGEST ORFs
longest_orf_fwd_frame = tk.LabelFrame(frame_inside, text="LONGEST ORF - Forward Strand", bg="lightgreen")
longest_orf_fwd_frame.pack(fill="x", padx=10, pady=5)
longest_orf_fwd_box = scrolledtext.ScrolledText(longest_orf_fwd_frame, height=12, font=("Courier",11), bg="honeydew")
longest_orf_fwd_box.pack(fill="x")

longest_orf_rev_frame = tk.LabelFrame(frame_inside, text="LONGEST ORF - Reverse Strand", bg="lightblue")
longest_orf_rev_frame.pack(fill="x", padx=10, pady=5)
longest_orf_rev_box = scrolledtext.ScrolledText(longest_orf_rev_frame, height=12, font=("Courier",11), bg="azure")
longest_orf_rev_box.pack(fill="x")

# Frames codons
orf_frame_fwd = tk.LabelFrame(frame_inside, text="Forward Frames")
orf_frame_fwd.pack(fill="x", padx=10, pady=5)
orf_box_fwd = scrolledtext.ScrolledText(orf_frame_fwd, height=10, font=("Courier",12))
orf_box_fwd.pack(fill="x")

orf_frame_rev = tk.LabelFrame(frame_inside, text="Reverse Frames")
orf_frame_rev.pack(fill="x", padx=10, pady=5)
orf_box_rev = scrolledtext.ScrolledText(orf_frame_rev, height=10, font=("Courier",12))
orf_box_rev.pack(fill="x")

# Bars for coloring
bars = []
titles = ["Start/Stop codons", "Promoter -10 Forward", "Promoter -35 Forward",
          "Promoter -10 Reverse", "Promoter -35 Reverse", "Rho-independent terminator", "Shine-Dalgarno"]

for t in titles:
    frame = tk.LabelFrame(frame_inside, text=t)
    frame.pack(fill="x", padx=10, pady=5)
    bar = scrolledtext.ScrolledText(frame, height=3, font=("Courier",12))
    bar.pack(fill="x")
    bars.append(bar)

window.mainloop()