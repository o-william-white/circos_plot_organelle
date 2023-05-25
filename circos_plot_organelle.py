import argparse
import sys
import subprocess
import statistics
import os
import shutil

usage = """
 
"""

description = """
Python script to create a circos plot for an annotated organelle genome.
"""
# argparse
parser = argparse.ArgumentParser(usage=usage, description=description)
parser.add_argument("--fasta",  help = "Input fasta to plot", required=True)
parser.add_argument("--bed",    help = "BED file containing annotations", required=True)
parser.add_argument("--bam",    help = "BAM file for coverage estimation", required=True)
parser.add_argument("--window", help = "Window size for GC and repeat content", required=True, type=int)
args = parser.parse_args()

# check samtools in path
try:
    subprocess.call(["samtools"], stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print("Error: samtools not in path")
    sys.exit()

# create output directory
outdir = 'data'
if os.path.exists(outdir):
    print(f'Overwriting existing {outdir} directory')
    shutil.rmtree(outdir)
os.mkdir(outdir)

# function to read fasta
def read_fasta(filename):
    name, seq = None,""
    fasta = open(filename, "r")
    for line in fasta:
        if line.startswith(">") and name == None:
            name = line.rstrip("\n").replace(">","")
        else:
            if line.startswith(">") and name != None:
                yield [name, seq]
                name = line.rstrip("\n").replace(">","")
                seq = ""
            else:
                seq = seq + line.rstrip("\n")
    yield [name, seq]
    fasta.close()

# read fasta
f = list(read_fasta(args.fasta))
name = f[0][0]
sequence = f[0][1]
length = len(sequence)

# samtools depth
print("Running samtools coverage")
depth_list = []
for line in subprocess.check_output(["samtools", "depth",  "-a", args.bam]).decode().split("\n"):
    line = line.rstrip("\n").split("\t")
    if len(line) == 3:
        depth_list.append(line)

# mean depth per window
print(f"Calculating mean coverage in {args.window} bp windows")
cov = open(f"{outdir}/coverage.bed", "w")
for i in range(0, len(depth_list), args.window):
    start = i
    end   = i + args.window
    chunk = depth_list[start:end]
    # get mean coverage for window
    depth = []
    for x in chunk:
        name = x[0]
        depth.append(int(x[2]))
    m = statistics.mean(depth)
    # write to output
    cov.write(f"{name}\t{start}\t{end}\t{m}\n")
cov.close()

# function to calculate gc content
def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)
assert gc_content("ATCG") == 0.5

# gc content per window
print(f"Calculating GC content in {args.window} bp windows")
gc = open(f"{outdir}/gc.bed", "w")
for i in range(0, len(sequence), args.window):
    start = i
    end   = i + args.window
    chunk = sequence[start:end]
    g = round(gc_content(chunk), 2)
    gc.write(f"{name}\t{start}\t{end}\t{g}\n")
gc.close()

# function to calculate n content
def n_content(seq):
    return (seq.count("N") / len(seq))
assert n_content("NNCG") == 0.5

# trf
print(f"Running trf\n**********")
cmd_trf = f'trf {args.fasta} 2 7 7 80 10 50 500 -f -d -m -h'
subprocess.run(cmd_trf.split(' '))
print(f"\n\n**********\nFinished trf")

# read masked fasta
# the trf output is added to the present dir
fm = list(read_fasta(f"{os.path.basename(args.fasta)}.2.7.7.80.10.50.500.mask"))

# n content per window
name = fm[0][0]
sequence = fm[0][1]
print(f"Calculating N content in {args.window} bp windows")
rp = open(f"{outdir}/repeat.bed", "w")
for i in range(0, len(sequence), args.window):
    start = i
    end   = i + args.window
    chunk = sequence[start:end]
    n = round(n_content(chunk), 2)
    rp.write(f"{name}\t{start}\t{end}\t{n}\n")
rp.close()

# write karyotype file
print("Writing karyotype.txt")
with open("karyotype.txt", "w") as k:
    k.write(f"chr - {name} unused_label 0 {length} black\n")

# write ideogram conf
print("Writing ideogram.conf")
i_text = """
<ideogram>

<spacing>
default = 0u
break   = 0u
</spacing>

# Ideogram position, fill and outline
radius           = 0.90r
thickness        = 10p
fill             = yes
stroke_color     = black
stroke_thickness = 2p

# Minimum definition for ideogram labels.

show_label       = no
# see etc/fonts.conf for list of font names
label_font       = normal
label_size       = 30
label_parallel   = yes

</ideogram>

"""

with open("ideogram.conf", "w") as i:
    i.write(i_text)

# write ticks conf
print("Writing ticks.conf")
t_text = """
show_ticks          = yes
show_tick_labels    = yes


<ticks>

radius               = dims(ideogram,radius_outer)
colour               = black
label_multiplier     = 1e-3
orientation          = out
skip_first_label     = yes

<tick>
spacing = 100u
size = 10p
thickness = 2p
</tick>

<tick>
spacing = 1000u
size = 25p
thickness = 5p
show_label     = yes
label_size     = 40p
label_offset   = 10p
suffix         = " kb"
</tick>

</ticks>
"""

with open("ticks.conf", "w") as t:
    t.write(t_text)

# write circos conf
print("Writing circos.conf")
c_text = """
karyotype = karyotype.txt

<highlights>
<highlight>
file  = data/annotations_genes.txt
r1    = 0.90r
r0    = 0.85r
stroke_color = black
stroke_thickness = 2
</highlight>
</highlights>

<plots>

<plot>
type = text
file = data/annotations_labels.txt
r1   = 0.98r
r0   = 0.82r
label_font = normal
label_size = 30p
</plot>

<plot>
type = line
file = data/coverage.bed
r1   = 0.65r
r0   = 0.55r
fill_color = grey
thickness = 0p
extend_bin = no
</plot>

<plot>
type = line
file = data/gc.bed
r1   = 0.50r
r0   = 0.40r
fill_color = vdgrey
thickness = 0p
extend_bin = no
#orientation = in
</plot>

<plot>
type = line
file = data/repeat.bed
r1   = 0.35r
r0   = 0.25r
fill_color = vvdgrey
thickness = 0p
extend_bin = no
#orientation = in
</plot>

</plots>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>

<<include etc/housekeeping.conf>>
"""

with open("circos.conf", "w") as c:
    c.write(c_text)

# write annotations labels and genes file 
print("Writing annotations lables and gene files")
ann = open(f"{args.bed}", "r")
gen = open(f"{outdir}/annotations_genes.txt", "w")
lab = open(f"{outdir}/annotations_labels.txt", "w")

for line in ann: 
    line = line.rstrip("\n").split("\t")
    chr, start, end, gene, orien = line[0], line[1], line[2], line[3], line[5]
    col="blue"
    if gene.startswith("cox"):
        col="red"
    if gene.startswith("atp"):
        col="green"
    if gene.startswith("rrn"):
        col="yellow"
    if gene.startswith("cob"):
        col="vdgreen"
    if gene.startswith("trn"):
        col="purple"
    if orien == "+":
        r1 = "0.80r"
        r0 = "0.75r"
    else: 
        r1 = "0.75r"
        r0 = "0.70r"
    if int(start) < int(end):
        gen.write(f"{chr}\t{start}\t{end}\t{gene}\tfill_color={col},r0={r0},r1={r1}\n")
        lab.write(f"{chr}\t{start}\t{end}\t{gene}\n")
    
ann.close()
gen.close()
lab.close()

# run circos
print("Running circos\n**********\n")
cmd_circos = ["circos", "--conf", "circos.conf"]
subprocess.run(cmd_circos)
print(f"\n**********\nFinished circos")

# Complete
print("Complete!")


