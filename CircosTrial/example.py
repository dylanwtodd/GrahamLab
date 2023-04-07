from pycirclize import Circos
from pycirclize.parser import Gff
import pandas as pd
import numpy as np

# Load GFF file

gff_file = 'dromelmitofeatures.gff'
gff = Gff(gff_file)

circos = Circos(sectors={gff.name: gff.range_size})
circos.text("Drosophila Melanogaster Mitochondrial Genome", size=8)

sector = circos.sectors[0]

cds_track = sector.add_track((90, 100))
cds_track.axis(fc="#EEEEEE", ec="none")

berg_track = sector.add_track((70, 90))
dgrp_track = sector.add_track((50, 70)) 

# Plot forward CDS
cds_track.genomic_features(
    gff.extract_features("CDS", target_strand=1),
    plotstyle="arrow",
    r_lim=(95, 100),
    fc="salmon",
)
# Plot reverse CDS
cds_track.genomic_features(
    gff.extract_features("CDS", target_strand=-1),
    plotstyle="arrow",
    r_lim=(90, 95),
    fc="skyblue",
)
# Extract CDS product labels
pos_list, labels = [], []
for f in gff.extract_features("CDS"):
    start, end = int(str(f.location.end)), int(str(f.location.start))
    pos = (start + end) / 2
    label = f.qualifiers.get("product", [""])[0]
    if label == "" or label.startswith("hypothetical"):
        continue
    if len(label) > 20:
        label = label[:20] + "..."
    pos_list.append(pos)
    labels.append(label)

# Plot CDS product labels on outer position
cds_track.xticks(
    pos_list,
    labels,
    label_orientation="vertical",
    show_bottom_line=True,
    label_size=6,
    line_kws=dict(ec="grey"),
)
# Plot xticks & intervals on inner position
cds_track.xticks_by_interval(
    interval=1000,
    outer=True,
    show_bottom_line=True,
    label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
    label_orientation="vertical",
    line_kws=dict(ec="grey"),
)

# Load Variant Data

df1 = pd.read_csv('Bergmancount.csv')
x1 = df1['position']
y1 = df1['count'] 

df2 = pd.read_csv('DGRPcount.csv')
x2 = df2['position']
y2 = df2['count'] 

# Plot Variant Data
berg_track.bar(x1, y1, width = 50)
dgrp_track.bar(x2, y2, width = 50 )

circos.savefig("example01.png")