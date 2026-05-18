# Data Format Specifications

This directory should contain your input data files. The pipeline expects:

## otutab.txt (TSV format)

A tab-separated OTU/ASV abundance table.

```
OTUID	S1	S2	S3	S4	...
ASV_1	12	0	45	67	...
ASV_2	0	89	0	12	...
```

- **First row**: sample IDs as column headers
- **OTUID**: unique ASV/OTU identifier
- **Cells**: integer read counts

## taxonomy.txt (TSV format)

A tab-separated taxonomy annotation file.

```
OTUID	Kingdom	Phylum	Class	Order	Family	Genus	Species
ASV_1	Bacteria	Firmicutes	Bacilli	Lactobacillales	Carnobacteriaceae	Trichococcus	Unassigned
ASV_2	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae_1	Clostridium_sensu_stricto	Unassigned
```

## Example Data

For testing, place small example files here. The repository `.gitignore` excludes actual data files — users should replace with their own data.
