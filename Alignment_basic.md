```markdown
# Multiple Alignment of Protein Sequences Script

## Overview
This is the script created on **26/11/2024** for multiple alignment of protein sequences fetched by accession numbers from PubMed.

### Cool about it:
- Fetches multiple sequences by accession numbers.
- Converts them into FASTA format.
- Aligns and generates pretty images.

### Lousy about it:
- Uses both R and Python. 
- Python fetches and creates a file, R performs the alignment, and Python draws. 
- This is due to difficulties in managing "command line Clustal" for Python. Needs improvement!

## Used Packages
- **Images**: [pyMSAviz Documentation](https://moshi4.github.io/pyMSAviz/getting_started/#2-customized-visualization)
- **FASTA tools for Python**: [Reading and Writing FASTA Files](https://github.com/zaneveld/full_spectrum_bioinformatics/blob/master/content/06_biological_sequences/reading_and_writing_fasta_files.ipynb)
- **MSA tools for R**: [Bioconductor MSA Documentation](https://bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf)

## Align AAV Capsid Proteins

### Python Fetching and Imaging

```python
from pymsaviz import MsaViz, get_msa_testdata
import os
from Bio import AlignIO
from Bio import Entrez
from Bio.Align.Applications import ClustalwCommandline

os.chdir(r'C:\Users\Елена\Downloads')

def fetch_protein_sequence(protein_id):
    Entrez.email = "elena.n.kozhevnikova@gmail.com"
    try:
        # Fetch the protein record from the NCBI protein database
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        sequence_data = handle.read()
        handle.close()
        
        return sequence_data
    except Exception as e:
        print(f"An error occurred while fetching {protein_id}: {e}")
        return None

# List of protein IDs to fetch
protein_ids = [
    'NP_049542.1', 'YP_680426.1', 'NP_043941.1', 'NP_044927.1',
    'YP_068409.1', 'AAB95450.1', 'YP_077178.1', 'YP_077180.1',
    'AAS99264.1'
]
# Output file to store the sequences
output_file = "protein_sequences.fasta"

# Open the output file in write mode
with open(output_file, 'w') as f:
    for protein_id in protein_ids:
        sequence = fetch_protein_sequence(protein_id)
        if sequence:
            f.write(sequence)  # Write the fetched sequence to the file

print(f"Protein sequences have been saved to {output_file}.")
```

:warning: **This part is interrupted and used further with R-made file**: Be very careful here!

```python
alignment = AlignIO.read("out.fasta", "fasta")

mv = MsaViz(alignment, wrap_length=60, color_scheme="Flower", show_grid=True, show_consensus=True)
fig = mv.plotfig()

fig.savefig('alignment_plot.pdf', bbox_inches='tight')  # Change the filename and format as needed
```

### R Alignment

```r
library(msa)
library(seqinr)

setwd("C:/Users/Елена/Downloads")

alignment2Fasta <- function(alignment, filename) {
    sink(filename)
  
    n <- length(rownames(alignment))
    for(i in seq(1, n)) {
        cat(paste0('>', rownames(alignment)[i]))
        cat('\n')
        the.sequence <- toString(unmasked(alignment)[[i]])
        cat(the.sequence)
        cat('\n')  
    }
  
    sink(NULL)
}

mySeqs <- readAAStringSet("protein_sequences.fasta")
myAlignment <- msa(mySeqs)
alignment2Fasta(myAlignment, 'out.fasta')
```
```

### Instructions for Use
1. Ensure all required packages are installed in both Python and R environments.
2. Modify the script as necessary for your specific requirements.
3. Update file paths according to your system structure.
4. Run the Python script first to fetch the sequences, then run the R script for alignment.

### Output Plot

The alignment plot has been generated and saved as `alignment_plot.pdf`. You can view or download it using the link below:

[View the alignment plot](alignment_plot.pdf)


