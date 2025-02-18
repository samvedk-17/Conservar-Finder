**How Conservar Finder Works?**

Conservar Finder analyzes aligned FASTA files to identify conserved and variable regions in protein or nucleotide sequences. By examining each column of the alignment relative to the **Reference ID**, the tool classifies positions as **Conserved** (if all amino acids or nucleotides in a column are identical) or **Variable** (if differences exist). Once variable positions are identified, GenBank IDs associated with these variations are retrieved from the uploaded FASTA file. For example, if Tryptophan (W) in the reference sequence is replaced by Valine (V) at position 39, Conservar Finder finds the GenBank IDs containing this variation and fetches real-time metadata from [**Entrez**](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html), including **Date of Collection, Location, Strain, Isolate, and Genotype**. The tool also calculates the frequency and percentage occurrence of each mutation. Results are displayed in a clear tabular format, showing conserved and variable positions alongside respective amino acids or nucleotides, frequency percentages, GenBank IDs, and additional metadata, all of which can be downloaded in .csv format. To enhance visualization, Conservar Finder also maps the geographic origins of sequence variations on a global map and a GenBank ID-wise variation distribution plot can also be generated.

![Input-hz](https://github.com/user-attachments/assets/5a583718-defd-4c45-84f9-df7286df027f)

**How to use Conservar Finder?**
1. Upload your aligned FASTA file. You can use tools such as **Clustal Omega, MUSCLE, MAFFTT** to perform Multiple Sequence Alignment on your FASTA sequences. For Pairwise Alignment you can use **EMBOSS NEEDLE**. You can also use **Jalview** to edit, visualise and analyse aligned sequences and to save them into .fasta or .fa format. Make sure that GenBank IDs are present in your input file.
2. Mark the checkbox if you want to omit sequences with missing data (X) or gaps (-) from your input file.
3. Provide your E-mail ID which necessary to fetch metadata from Entrez.
4. Click on "Upload" button.
5. You can also download results in .csv format.

**User Interface:**
![GUI](https://github.com/user-attachments/assets/4fcfc0c1-285f-45a9-811c-106ded57a82e)

**Results**

**1. Conserved positions table (Reference GenBank ID: AF266290.1):** 
![consv_table](https://github.com/user-attachments/assets/3b767177-452e-467e-bdef-f7e9dc78d506)
**2. Variable positions table:** 
![variable_pos_table](https://github.com/user-attachments/assets/b3572d4b-3a8d-4548-a321-bb2033218136)
**3. Metadata table:**
![metadata_table](https://github.com/user-attachments/assets/378ae1cd-67a9-45a6-8603-158157311684)
**4.Geographic Map of Variable Positions:**
![maps](https://github.com/user-attachments/assets/a84f4a04-77cf-43c6-bd1e-8223efbed9a2)
**5. GenBank ID Variation Distribution:**
![genbank_var_plot](https://github.com/user-attachments/assets/31e809b5-134b-48eb-a9d9-a4b1fb33e441)



