# Conservar-Finder
What is Conservar Finder?

Conservar Finder is a tool developed to find and analyze conserved and variable positions within aligned FASTA sequences.

How Conservar Finder Works?

Conservar Finder analyzes aligned FASTA files to identify conserved and variable regions in protein or nucleotide sequences by examining each column of the alignment. It then classifies positions as either Conserved (if all amino acids or nucleotides in a column are identical) or Variable (if there are differences). Once variable positions are identified, GenBank IDs associated with the variations are retrieved from the uploaded FASTA file. For example, if Tryptophan (W) is replaced by Valine (V) at position 39, it finds the GenBank IDs that have Valine (V) at this position and fetches real-time metadata from Entrez, including Date of Collection, Location, Strain, Isolate, and Genotype. This data is displayed in a clear, tabular format showing conserved and variable positions alongside respective amino acids and additional metadata. The geographic origins of sequence variations are also marked on a global map, providing a visual overview of variations.

How to use Conservar Finder?
1. Upload your aligned FASTA file. You can use tools such as Clustal Omega, MUSCLE , MAFFTT to perform
Multiple Sequence Alignment on your FASTA sequences. For Pairwise Alignment you can use EMBOSS NEEDLE . You can also use Jalview to edit , visualise and analyse aligned sequences and to save them into .fasta or .fa format. Make sure that GenBank IDs are present in your input file.
2. Mark the checkbox if you want to omit sequences with missing data (X) or gaps (-) from your input file.
3. Provide your E-mail ID which necessary to fetch metadata from Entrez.
4. Click on "Upload" button.
5. You can also download results in .csv format.

Further areas of development: This tool takes a lot of time which increases with the size of dataset to conduct the abovementioned processes. To reduce this, time complexity of the code can be reduced or more flexibility can be provided to user to choose what processes he wants to perform on aligned FASTA files. For e.g. If a user wants only conserved and variable positions then he can mark respective checkboxes. Choice of functions to be performed can be helpful to deal with time complexity issue.
