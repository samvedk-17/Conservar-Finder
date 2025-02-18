{% extends 'base.html' %}
{% block title %}Conservar Finder - About{% endblock %}

{% block body %}
    <div class="container my-4" style="padding-left: 3px;">
        <h3 class="text-left">How Conservar Finder Works?</h3>
    </div>
    <div class="container my-4" style="padding-left: 3px;">
        <p class="text-align">
            Conservar Finder analyzes aligned FASTA files to identify conserved and variable regions in protein or nucleotide sequences. By examining each column of the alignment relative to the <b>Reference ID</b>, the tool classifies positions as <b>Conserved</b> (if all amino acids or nucleotides in a column are identical) or <b>Variable</b> (if differences exist). Once variable positions are identified, GenBank IDs associated with these variations are retrieved from the uploaded FASTA file. For example, if Tryptophan (W) in the reference sequence is replaced by Valine (V) at position 39, Conservar Finder finds the GenBank IDs containing this variation and fetches real-time metadata from <a href="https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html">Entrez</a>, including <b>Date of Collection, Location, Strain, Isolate, and Genotype</b>. The tool also calculates the frequency and percentage occurrence of each mutation. Results are displayed in a clear tabular format, showing conserved and variable positions alongside respective amino acids or nucleotides, frequency percentages, GenBank IDs, and additional metadata, all of which can be downloaded in .csv format. To enhance visualization, Conservar Finder also maps the geographic origins of sequence variations on a global map and a GenBank ID-wise variation distribution plot can also be generated.
        </p>
    </div>
    <div class="container my-4">
        <div class="row">
            <!-- Column for the image -->
            <div class="d-flex justify-content-center align-items-center">
                <img src="{{ url_for('static', filename='Input-hz.png') }}" alt="Conservar Finder Workflow" class="img-fluid" style="width: 800px; height: auto;">
            </div>
        </div>
    </div>
    <div class="container my-4" style="padding-left: 3px;">
        <h3 class="text-left">How to use Conservar Finder?</h3>
    </div>
    <div class="container my-4" style="padding-left: 3px;">
        <ol class="text-align">
            <li>Upload your aligned FASTA file. You can use tools such as <a href="https://www.ebi.ac.uk/jdispatcher/msa/clustalo">Clustal Omega</a>, <a href="https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein">MUSCLE</a>, <a href="https://mafft.cbrc.jp/alignment/server/index.html">MAFFTT</a> to perform <a href="https://en.wikipedia.org/wiki/Multiple_sequence_alignment">Multiple Sequence Alignment</a> on your FASTA sequences. For <a href="https://en.wikipedia.org/wiki/Sequence_alignment#Pairwise_alignment">Pairwise Alignment</a> you can use <a href="https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle">EMBOSS NEEDLE</a>. You can also use <a href="https://www.jalview.org/">Jalview</a> to edit, visualise and analyse aligned sequences and to save them into .fasta or .fa format. Make sure that GenBank IDs are present in your input file.</li>
            <li>Mark the checkbox if you want to omit sequences with missing data (X) or gaps (-) from your input file.</li>
            <li>Provide your E-mail ID which necessary to fetch metadata from Entrez.</li>
            <li>Click on "Upload" button.</li>
            <li>You can also download results in .csv format.</li>
        </ol>
    </div>

    <div class="container my-4" style="padding-left: 3px;">
        <h3 class="text-left">Results</h3>
        <div class="results-section my-4">
            <div class="result-item mb-5">
                <h4>1. Conserved positions table (Reference GenBank ID: AF266290.1):</h4>
                <div class="text-center">
                    <img src="{{ url_for('static', filename='consv_table.png') }}" alt="conserved positions table" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                </div>
            </div>
            
            <div class="result-item mb-5">
                <h4>2. Variable positions table:</h4>
                <div class="text-center">
                    <img src="{{ url_for('static', filename='variable_pos_table.png') }}" alt="variable positions table" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                    <img src="{{ url_for('static', filename='csv_down.png') }}" alt="csv download" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                </div>
            </div>
            
            <div class="result-item mb-5">
                <h4>3. Metadata table:</h4>
                <div class="text-center">
                    <img src="{{ url_for('static', filename='metadata_table.png') }}" alt="metadata table" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                </div>
            </div>
            
            <div class="result-item mb-5">
                <h4>4. Geographic Map of Variable Positions:</h4>
                <div class="text-center">
                    <img src="{{ url_for('static', filename='maps.png') }}" alt="geographic Map of Variable Positions" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                </div>
            </div>
            
            <div class="result-item mb-5">
                <h4>5. GenBank ID Variation Distribution:</h4>
                <div class="text-center">
                    <img src="{{ url_for('static', filename='genbank_var_plot.png') }}" alt="genBank ID Variation Distribution" class="img-fluid mb-3" style="width: 1500px; height: auto;">
                </div>
            </div>
        </div>
    </div>
{% endblock body %}
      

    

