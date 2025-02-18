from flask import Flask, render_template, request, redirect, url_for, send_file, render_template_string
import os
import pandas as pd
from Bio import AlignIO, Entrez, SeqIO
from Bio.Align import MultipleSeqAlignment
import folium
from geopy.geocoders import Nominatim
import geopandas as gpd
from folium.plugins import MarkerCluster
from werkzeug.utils import secure_filename
from concurrent.futures import ThreadPoolExecutor
import io
from datetime import datetime
from functools import lru_cache
from collections import Counter, defaultdict
import csv
import time
import plotly.graph_objects as go


app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
app.secret_key = os.urandom(24)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


# Global variables with proper initialization
class GlobalState:
    def __init__(self):
        self.conserved_pos = []
        self.variable_pos_data = {}
        self.metadata_df = pd.DataFrame()
        self.uploaded_filename = ''

global_state = GlobalState()

Entrez.email = ''  # User's email for Entrez API

# Optimized function to find conserved and variable positions with parallelized filtering
def find_variations_and_conserved(input_fasta, ref_seq_id, omit_gaps=False):
    """Optimized function to find variations and conserved positions with proper gap handling"""
    alignment = AlignIO.read(input_fasta, 'fasta')
    
    ref_seq = next((record.seq for record in alignment if ref_seq_id in record.id), None)
    if ref_seq is None:
        raise ValueError(f"Reference sequence with ID {ref_seq_id} not found")
    
    variations = {}
    conserved = []
    position_cache = defaultdict(Counter)
    for record in alignment:
        # Skip the reference sequence when counting
        if ref_seq_id in record.id:
            continue
        for i, aa in enumerate(record.seq):
            # Skip gaps and 'X' if omit_gaps is True
            if omit_gaps and (aa == '-' or aa == 'X'):
                continue
            position_cache[i][aa] += 1
    
    # Analyze positions
    for i in range(len(ref_seq)):
        pos = i + 1
        ref_aa = ref_seq[i]
        # Skip positions where reference has a gap if omit_gaps is True
        if omit_gaps and (ref_aa == '-' or ref_aa == 'X'):
            continue
        counts = position_cache[i]
        # Calculate total count excluding gaps if omit_gaps is True
        if omit_gaps:
            counts['-'] = 0
            counts['X'] = 0
        total_count = sum(counts.values())
        if total_count == 0:
            continue
        
        unique_aa = {aa for aa in counts.keys() if not (omit_gaps and aa in ['-', 'X'])}
        # Calculate percentages based on total count
        percentages = {aa: round((counts[aa]/total_count)*100, 1) 
                      for aa in unique_aa}
        if len(unique_aa) == 1 and list(unique_aa)[0] == ref_aa:
            conserved.append(f"{pos}: {ref_aa}")
        else:
            variations[pos] = {
                'ref_aa': ref_aa,
                'amino_acids': {aa: count for aa, count in counts.items() 
                              if not (omit_gaps and aa in ['-', 'X'])},
                'percentages': percentages,
                'mutations': [f"{ref_aa}{pos}{aa}" for aa in unique_aa 
                            if aa != ref_aa and not (omit_gaps and aa in ['-', 'X'])],
                'ids': defaultdict(list)
            }
            # Collect IDs for variations, respecting gap handling
            for record in alignment:
                # Skip reference sequence when collecting IDs
                if ref_seq_id in record.id:
                    continue
                aa = record.seq[i]
                if aa != ref_aa and not (omit_gaps and aa in ['-', 'X']):
                    variations[pos]['ids'][aa].append(record.id)
    
    return conserved, variations, alignment

def fetch_genbank_metadata(genbank_ids, email, variable_positions, batch_size=100):
    """Optimized metadata fetching with thread pooling"""
    if not email:
        print("Error: No email provided for Entrez API")
        return pd.DataFrame()
        
    Entrez.email = email
    records_data = []
    print(f"Fetching metadata for {len(genbank_ids)} GenBank IDs")
    
    def process_batch(ids_batch):
        try:
            with Entrez.efetch(db="nucleotide", id=','.join(ids_batch), 
                             rettype="gb", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "genbank"))
                print(f"Fetched {len(records)} records for batch")
                batch_data = []
                
                for record in records:
                    features = extract_features(record)
                    base_id = record.id.split('_prot_')[0].split('|')[-1]
                    
                    for pos, var_data in variable_positions.items():
                        if base_id in str(var_data['ids']):
                            batch_data.append(create_record_data(pos, var_data, base_id, features))
                            
                return batch_data
        except Exception as e:
            app.logger.error(f"Error in batch processing: {e}")
            if "HTTP Error 503" in str(e) or "HTTP Error 429" in str(e):
                print("\n⚠️  NCBI SERVER IS CURRENTLY EXPERIENCING ISSUES OR IS DOWN")
                print("Error details:", str(e))
                print("Please try again later.\n")
            
            elif "HTTP Error 400" in str(e):
                print("\n⚠️ INVALID REQUEST TO NCBI SERVER")
                print("Error details:", str(e))
                print("Please check your input data and try again.\n")

            else:
                print(f"\n⚠️ Error fetching metadata from NCBI: {e}\n")    
            return []
    
    # Process batches using thread pool
    batches = [genbank_ids[i:i + batch_size] for i in range(0, len(genbank_ids), batch_size)]
    with ThreadPoolExecutor(max_workers=min(len(batches), 10)) as executor:
        batch_results = list(executor.map(process_batch, batches))
        
    # Merge results
    for batch in batch_results:
        records_data.extend(batch)

    df = pd.DataFrame(records_data).sort_values(by=['Position', 'Variant AA'])
    print("DataFrame First Few Rows:\n", df.head())

    return df

def extract_features(record):
    """Helper function to extract features from GenBank record"""
    features = {
        'Location': "N/A",
        'Collection Date': "N/A",
        'Strain': "N/A",
        'Isolate': "N/A",
        'Genotype': "N/A"
    }
    
    for feature in record.features:
        if feature.type == "source":
            base_id = record.id.split('_prot_')[0].split('|')[-1]
            qualifiers = feature.qualifiers
            features.update({
                'Location': qualifiers.get("geo_loc_name", ["N/A"])[0],
                'Collection Date': qualifiers.get("Collection Date", ["N/A"])[0],
                'Strain': qualifiers.get("Strain", ["N/A"])[0],
                'Isolate': qualifiers.get("Isolate", ["N/A"])[0],
                'Genotype': qualifiers.get("note", ["N/A"])[0]
            })
            break
    print(f"Processing record: {base_id}")
    
    return features

def create_record_data(pos, var_data, base_id, features):
    """Helper function to create record data dictionary"""
    mutation = "N/A"
    variant_aa = "N/A"
    
    # Find the specific mutation for the current GenBank ID
    for aa, ids in var_data['ids'].items():
        if base_id in str(ids):
            mutation = f"{var_data['ref_aa']}{pos}{aa}"
            variant_aa = aa
            break

    return {
        "Position": pos,
        "Reference AA": var_data['ref_aa'],
        "Variant AA": variant_aa,
        "Mutation": mutation,
        "GenBank ID": base_id,
        **features
    }

# Updated function for generating map with optimized coordinate handling
def generate_map(df):
    geolocator = Nominatim(user_agent="geoapi")
    
    @lru_cache(maxsize=1024)
    def get_coordinates(location):
        location = str(location)
        try:
            loc = geolocator.geocode(location)
            if loc:
                return loc.latitude, loc.longitude
            else:
                simplified_location = location.split(":")[0]
                loc = geolocator.geocode(simplified_location)
                return (loc.latitude, loc.longitude) if loc else None
        except Exception:
            return None

    # Make sure we're using the correct column name from the metadata
    # Filter out rows where location is missing or N/A
    df = df[df['Location'] != 'N/A']  
    df['Coordinates'] = df['Location'].apply(get_coordinates)  
    df = df[df['Coordinates'].notnull()]
    df_unique = df.groupby(['GenBank ID', 'Location']).first().reset_index()  

    gdf = gpd.GeoDataFrame(df_unique, geometry=gpd.points_from_xy(
        df_unique['Coordinates'].apply(lambda x: x[1]),
        df_unique['Coordinates'].apply(lambda x: x[0])
    ))

    folium_map = folium.Map(location=[20, 0], zoom_start=2, tiles='CartoDB Positron')
    marker_cluster = MarkerCluster().add_to(folium_map)

    for _, row in gdf.iterrows():
        popup_text = f"""
        <b>GenBank ID:</b> {row['GenBank ID']}<br>
        <b>Genotype:</b> {row['Genotype']}<br>  
        <b>Strain:</b> {row['Strain']}<br>
        <b>Location:</b> {row['Location']}<br>
        """
        folium.Marker(
            location=[row['Coordinates'][0], row['Coordinates'][1]],
            popup=popup_text,
            tooltip=f"GenBank ID: {row['GenBank ID']}"
        ).add_to(marker_cluster)
        
    folium_map.save("static/map.html")

#Helper function to get GenBank IDs.
def get_base_genbank_id(full_id):
    """Extract base GenBank ID from the full identifier."""
    # Split on '_prot_' and take the first part
    base_id = full_id.split('_prot_')[0]
    # If there are pipe characters, take the last segment
    if '|' in base_id:
        base_id = base_id.split('|')[-1]
    return base_id

#Data formating
def format_variable_positions(variable_pos_data):
    """Format variable positions data for template rendering"""
    if not isinstance(variable_pos_data, dict):
        return {}

    formatted_data = {}
    for position, data in variable_pos_data.items():
        variations = []
        genbank_ids = {}
        frequencies = {}

        total_sequences = sum(data['amino_acids'].values())

        for aa, count in data['amino_acids'].items():
            if aa != data['ref_aa']:
                percentage = data['percentages'].get(aa, 0)
                variations.append(f"{aa}({percentage:.1f}%)")
                frequencies[aa] = count
                base_ids = [get_base_genbank_id(id_str) for id_str in data['ids'].get(aa, [])]
                genbank_ids[aa] = ', '.join(base_ids)
                
        if variations:
            formatted_data[position] = {
                'ref_aa': data['ref_aa'],
                'variations': ', '.join(variations),
                'mutations': data.get('mutations', []),
                'frequencies': frequencies,
                'genbank_ids': genbank_ids
            }

    return formatted_data
        
                               
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
@app.route('/results', methods=['POST'])
def upload_file():
    start_time = time.time()
    step_times = []

    try:
        uploaded_file = request.files['file']
        email = request.form['email']
        ref_seq_id = request.form['ref_seq_id']
        omit_gaps = bool(request.form.get('omit_missing'))  # Properly get the checkbox value

        if uploaded_file.filename and email and ref_seq_id:
            filename = secure_filename(uploaded_file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            uploaded_file.save(file_path)
            global_state.uploaded_filename = filename

            # Process data with explicit omit_gaps parameter
            step1_start = time.time()
            global_state.conserved_pos, global_state.variable_pos_data, filtered_alignment = \
                find_variations_and_conserved(file_path, ref_seq_id, omit_gaps=omit_gaps)
            step_times.append(("Find variations and conserved positions", time.time() - step1_start))

            # Get unique GenBank IDs, excluding those from gap positions if omit_gaps is True
            genbank_ids = set()
            for pos_data in global_state.variable_pos_data.values():
                for aa, ids in pos_data['ids'].items():
                    if not (omit_gaps and aa in ['-', 'X']):
                        genbank_ids.update([id.split('_prot_')[0].split('|')[-1] for id in ids])
            print("Collected Variable IDs for metadata fetching:")
            print(genbank_ids)
       
            # Fetch metadata
            metadata_start_time = time.time()
            global_state.metadata_df = fetch_genbank_metadata(
                list(genbank_ids), email, global_state.variable_pos_data)
            step_times.append(("Fetch metadata", time.time() - metadata_start_time))

            # Format data for template
            variable_pos_data_formatted = format_variable_positions(global_state.variable_pos_data)
            metadata_data = global_state.metadata_df.to_dict('records') if not global_state.metadata_df.empty else []
            print("\nMetadata Fetched:")
            print(global_state.metadata_df)

            # Generate map
            try:
                map_start_time = time.time()
                generate_map(global_state.metadata_df)
                step_times.append(("Generate map", time.time() - map_start_time))
            except Exception as e:
                app.logger.error(f"Error generating map: {str(e)}")

            # Log total time
            total_time = time.time() - start_time
            step_times.append(("Total process time", total_time))

            # Log the timing to the console (optional)
            print("\nTiming Breakdown:")
            for step, duration in step_times:
                print(f"{step}: {duration:.2f} seconds")
        
            return render_template('results.html',
                                conserved=global_state.conserved_pos,
                                variable=variable_pos_data_formatted,
                                metadata=metadata_data,
                                date=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
                                timing=step_times,
                                omit_gaps=omit_gaps)  # Pass omit_gaps to template

    except Exception as e:
        import traceback
        print("Debug - Full traceback:")
        print(traceback.format_exc())
        app.logger.error(f"Error in processing file upload: {str(e)}")
        return f"<h3>Error occurred during file processing: {str(e)}</h3>", 500

    return redirect(url_for('index'))

@app.route('/About')
def about():
    return render_template('about.html')

@app.route('/Contact')
def contact():
    return render_template('contact.html')

@app.route('/download_consv')
def download_consv():
    if not global_state.conserved_pos:
        return "<h3>No conserved positions data available</h3>", 404
        
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(['Position', 'Conserved_AA'])
    
    for pos in global_state.conserved_pos:
        if ':' in pos:
            position, aa = pos.split(':', 1)
            writer.writerow([position.strip(), aa.strip()])
            
    output.seek(0)
    csv_filename = f"{os.path.splitext(global_state.uploaded_filename)[0]}_conserved_positions.csv"
    
    return send_file(
        io.BytesIO(output.getvalue().encode('utf-8')),
        mimetype='text/csv',
        as_attachment=True,
        download_name=csv_filename
    )

@app.route('/download_var')
def download_variable():
    if not global_state.variable_pos_data:
        return "<h3>No variable positions data available</h3>", 404
        
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(['Position', 'Reference_AA', 'Variant_AA', 'Mutation', 'Count', 'Percentage', 'GenBank_IDs'])
    
    for position, data in global_state.variable_pos_data.items():
        for aa, count in data['amino_acids'].items():
            if aa != data['ref_aa']:
                percentage = data['percentages'][aa]
                mutation = f"{data['ref_aa']}{position}{aa}"
                base_ids = [get_base_genbank_id(id_str) for id_str in data['ids'].get(aa, [])]
                
                writer.writerow([
                    position,
                    data['ref_aa'],
                    aa,
                    mutation,
                    count,
                    f"{percentage:.1f}%",
                    ','.join(base_ids)
                ])
            
    output.seek(0)
    csv_filename = f"{os.path.splitext(global_state.uploaded_filename)[0]}_variable_positions.csv"
    
    return send_file(
        io.BytesIO(output.getvalue().encode('utf-8')),
        mimetype='text/csv',
        as_attachment=True,
        download_name=csv_filename
    )
    
@app.route('/full_variable_positions')
def full_variable_positions():
    if not global_state.variable_pos_data:
        return "<h3>No variable positions data available</h3>"
    
    data = []
    for position, data_dict in global_state.variable_pos_data.items():
        variations = []
        frequencies = []
        genbank_ids = []

        for aa, count in data_dict['amino_acids'].items():
            if aa != data_dict['ref_aa']:  # Only include if different from reference
                percentage = data_dict['percentages'].get(aa, 0)
                variations.append(f"{aa}({percentage:.1f}%)")
                frequencies.append(f"{aa}:{count}")
                base_ids = [get_base_genbank_id(id_str) for id_str in data_dict['ids'].get(aa, [])]
                genbank_ids.append(f"{aa}:{','.join(base_ids)}")
                
        if variations:
            row = {
                'Position': position,
                'Reference AA': data_dict['ref_aa'],
                'Variations': ', '.join(variations),
                'Mutations': ', '.join(data_dict.get('mutations', [])),
                'Frequency': ', '.join(frequencies),
                'GenBank IDs': '; '.join(genbank_ids)
            }
            data.append(row)
    
    df = pd.DataFrame(data).sort_values('Position')
    
    table_html = df.to_html(
        classes='table table-bordered table-striped table-hover table-info',
        index=False,
        escape=False
    )
    
    return render_template_string("""
        <!DOCTYPE html>
        <html>
            <head>
                <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/4.5.2/css/bootstrap.min.css">
                <title>Variable Positions - Complete Table</title>
                <style>
                    .table { margin: 20px; }
                    .table-info thead th { 
                        background-color: #a4d5f8;
                        position: sticky;
                        top: 0;
                        z-index: 1;
                    }
                    .container-fluid {
                        padding: 20px;
                    }
                    /* Add styles for better text wrapping in cells */
                    .table td {
                        white-space: normal;
                        word-wrap: break-word;
                        max-width: 300px;
                    }
                    /* Ensure the GenBank IDs column doesn't get too wide */
                    .table td:last-child {
                        max-width: 250px;
                    }
                </style>
            </head>
            <body>
                <div class="container-fluid">
                    <h2 class="mb-4">Variable Positions</h2>
                    <div class="table-responsive">
                        {{ table_html|safe }}
                    </div>
                    <a href="/" class="btn btn-primary mt-3">Back to Home</a>
                </div>
            </body>
        </html>
    """, table_html=table_html)

@app.route('/full_metadata')
def full_metadata():
    if global_state.metadata_df.empty:
        return "<h3>No Metadata Available</h3>"
    
    df = global_state.metadata_df.copy()
    df = df.sort_values(['Position', 'GenBank ID'])
    
    table_html = df.to_html(
        classes='table table-bordered table-striped table-hover table-info',
        index=False,
        escape=False
    )
    
    return render_template_string("""
        <!DOCTYPE html>
        <html>
            <head>
                <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/4.5.2/css/bootstrap.min.css">
                <title>Metadata - Complete Table</title>
                <style>
                    .table { margin: 20px; }
                    .table-info thead th { 
                        background-color: #a4d5f8;
                        position: sticky;
                        top: 0;
                        z-index: 1;
                    }
                    .container-fluid {
                        padding: 20px;
                    }
                </style>
            </head>
            <body>
                <div class="container-fluid">
                    <h2 class="mb-4">Complete Metadata</h2>
                    <div class="table-responsive">
                        {{ table_html|safe }}
                    </div>
                    <a href="/" class="btn btn-primary mt-3">Back to Home</a>
                </div>
            </body>
        </html>
    """, table_html=table_html)

# Modified GenBank cluster visualization route
@app.route('/genbank_cluster')
def genbank_cluster():
    if global_state.metadata_df.empty:
        return "<h3>No Metadata Available for Cluster Plot</h3>"
    
    df = global_state.metadata_df.copy()
    
    # Group by GenBank ID and count variations
    variation_counts = df.groupby('GenBank ID').size().reset_index(name='Variation Count')
    variation_counts = variation_counts.sort_values('Variation Count', ascending=False)
    
    # Create hover text with truncated mutation details
    hover_texts = []
    
    # Calculate color values for the colorscale
    max_variations = variation_counts['Variation Count'].max()
    min_variations = variation_counts['Variation Count'].min()
    
    for gb_id in variation_counts['GenBank ID']:
        gb_data = df[df['GenBank ID'] == gb_id]
        mutations = gb_data['Mutation'].tolist()
        positions = gb_data['Position'].tolist()
        count = len(mutations)
        
        # Limit the number of mutations shown in tooltip
        MAX_MUTATIONS_DISPLAY = 30
        if count <= MAX_MUTATIONS_DISPLAY:
            mutation_text = "<br>".join(f"Position {pos}: {mut}" 
                                      for pos, mut in zip(positions, mutations))
        else:
            mutation_text = "<br>".join(f"Position {pos}: {mut}" 
                                      for pos, mut in zip(positions[:MAX_MUTATIONS_DISPLAY], 
                                                        mutations[:MAX_MUTATIONS_DISPLAY]))
            mutation_text += f"<br>... and {count - MAX_MUTATIONS_DISPLAY} more mutations"
        
        hover_text = (
            f"<b>GenBank ID:</b> {gb_id}<br>"
            f"<b>Total Variations:</b> {count}<br>"
            "<b>Mutations:</b><br>" + mutation_text
        )
        hover_texts.append(hover_text)
    
    # Create the plot
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=variation_counts['Variation Count'],
        y=variation_counts['GenBank ID'],
        mode='markers',
        marker=dict(
            size=12,
            color=variation_counts['Variation Count'],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title='Variation Count')
        ),
        text=hover_texts,
        hoverinfo='text'
    ))
    
    fig.update_layout(
        title='GenBank ID Variation Distribution',
        xaxis_title='Number of Variations',
        yaxis_title='GenBank ID',
        height=max(600, len(variation_counts) * 30),
        margin=dict(l=150, r=50, t=50, b=50),
        hovermode='closest',
        yaxis={'autorange': 'reversed'}  
    )
    
    return render_template_string("""
        <!DOCTYPE html>
        <html>
            <head>
                <title>GenBank ID Cluster Plot</title>
                <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/4.5.2/css/bootstrap.min.css">
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            </head>
            <body>
                <div class="container-fluid">
                    <div class="row justify-content-center">
                        <div class="col-12">
                            <h2 class="text-center my-4">GenBank ID Variation Distribution</h2>
                            <div class="card">
                                <div class="card-body">
                                    {{ plot_html|safe }}
                                </div>
                            </div>
                            <div class="text-center mt-3">
                                <p class="text-muted">
                                    Distribution of variations across GenBank IDs. 
                                    Hover over points to see detailed mutation information.
                                </p>
                                <a href="/" class="btn btn-primary">Back to Home</a>
                            </div>
                        </div>
                    </div>
                </div>
            </body>
        </html>
    """, plot_html=fig.to_html(full_html=False))
if __name__ == '__main__':
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=False)
