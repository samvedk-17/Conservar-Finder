from flask import Flask, render_template, request, redirect, url_for, send_file , session
import os
import pandas as pd
from Bio import AlignIO, Entrez, SeqIO
import folium
from geopy.geocoders import Nominatim
import geopandas as gpd
from folium.plugins import MarkerCluster
from werkzeug.utils import secure_filename
from concurrent.futures import ThreadPoolExecutor
import io
from datetime import datetime

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
app.secret_key = os.urandom(24) 
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


Entrez.email = ''  # User's email for Entrez API

# Optimized function to find conserved and variable positions with parallelized filtering
def find_conserved_and_variable_pos(input_fasta, omit_missing):
    alignment = AlignIO.read(input_fasta, 'fasta')
    
    with ThreadPoolExecutor() as executor:
       filtered_alignment = list(executor.map(lambda record: record if not (omit_missing and ('-' in record.seq or 'X' in record.seq)) else None, alignment))
       filtered_alignment = [record for record in filtered_alignment if record]
      

    conserved = []
    variable = {}
    alignment_length = len(filtered_alignment[0].seq) if filtered_alignment else 0

    for i in range(alignment_length):
        column = [record.seq[i] for record in filtered_alignment]
        unique_aa = set(column)
        position = i + 1  # 1-based index

        if len(unique_aa) == 1:
            conserved.append(f"{position}: {column[0]}")
        else:
            variable[position] = {'amino_acids': unique_aa, 'ids': [record.id for record in filtered_alignment if record.seq[i] in unique_aa]}

    return conserved, variable

# Optimized function for extracting GenBank IDs from sequence headers
def extract_fasta_ids(fasta_ids):
    return [entry.split('|')[1].split('_prot_')[0] for entry in fasta_ids]

# Optimized batch metadata fetching function
def fetch_genbank_metadata(genbank_ids, email, variable_positions):
    Entrez.email = email
    batch_size = 50
    records_data = []
    
    def fetch_batch(ids_batch):
        handle = Entrez.efetch(db="nucleotide", id=','.join(ids_batch), rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "genbank")
        data_batch = []
        
        for record in records:
            for feature in record.features:
                if feature.type == "source":
                    location = feature.qualifiers.get("geo_loc_name", ["N/A"])[0]
                    collection_date = feature.qualifiers.get("collection_date", ["N/A"])[0]
                    strain = feature.qualifiers.get("strain", ["N/A"])[0]
                    isolate = feature.qualifiers.get("isolate", ["N/A"])[0]
                    genotype = feature.qualifiers.get("note", ["N/A"])[0]

                    for pos, amino_acid_set in variable_positions.items():
                        for amino_acid in amino_acid_set['amino_acids']:
                            data_batch.append({
                                "Position": pos,
                                "Amino Acids": amino_acid,
                                "GenBank ID": record.id,
                                "Collection Date": collection_date,
                                "Location": location,
                                "Strain": strain,
                                "Isolate": isolate,
                                "Genotype": genotype
                            })
        handle.close()
        return data_batch

    with ThreadPoolExecutor() as executor:
        for i in range(0, len(genbank_ids), batch_size):
            records_data.extend(executor.submit(fetch_batch, genbank_ids[i:i + batch_size]).result())

    return pd.DataFrame(records_data).sort_values(by=['Position', 'Amino Acids'])

# Updated function for generating map with optimized coordinate handling
def generate_map(df):
    geolocator = Nominatim(user_agent="geoapi")

    def get_coordinates(location):
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

    df['Coordinates'] = df['Location'].apply(get_coordinates)
    df = df[df['Coordinates'].notnull()]
    df_unique = df.groupby(['GenBank ID', 'Location']).first().reset_index()

    gdf = gpd.GeoDataFrame(df_unique, geometry=gpd.points_from_xy(
        df_unique['Coordinates'].apply(lambda x: x[1]),  # Longitude
        df_unique['Coordinates'].apply(lambda x: x[0])   # Latitude
    ))

    folium_map = folium.Map(location=[20, 0], zoom_start=2, tiles='CartoDB Positron')
    marker_cluster = MarkerCluster().add_to(folium_map)

    for _, row in gdf.iterrows():
        popup_text = f"""
        <b>GenBank ID:</b> {row['GenBank ID']}<br>
        <b>Strain:</b> {row['Strain']}<br>
        <b>Location:</b> {row['Location']}<br>
        """
        tooltip_text = f"""
        <b>GenBank ID:</b> {row['GenBank ID']}<br>
        <b>Strain:</b> {row['Strain']}<br>
        <b>Location:</b> {row['Location']}<br>
        """
        folium.Marker(
            location=[row['Coordinates'][0], row['Coordinates'][1]],
            popup=popup_text,
            tooltip=folium.Tooltip(tooltip_text)
        ).add_to(marker_cluster)

    folium_map.save("static/map.html")

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
@app.route('/results', methods=['POST'])
def upload_file():
    global conserved_pos , variable_pos_data , metadata_df
    conserved_pos = [] 
    variable_pos_data = {}
    metadata_df = {}
    uploaded_file = request.files['file']
    email = request.form['email']
    omit_missing = 'omit_missing' in request.form

    if uploaded_file.filename != '' and email != '':
        filename = secure_filename(uploaded_file.filename)
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        uploaded_file.save(file_path)

        conserved_pos, variable_pos_data = find_conserved_and_variable_pos(file_path, omit_missing)
        session['uploaded_filename'] = filename
        base_filename = os.path.splitext(filename)[0] 

        variable_ids = set()
        for position_data in variable_pos_data.values():
            fasta_ids = position_data['ids']
            extracted_ids = extract_fasta_ids(fasta_ids)
            variable_ids.update(extracted_ids)

        genbank_ids = list(variable_ids)
        metadata_df = fetch_genbank_metadata(genbank_ids, email, variable_pos_data) if genbank_ids else pd.DataFrame()

        if not metadata_df.empty:
            generate_map(metadata_df)


        date = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")  # Get current UTC date and time
        return render_template('results.html', conserved=conserved_pos, variable=variable_pos_data, metadata=metadata_df.to_dict(orient="records") , date=date)

    return redirect(url_for('index'))

@app.route('/About')
def about():
    return render_template('about.html')

@app.route('/Contact')
def contact():
    return render_template('contact.html')

###Download Data##
@app.route('/download_consv')
def download_consv():
    global conserved_pos
    if not conserved_pos:  # Initialize with empty list if undefined
        conserved_pos = []
    uploaded_filename = session.get('uploaded_filename', 'results')  # Get uploaded filename or default name    

    csv_data = "Conserved Position, Conserved Amino Acid\n"
    csv_data += "\n".join([f"{pos.split(': ')[0]},{pos.split(': ')[1]}" for pos in conserved_pos])
    
    output = io.BytesIO()
    output.write(csv_data.encode('utf-8'))
    output.seek(0)
    
    # Use the uploaded filename to create the CSV name
    csv_filename = f"{os.path.splitext(uploaded_filename)[0]}_conserved_positions.csv"
    
    return send_file(output, as_attachment=True, download_name= csv_filename , mimetype='text/csv')

@app.route('/download_var')
def download_var():
    global variable_pos_data
    if not variable_pos_data:  # Initialize with empty list if undefined
        variable_pos_data = {}
    uploaded_filename = session.get('uploaded_filename', 'results')  # Get uploaded filename or default name    

    csv_data = "Variable Position, Variable Amino Acid(s)\n"
    csv_data += "\n".join([f"{position},\"{', '.join(data['amino_acids'])}\"" for position, data in variable_pos_data.items()])
    output = io.BytesIO()
    output.write(csv_data.encode('utf-8'))
    output.seek(0)
    
    # Use the uploaded filename to create the CSV name
    csv_filename = f"{os.path.splitext(uploaded_filename)[0]}_variable_positions.csv"
    
    return send_file(output, as_attachment=True, download_name= csv_filename , mimetype='text/csv')

@app.route('/download_metadata')
def download_metadata():
    global variable_pos_data , metadata_df
    if not variable_pos_data and metadata_df:  # Initialize with empty list if undefined
        variable_pos_data = {}
        metadata_df = {}

    uploaded_filename = session.get('uploaded_filename', 'results')  # Get uploaded filename or default name    

    csv_data = "Variable Position, Variable Amino Acid , Genbank id , Collection date , Loaction , Strain , Isolate , Genotype\n"
    csv_data += "\n".join([
        f"{data['Position']},\"{(data['Amino Acids'])}\","
        f"{data['GenBank ID']},"
        f"{data['Collection Date']},"
        f"{data['Location']},"
        f"{data['Strain']},"
        f"{data['Isolate']},"
        f"{data['Genotype']}"
        for data in metadata_df.to_dict(orient='records')  # Convert DataFrame to a list of dictionaries
    ])
   

    
    output = io.BytesIO()
    output.write(csv_data.encode('utf-8'))
    output.seek(0)
    
    # Use the uploaded filename to create the CSV name
    csv_filename = f"{os.path.splitext(uploaded_filename)[0]}_metadata.csv"
    
    return send_file(output, as_attachment=True, download_name= csv_filename , mimetype='text/csv')
        
    
if __name__ == '__main__':
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)
