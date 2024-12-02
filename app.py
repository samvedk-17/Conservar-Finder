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
from functools import lru_cache
import csv

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
app.secret_key = os.urandom(24) 
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


Entrez.email = ''  # User's email for Entrez API

# Optimized function to find conserved and variable positions with parallelized filtering
@lru_cache(maxsize=32)
def find_conserved_and_variable_pos(input_fasta, omit_missing):
    input_fasta = str(input_fasta)
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

        total_count = len(column)
        frequencies = {aa:column.count(aa) for aa in unique_aa}
        percentages = {aa:round((count/total_count)*100 , 1) for aa , count in frequencies.items()}

        if len(unique_aa) == 1:
            conserved.append(f"{position}: {column[0]}")
        else:
            variable[position] = {
                'amino_acids': {aa: freq for aa, freq in percentages.items() if freq <= 50.0},  # Filter
                'percentages': percentages,
                'ids': {aa: [
                    record.id for record in filtered_alignment if record.seq[i] == aa and percentages[aa] <= 50.0
                ] for aa in unique_aa if percentages[aa] <= 50.0}
            }
            

    return conserved, variable , filtered_alignment


# Optimized batch metadata fetching function
def fetch_genbank_metadata(genbank_ids, email, variable_positions):
    # Validate email
    if not email:
        print("Error: No email provided for Entrez API")
        return pd.DataFrame()

    Entrez.email = email
    batch_size = 50
    records_data = []

    print(f"Fetching metadata for {len(genbank_ids)} GenBank IDs")

    @lru_cache(maxsize=512)
    def fetch_batch(ids_batch):
        ids_batch = tuple(ids_batch)
        try:
            handle = Entrez.efetch(db="nucleotide", id=','.join(ids_batch), rettype="gb", retmode="text")
            records = list(SeqIO.parse(handle, "genbank"))
            print(f"Fetched {len(records)} records for batch")
            
            data_batch = []
            
            for record in records:
                # Extract the base GenBank ID without additional qualifiers
                base_id = record.id.split('_prot_')[0].split('|')[-1]
                
                location = "N/A"
                collection_date = "N/A"
                strain = "N/A"
                isolate = "N/A"
                genotype = "N/A"

                for feature in record.features:
                    if feature.type == "source":
                        location = feature.qualifiers.get("geo_loc_name", ["N/A"])[0]
                        collection_date = feature.qualifiers.get("collection_date", ["N/A"])[0]
                        strain = feature.qualifiers.get("strain", ["N/A"])[0]
                        isolate = feature.qualifiers.get("isolate", ["N/A"])[0]
                        genotype = feature.qualifiers.get("note", ["N/A"])[0]

                # Debug print
                print(f"Processing record: {base_id}")

                # Find positions and amino acids for this specific record
                for pos, aa_data in variable_positions.items():
                    for aa, percentage in aa_data['percentages'].items():
                        # Changed condition: include amino acids with percentages less than or equal to 50%
                        if percentage <= 50.0:
                            # Check if this record's ID matches any in the IDs for this amino acid
                            matching_records = [r_id for r_id in aa_data['ids'].get(aa, []) if base_id in r_id]
                            
                            if matching_records:
                                data_batch.append({
                                    "Position": pos,
                                    "Amino Acid": aa,
                                    "Percentage": percentage,
                                    "GenBank ID": base_id,
                                    "Collection Date": collection_date,
                                    "Location": location,
                                    "Strain": strain,
                                    "Isolate": isolate,
                                    "Genotype": genotype
                                })

            handle.close()
            return data_batch
        except Exception as e:
            print(f"Error fetching metadata: {e}")
            return []

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(fetch_batch, tuple(genbank_ids[i:i + batch_size])) for i in range(0, len(genbank_ids), batch_size)]
        
        for future in futures:
            try:
                batch_results = future.result()
                records_data.extend(batch_results)
            except Exception as e:
                print(f"Error processing future: {e}")

    print(f"Total metadata records fetched: {len(records_data)}")
    
    # Convert to DataFrame and add more detailed print statements
    df = pd.DataFrame(records_data).sort_values(by=['Position', 'Amino Acid']) if records_data else pd.DataFrame()
    
    #print("DataFrame Column Names:", list(df.columns))
    print("DataFrame First Few Rows:\n", df.head())
    
    return df

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
    metadata_df = pd.DataFrame()

    uploaded_file = request.files['file']
    email = request.form['email']
    omit_missing = 'omit_missing' in request.form

    if uploaded_file.filename != '' and email != '':
        filename = secure_filename(uploaded_file.filename)
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        uploaded_file.save(file_path)

        conserved_pos, variable_pos_data , filtered_alignment = find_conserved_and_variable_pos(file_path, omit_missing)
        session['uploaded_filename'] = filename
        base_filename = os.path.splitext(filename)[0] 
        # Debug: Print variable positions data
        print("Variable Positions Data:")
        for pos, data in variable_pos_data.items():
            print(f"Position {pos}:")
            print(f"  Amino Acids: {data['amino_acids']}")
            print(f"  Percentages: {data['percentages']}")
            #print(f"  IDs: {data['ids']}")

        variable_ids = set()
        for position, data in variable_pos_data.items():
            #print(f"\nProcessing Position {position}:")
            #print(f"Amino Acids: {data['amino_acids']}")
            #print(f"Percentages: {data['percentages']}")
         
            for aa, percentage in data['percentages'].items():
                if percentage <= 50.0:
                    ids_for_aa = [
                        record.id for record in filtered_alignment 
                        if record.seq[position - 1] == aa
                    ]
                    print(f"  Amino Acid {aa} (Percentage: {percentage:.1f}%):")
                    print(f"  IDs: {ids_for_aa}")
                    variable_ids.update(ids_for_aa)
                 
        print("Collected Variable IDs for metadata fetching:")
        print(variable_ids)

        # Fetch metadata
        
        metadata_df = pd.DataFrame()
        if variable_ids:
            clean_ids = []
            for id_str in variable_ids:
                if '|' in id_str:
                    clean_id = id_str.split('_prot_')[0].split('|')[-1]
                elif '_prot_' in id_str:
                    clean_id = id_str.split('_prot_')[0]
                else:
                    clean_id = id_str
                clean_ids.append(clean_id)

            # Remove duplicates while preserving order
            clean_ids = list(dict.fromkeys(clean_ids))
            metadata_df = fetch_genbank_metadata(clean_ids, email, variable_pos_data)
            metadata_data = metadata_df.to_dict(orient="records") if not metadata_df.empty else []
            
            print("\nMetadata Fetched:")
            print(metadata_df)
        else:
            metadata_data = []


        variable_pos_data_formatted = {}
        
        for position, data in variable_pos_data.items():
            total_count = len(data['ids'])  # Total number of sequences at this position
            formatted_aa = []
            
            for aa in data['amino_acids']:
                count = sum(
                    1 for record in filtered_alignment 
                    if record.id in data['ids'] and record.seq[position - 1] == aa
                )
                percentage = (count / total_count) * 100
                formatted_aa.append(f"{aa} ({percentage:.1f}%)")

            variable_pos_data_formatted[position] = ", ".join(formatted_aa)
            
        if not metadata_df.empty:
            generate_map(metadata_df)

        variable_pos_data_formatted = dict(list(variable_pos_data.items())[:10])
        date = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")  # Get current UTC date and time
        return render_template('results.html', conserved=conserved_pos, variable= variable_pos_data_formatted, metadata=metadata_data , date=date)

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
    csv_data += "\n".join([f"{position},\"{', '.join([f'{aa}({data['percentages'][aa]:.1f}%)' for aa in data['percentages']])}\""
                           for position, data in variable_pos_data.items()])
   
    output = io.BytesIO()
    output.write(csv_data.encode('utf-8'))
    output.seek(0)
    
    # Use the uploaded filename to create the CSV name
    csv_filename = f"{os.path.splitext(uploaded_filename)[0]}_variable_positions.csv"
    
    return send_file(output, as_attachment=True, download_name= csv_filename , mimetype='text/csv')

@app.route('/download_metadata')
def download_metadata():
    global metadata_df
    if metadata_df.empty:
        return "<h3>No Metadata Available</h3>"

    output = io.StringIO()
    csv_writer = csv.writer(output)

    # Write header
    csv_writer.writerow([
        "Position", "Variable Amino Acid", "GenBank ID", 
        "Collection Date", "Location", "Strain", "Isolate", "Genotype"
    ])

    # Write data rows
    for _, row in metadata_df.iterrows():
        csv_writer.writerow([
            row.get("Position", ""), 
            row.get("Amino Acid", ""), 
            row.get("GenBank ID", ""), 
            row.get("Collection Date", ""), 
            row.get("Location", ""), 
            row.get("Strain", ""), 
            row.get("Isolate", ""), 
            row.get("Genotype", "")
        ])

    # Prepare the output for download
    output.seek(0)
    uploaded_filename = session.get('uploaded_filename', 'results')
    csv_filename = f"{os.path.splitext(uploaded_filename)[0]}_metadata.csv"

    return send_file(
        io.BytesIO(output.getvalue().encode('utf-8')),
        as_attachment=True,
        download_name=csv_filename,
        mimetype='text/csv'
    )


@app.route('/full_metadata')
def full_metadata():
    global metadata_df  # Assuming metadata_df is a global variable holding the data
    if metadata_df.empty:
        return "<h3>No Metadata Available</h3>"

    # Drop the 'Coordinates' column from the DataFrame
    filtered_df = metadata_df.drop(columns=['Coordinates'], errors='ignore')

    # Convert the filtered DataFrame to HTML with Bootstrap table classes for styling
    table_html = filtered_df.to_html(classes='table table-bordered table-striped', index=False)
    return f"""
    <html>
        <head>
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/4.5.2/css/bootstrap.min.css">
            <title>Metadata-Complete Table</title>
        </head>
        <body>
            <div class="container mt-4">
                <h2>Metadata on variable positions</h2>
                {table_html}
                <a href="/" class="btn btn-primary mt-3">Back to Home</a>
            </div>
        </body>
    </html>
    """

@app.route('/full_variable_positions')
def full_variable_positions():
    global variable_pos_data

    if 'variable_pos_data' not in globals() or not variable_pos_data:
        return "<h3>No Variable Positions Data Available</h3>"

    data = []
    for position, data_dict in variable_pos_data.items():
        total = len(data_dict['ids'])
        amino_acids_with_percentages = ", ".join(
            f"{aa}({(data_dict['percentages'][aa]):.1f}%)"
            for aa in data_dict['percentages']
        )

        data.append([position, amino_acids_with_percentages])


    df = pd.DataFrame(data, columns=['Variable Position', 'Variable Amino Acid(s) with Percentages'])

    # Render as HTML table
    table_html = df.to_html(
        classes='table table-striped table-bordered',
        index=False  # Hide index
    )

    # Return the formatted HTML page
    return f"""
    <!DOCTYPE html>
    <html>
        <head>
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/4.5.2/css/bootstrap.min.css">
            <title>Variable Positions- Complete Table</title>
        </head>
        <body>
            <div class="container mt-4">
                <h2 class="mb-4">Variable Positions</h2>
                {table_html}
                <a href="/" class="btn btn-primary mt-3">Back to Home</a>
            </div>
        </body>
    </html>
    """



    
if __name__ == '__main__':
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)
