## This file scans the PA table, uses the GCA (contig) accession number to lookup on the ENA,
## then matches using the ENA name ID to match with records from the GPS metadata table.
import pandas as pd
import xmltodict, json
import requests ## for getting the xml from ENA
import time

## Metadata from GPS
df = pd.read_csv("outputs_from_GPS/monocle-metadata.csv", index_col="Lane_id") ### strain metadata
pa = pd.read_csv("outputs_from_GPS/labeled_pa_table_85pid_blpA_cleaned.csv", delimiter="\t", index_col="ID") ## table with gene info we want
full_pa = pd.read_csv("outputs_from_GPS/labeled_pa_table.csv", delimiter= ",", index_col="Strain") ## table with strain and GCA IDs
pa_gcas = full_pa[["Isolate","GCA_ID"]]
pa = pa.join(pa_gcas)
## Creating a GCA to ENA ID access
## iterating over each row, finding and storing the ID from ENA

ENA_names = []
GCAs = []
i=0

ENA_tot = len(set(pa.GCA_ID))
for GCA_id in set(pa.GCA_ID):

    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + GCA_id
    response = requests.get(url)
    decoded = response.content.decode('utf-8')
    ENA_name = xmltodict.parse(decoded)["ASSEMBLY_SET"]["ASSEMBLY"]["NAME"]
    ENA_names.append(ENA_name)
    GCAs.append(GCA_id)

    if i % 40 == 0:
        print("Now checking GCA ", i + 1, " of  ", ENA_tot, "   ", GCA_id) 
    i += 1
    #time.sleep(0.03) ## avoiding the 50 per second limit


ENA_lanes_df = pd.DataFrame(data = {"ENA_lane_id" : ENA_names, "GCA_ID" : GCAs})
pa = pa.merge(ENA_lanes_df, on = "GCA_ID", how = "left")
pa = pa.set_index("ENA_lane_id")
out_pa = pa.join(df)
out_pa.to_csv("outputs_from_GPS/PA_with_metadata.csv")
