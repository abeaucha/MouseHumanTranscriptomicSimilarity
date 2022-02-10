# downloadMicroarrayData.py 
#
# Download AHBA microarray data from the web API
# 
# Antoine Beauchamp
# Created: February 10th, 2022
# Edited: February 10th, 2022

#Packages
import os
import abagen

#Directory to store the data
data_dir = os.getcwd()+'/Data/microarray'

#If the "microarray" sub-directory does not exist, create it
#abagen will download to ~ otherwise
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

#Download microarray data
files = abagen.fetch_microarray(donors = 'all', data_dir = data_dir)

#Download AHBA hierarchical ontology as JSON
hierarchy_url = "http://api.brain-map.org/api/v2/structure_graph_download/10.json"
hierarchy_url_get = requests.get(hierarchy_url)
hierarchy_dict = json.loads(hierarchy_url_get.text)
with open("Data/AHBA_hierarchy_definitions.json", "w") as outfile: 
    json.dump(hierarchy_dict, outfile)