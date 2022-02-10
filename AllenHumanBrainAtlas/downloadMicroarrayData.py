#Packages
import os
import abagen

#Directory to store data
data_dir = os.getcwd()+'/Data/microarray'

#If the "microarray" directory does not exist, create it
#abagen will download to ~ if this does not exist
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
   
files = abagen.fetch_microarray(donors = 'all', data_dir = data_dir)