
#On MICe machines: Remove loaded modules
#Comment out if unnecessary
module purge

#Create the virtual environment
python3 -m venv AHBAEnv

#Activate the virtual environment
source AHBAEnv/bin/activate

#Installed necessary python packages
pip3 install -r PythonRequirements.txt

#Deactivate virtual environment
deactivate

