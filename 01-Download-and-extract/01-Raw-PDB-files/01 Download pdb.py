"""
Purpose of this script is to extract the list of pdb files required for download.
Output of this script are the downloaded pdb files from the chosen pdb website.
"""

# STEP 1
# Necessary packages
import os
from selenium import webdriver
import time
import csv

# STEP 2
main_folder = r'D:\PHML B factor estimation\02 Project'
# List of RNA names to download
folder = r'00 Basic information'
sources = ["Train_list", "Test_list"]
RNAproteins = []
for src in sources:
	pathway = main_folder + "\\" + folder + "\\" + src
	file = open(pathway, "r").read()
	file = file.split("\n")
	RNAproteins.extend(file)

# The PDB website only accept proteins name such as 1asy instead of the original 1asy_R.
# So the following step is to clean up the protein name.
RNAproteins = [protein for protein in RNAproteins if len(protein) > 1]
RNAproteins_clean = []
for protein in RNAproteins:
	temp = protein.split("_")[0]
	RNAproteins_clean.append(temp)


# STEP 3
# Web scraping
# The PDB website is open and the pdb file is downloaded for each of the RNA protein in the RNAproteins list.
# PDB website: https://www.rcsb.org/structure/<protein name>
website = r"https://files.rcsb.org/download/"
for protein in RNAproteins_clean:
	full_website = website + protein + ".pdb"
	driver = webdriver.Chrome()
	driver.get(full_website)
	time.sleep(15)
	driver.close()

# STEP 4
# Extract the list of downloaded pdb
folder = r'01 Download and extract\01 Raw PDB files'
downloaded_files = []
for file in os.listdir(main_folder + "\\" + folder):
    if file.endswith(".pdb"):
    	downloaded_files.append(file.split(".")[0])

# STEP 5
# Check for the missing files 
missing_files = [f for f in RNAproteins_clean if f not in downloaded_files]
with open("Missing pdb.csv", 'w') as myfile:
     writer = csv.writer(myfile, lineterminator='\n')
     for val in missing_files:
     	writer.writerow([val])
