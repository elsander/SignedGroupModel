GetTaxonomy.py
--------------
Create taxonomic records for a species list 
using the ITIS database.
   
written by Stefano Allesina (sallesina@uchicago.edu)
october 25, 2013. Version 0.9.0
   
Requirements:
--------------
python with the libraries 
- sys 
- sqlite3

Usage:
------
python GetTaxonomy.py SpeciesList.txt Where/Is/The/File/ITIS.sqlite

Input:
-------
A species list: text file with one record per line.
Each record can be:
- the latin binomial for the species
- or the genus name 
- or the family name

A path to the ITIS database (SQLite version)
This file can be downloaded from the page
http://www.itis.gov/downloads/index.html
Choose the SQLite version of the file, 
and unzip the file somewhere (about 1/2 Gb)

Output:
-------
A comma-separated text file .taxon containing all the taxonomic info
The file name is obtained by adding .taxon to the path of the input file

