# Resume Reposity

The following repository is a small amount of analysis done on the 
Colorado Air Monitoring Mobile Lab (CAMML) data found at: 

https://www.colorado.gov/airquality/tech_doc_repository.aspx#camml_data

The CAMMLclass.py file has been written to intake and clean the data found in 
the Data folder under Extraction_Livingston_Data. This data contains both values
and flags which are explained in the CAMML GC VOC Data Flag Key v2.pdf. The data 
has been cleaned and values below the method detection limit (MDL) have been 
replaced as according to Polissar et al. (1998) which has been used for MDLs within
Positive Matrix Factorization (PMF) analysis. The exact MDL values used within this
analysis can be found in the CAMML_Teir1_202008_MDLs.xlsx document. 

The CAMMLanalysis.ipynb document is where the analysis is done. This analysis walks 
through some short analysis of the "other" data -- which includes ozone, particulate
matter, nitrogen oxides, and others, VOC data, and merged data. 


### Dependencies
* pandas
* numpy
* seaborn
* matplotlib







