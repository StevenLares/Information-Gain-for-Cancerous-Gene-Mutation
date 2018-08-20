# -*- coding: utf-8 -*-
"""
Created on Wed Apr  19 01:22:47 2018

@author: Steven Lares
"""



import os
from pathlib import Path
import pandas as pd
import cx_Oracle
import numpy as np

saveDir = str(Path(os.getcwd())) + "\\" #can save at the current working directory

user = "REDACTED"
password = "REDACTED"
host = "REDACTED"
port = "REDACTED"
sid = "REDACTED"



#dsn and connection establish connection to server, should stay fixed.
dsn = cx_Oracle.makedsn(host = host, port = port, sid = sid)
connection = cx_Oracle.connect(user=user, password=password, dsn=dsn)



#Mostly used as a way to work with data in the host language
#query: String containing Oracle SQL query
def convertToDataframe(query):
    csr = connection.cursor() #creates a cursor
    csr.execute(query) #exectutes user given query
    data = [] #will hold tuples given by Oracle query
    data.append(csr.fetchall()) #adds the tuples to data list
    attributes = [] #will hold attributes from query
    
    #adds the attributes
    for desc in csr.description:
        attributes.append(desc[0])
    
    dataMatrix = np.array(data[0])
    
    #creates a dataframe holding the data and appropriate attributes
    #closes cursor, and returns the dataframe holding the data
    resultDF = pd.DataFrame(dataMatrix, columns = attributes)
    csr.close()
    return(resultDF)
    


#Returns expected inforation needed to encode a particular split
#Works regardless of number of splits
#Input: classCounts: Simply the number of entries in a particular class,
#like living or deceased
def Info(classCounts):
    total = sum(classCounts) #sums them all up, for easier calculation later
    
    info = 0 #will maintain running total
    
    for classCount in classCounts:
        ratio = classCount / total #the ratio that is used in calculations
        if ratio != 0:
            info = info + (-ratio * np.log2(abs(ratio))) #part of entropy calculation
            
    return(info) #entropy





IG_READY = convertToDataframe("SELECT * FROM IG_READY") #grabs data from database, just used for reference here
IG_READY = IG_READY.apply(pd.to_numeric, errors='ignore') #converts string values to ints where appropriate







#START of dictionary and table creation --------------------------------------------------------------


IG_READY_attributes = IG_READY.columns.values.tolist() #Gives list of columns

totalSamples = convertToDataframe("SELECT count(*) FROM IG_READY").iloc[0,0] #total number of samples in data


numDeceased = convertToDataframe("""
SELECT count(*)
FROM IG_READY
WHERE STATUS = 1                                 
""").iloc[0,0]


numLiving = convertToDataframe("""
SELECT count(*)
FROM IG_READY
WHERE STATUS = 0                                 
""").iloc[0,0]


Info_D = Info([numDeceased, numLiving]) #needed for first part of info gain.



#Below are the queries that will be used, per gene, to find CDT's
#They will be used in info gain calculations later
#A
true_positive_query = """
SELECT count(*) FROM
(
SELECT * FROM IG_READY
WHERE STATUS = 1

INTERSECT

SELECT * FROM IG_READY
WHERE [INSERT GENE] = 1
)                               
"""


#B
false_negative_query = """
SELECT count(*) FROM
(
SELECT * FROM IG_READY
WHERE STATUS = 1

MINUS

SELECT * FROM IG_READY
WHERE [INSERT GENE] = 1
)                               
"""


#C
false_positive_query = """
SELECT count(*) FROM
(
SELECT * FROM IG_READY
WHERE [INSERT GENE] = 1

MINUS

SELECT * FROM IG_READY
WHERE STATUS = 1
)                                                                
"""




#will hold all dataframes for each gene
CDT_dict = {}


#Holds formatting for problem A and B.
persistentDF = pd.DataFrame() 

#Used to make sure dataframes have correct index
indexer = 0

#goes through each gene in dataset, will create CDT's for each one.
#each CDT is saved as a dataframe, and stored in a dictionary for easy retrieving
for gene in IG_READY_attributes:
    if gene not in {"PATIENT_D", "STATUS"}:
       
        #A
        numTruePositives = convertToDataframe(true_positive_query.replace("[INSERT GENE]", gene)).iloc[0,0]
        
        
        #B
        numFalseNegatives = convertToDataframe(false_negative_query.replace("[INSERT GENE]", gene)).iloc[0,0]
        
        
        #C
        numFalsePositives = convertToDataframe(false_positive_query.replace("[INSERT GENE]", gene)).iloc[0,0]
        
        
        #D
        numTrueNegatives = totalSamples - (numTruePositives + numFalseNegatives + numFalsePositives)
        
        
        #Setting up CDT with appropriate labels
        geneCDT = pd.DataFrame(np.array([[numTruePositives, numFalsePositives],[numFalseNegatives, numTrueNegatives]]), columns = ["Deceased", "Living"], index = ["Mutation: 1","Mutation: 0"])
        
        CDT_dict.update({gene : geneCDT}) #Add this gene's CDT to dictionary for easy storage



        #Will calculate Info Gain for this gene while constructing
        #CDT dictionary in this loop
        #formula will be Gain(Gene) = Info_D - Info_Gene

        #mutationYRatio = ratio of positive mutation matches to total samples
        #mutationNRatio = ratio of negative mutation matches to total samples
        
        mutationYRatio = (geneCDT.loc["Mutation: 1"].sum()) / geneCDT.sum().sum()
        mutationNRatio = (geneCDT.loc["Mutation: 0"].sum()) / geneCDT.sum().sum()
        
        
        #Info_Gene is info needed to classify this gene
        Info_Gene =  mutationYRatio * Info(geneCDT.loc["Mutation: 1"]) + mutationNRatio * Info(geneCDT.loc["Mutation: 0"])
        
        #Application of Info Gain formula
        Info_Gain = Info_D - Info_Gene
        
        
        #create a one row dataframe for this gene. Contains data for problem A and B
        tempDF = pd.DataFrame(data = [gene, Info_Gain, numTruePositives], index = ["Gene ID", "IG", "O CNT"], columns = [indexer]).transpose()
        
        #union a persistent dataframe to have all data available at the end of this loop
        persistentDF = pd.concat([persistentDF, tempDF]) 
        
        
        #keeps the indexes of dataframe in consecutive order
        indexer = indexer + 1

#END of dictionary and table creation --------------------------------------------------------------







#Table for Problem B, sorted by IG
top5withO_CNT = persistentDF.sort_values(by = "IG", ascending = False).head(5)

#Created for Project Write Up as response to Problem B, sorted by O CNT and commented out
#top5sortedByO_CNT = persistentDF.sort_values(by = "O CNT", ascending = False).head(5)

#Table for Problem A, sorted by IG
top5withoutO_CNT = top5withO_CNT[["Gene ID", "IG"]]



        
#save both to csv

top5withoutO_CNT.to_csv(saveDir + "Problem A Table.csv", index = False)

top5withO_CNT.to_csv(saveDir + "Problem B Table.csv", index = False)




#Finally close the connection
connection.close()
