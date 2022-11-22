import streamlit as st
import sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem import AllChem
from rdkit import DataStructs

import random
import string
import pandas as pd
import csv
from rdkit.Chem import rdMolDescriptors as Desc

name = ''.join(random.choices(string.ascii_uppercase + string.ascii_uppercase, k = 4))
st.markdown("<h1 style='text-align: center'>Molecular properties predictor</h1>", unsafe_allow_html=True)
def get2D(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    structure =Draw.MolToImageFile(mol=mol, filename="pictures/"+name+".jpeg")
    return name+'.jpeg'

def getDescriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    tpsa = Desc.CalcTPSA(mol)
    acceptors = Descriptors.NumHAcceptors(mol)
    donors  = Descriptors.NumHDonors(mol)
    logp = Descriptors.MolLogP(mol)
    return mw, acceptors, tpsa,  donors, logp

@st.cache
def generateCsv(mw, acceptors, donor, logp, tpsa):
    
    data = {
        "Molecular weight": mw,
        "acceptors": acceptors,
        "donor": donor,
        "logP": logp,
        "tpsa": tpsa
    }
    df = pd.DataFrame(data, index=["1"])
    csv = df.to_csv(index=False).encode('utf-8')
    return csv


form =st.form("form")
input = form.text_input("Enter the smiles of the compound you want to see its properties: ")
submitted = form.form_submit_button("Submit")

if submitted:
    filename = get2D(input)
    st.image("pictures/"+filename)
    mw,acceptors, tpsa, donors, logp = getDescriptors(input)

    
    col1, col2 = st.columns(2)
    col1.write(f"Molecular Weight: {mw}")
    col2.write(f"Number of H Acceptors: {acceptors}")
    col3, col4 = st.columns(2)
    col3.write(f"Number of H Donors: {donors}")
    col3.write(f"LogP: {logp}")
    col3.write(f"TPSA: {tpsa}")
    csv = generateCsv(mw, acceptors, donors, logp, tpsa)
    st.download_button("Download", file_name=name+".csv", data=csv, mime='text/csv')
   
    


