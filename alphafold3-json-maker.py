## use this code to turn a fastafile into a bunch of individual jsons for Alphafold analysis

import json
from Bio import SeqIO
from pathlib import Path

fasta_file = ""
fasta_dict =  SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
list_of_dicts_for_json = []

def make_jsons(fasta_dict):
    for key, value in fasta_dict.items():
        new_dict = {
            'name':fasta_dict[key].name,
            'sequences':[
                {
                    'protein':{
                    'id':"A",
                    'sequence': 
                        str(fasta_dict[key].seq.replace('*',''))
                  }
                }
            ],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1
        }
            
        with open(fasta_dict[key].name+".json", mode="w", encoding="utf-8") as write_file:
            json.dump(new_dict, write_file)
        write_file.close()
    return
