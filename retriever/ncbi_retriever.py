import sys
import time
import random
from Bio import Entrez
from setup import update_status

# BACKEND FUNCTIONS
def get_gene_data(gene):
    """Fetch Gene Summary + Aliases (Synonyms)"""
    gene = gene.strip().upper()
    try:
        time.sleep(random.uniform(1, 3))
        handle = Entrez.esearch(db="gene", term=f"{gene}[Gene Name] AND Homo sapiens[Organism]")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return None, []
        
        time.sleep(random.uniform(1, 3))

        fetch = Entrez.efetch(db="gene", id=record["IdList"][0], retmode="xml")
        data = Entrez.read(fetch)
        fetch.close()

        summary = data[0].get("Entrezgene_summary", "")
        aliases = []
        try:
            aliases = data[0]["Entrezgene_gene"]["Gene-ref"].get("Gene-ref_syn", [])
        except: pass
        return summary, aliases
    except Exception as e:
        return None, []

def ncbi_retrieve(gene):
    
    try:
        update_status(f"Validating Gene '{gene}' via NCBI...")
        context = ""
        is_gene_valid = False
        aliases = []

        if gene != "NONE":
            summary, aliases = get_gene_data(gene)
            if summary:
                context += f"\n[NCBI GENE SUMMARY]\n{summary}\n"
                is_gene_valid = True
            else:
                # Gene invalid, stop pipeline with error
                print(f"\n[ERROR] Gene '{gene}' not found in NCBI. Please check the input.")
                sys.exit()
                    
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
    
    return context, aliases, is_gene_valid
            
