import sys
import time
import random
from Bio import Entrez
from setup import update_status

#BACKEND FUNCTION
def get_pubmed_data_structured(topic, year=None, retmax=10):
    """Fetch structured PubMed info (titles, abstracts, publication year)"""
    try:
        term = topic
        if year and year.isdigit():
            term += f" AND {year}[PDAT]"
            
        time.sleep(random.uniform(1, 3))

        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return []
        
        time.sleep(random.uniform(1, 3))

        fetch = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="xml")
        articles = Entrez.read(fetch)
        fetch.close()

        structured = []
        for art in articles["PubmedArticle"]:
            try:
                article_meta = art.get("MedlineCitation", {}).get("Article", {})
                # Title
                title = article_meta.get("ArticleTitle", "")
                # Abstract
                abstract_list = article_meta.get("Abstract", {}).get("AbstractText", [""])
                abstract = " ".join(abstract_list) if abstract_list else ""
                
                # Publication year logic
                pub_year = ""
                journal = article_meta.get("Journal", {})
                pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
                
                if "Year" in pub_date:
                    pub_year = pub_date["Year"]
                elif "MedlineDate" in pub_date:
                    pub_year = pub_date["MedlineDate"].split()[0]
                else:
                    pub_year = art.get("MedlineCitation", {}).get("DateCompleted", {}).get("Year", "")
                    if not pub_year:
                        pub_year = art.get("MedlineCitation", {}).get("DateCreated", {}).get("Year", "")

                structured.append({"title": title, "abstract": abstract, "year": pub_year})
            except:
                continue
        return structured
    except:
        return []
    

def pubmed_retrieve(gene, topic, year, is_gene_valid):
    
    try:
        update_status("Retrieving literature from PubMed...")
        final_topic = f"{gene} AND {topic}" if is_gene_valid else topic
        articles = get_pubmed_data_structured(final_topic, year, retmax=20 if is_gene_valid else 10)

        if not articles:
            print(f"\n[ERROR] No articles found for {final_topic}.")
            sys.exit()
                    
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        
    return articles

