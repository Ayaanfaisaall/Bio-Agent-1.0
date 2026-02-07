from setup import perform_setup
from query_extractor import extract_query
from retriever import ncbi_retrieve, pubmed_retrieve
from validator import validate_pubmed
from generator import generate_report
from safety_auditor import audit_safety
from save_file import save_report

def main():
    
    ask_llm = perform_setup()
    
    print("\nBio-Agent 1.0")
    print("------------------------------------------------------------")
    print("Evidence-Based Bioinformatics Data Retrieval")
    print("------------------------------------------------------------")

    query = input("\nEnter your promt containing gene name (e.g., TP53), biological topic or particular year:\n")

    if query:
        
        gene, year, topic = extract_query(query, ask_llm)
        
        if gene == "NONE" and topic is None:
            print("Provide a gene name or topic name or both to start!")
            return

        context, aliases, is_gene_valid = ncbi_retrieve(gene)
        
        articles = pubmed_retrieve(gene, topic, year, is_gene_valid)
        
        validation_note = validate_pubmed(articles, gene, aliases, topic, year)
        
        report = generate_report(query, ask_llm, context, validation_note, articles)
        
        audit_safety(context, report, ask_llm)
        
        save_report(query, report)

    else:
        print("\nReady to retrieve real-time data from NCBI & PubMed.")

if __name__ == "__main__":
    main()