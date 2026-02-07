import string
from setup import update_status

#BACKEND FUNCTION
def validate_pubmed_structured(articles, gene, aliases, topic, year):
    """
    Soft Validation: Checks for keywords but returns a list of warnings 
    instead of blocking the process.
    """
    missing_criteria = []
    
    # Clean topic: Remove punctuation and split into significant words
    translator = str.maketrans('', '', string.punctuation)
    clean_topic = topic.translate(translator).lower()
    # Filter out small stop words
    topic_words = [w for w in clean_topic.split() if len(w) > 3]

    genes_to_check = [g.lower() for g in [gene] + aliases if g != "NONE"]
    
    gene_found = False
    topic_found = False
    year_found = False

    for art in articles:
        # Combine title and abstract for search
        text = (art["title"] + " " + art["abstract"]).lower()
        
        # 1. Flexible Gene Check
        if gene != "NONE":
            # Check if any gene alias exists in the text
            if any(g in text for g in genes_to_check):
                gene_found = True
        else:
            gene_found = True 
        
        # 2. Flexible Topic Check
        if not topic_words or any(w in text for w in topic_words):
            topic_found = True

        # 3. Year Check
        if year:
            if str(year) in str(art["year"]):
                year_found = True
        else:
            year_found = True

    # Compile warnings instead of blocking
    if gene != "NONE" and not gene_found:
        missing_criteria.append(f"Gene '{gene}' or aliases not explicitly mentioned.")
    
    if not topic_found and topic_words:
        missing_criteria.append(f"Topic keywords '{topic}' not explicitly found.")
        
    if year and not year_found:
        missing_criteria.append(f"Publication Year '{year}' not found.")

    return missing_criteria



def validate_pubmed(articles, gene, aliases, topic, year):
    
    validation_note = ""
    
    try:
        missing_criteria = validate_pubmed_structured(articles, gene, aliases, topic, year)
            
        if missing_criteria:
            warning_msg = "\n".join(missing_criteria)
            print(f"\n[NOTE ON RELEVANCE]")
            print("The retrieved papers might be indirectly related. We couldn't find exact matches for:")
            print(f"{warning_msg}")
            print("Proceeding with analysis based on available context.\n")
            validation_note = f"\n[SYSTEM NOTE]: The user searched for terms that were not explicitly found in these abstracts. Please mention in 'Study Limitations' if the evidence is indirect.\n"
        else:
            update_status("Validation passed. Generating report...")
                    
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        
    return validation_note    