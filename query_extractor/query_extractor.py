from setup import update_status

def extract_query(query, ask_llm):

    try:
        # PARSE QUERY
        update_status("Parsing query & identifying gene symbols...")
        plan_prompt = f"""
        Extract ONLY explicit gene symbols. If none, return NONE.
        Extract publication year if explicitly stated.
        Format: GENE|YEAR|TOPIC
        Query: "{query}"
        """
        plan = ask_llm(plan_prompt)
        
        try:
            parts = [p.strip() for p in plan.split("|")]
            gene = parts[0] if len(parts) > 0 else "NONE"
            year = parts[1] if len(parts) > 1 and parts[1] != "NONE" else None
            topic = parts[2] if len(parts) > 2 else query
        except:
            gene, year, topic = "NONE", None, query
                
        return gene, year, topic

    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        
        return "NONE", None, None
            