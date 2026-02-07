from setup import update_status

def generate_report(query, ask_llm, context, validation_note, articles):
    
    report = "Report generation failed."
    
    try:
        # Prepare context
        context += validation_note
        for art in articles:
            context += f"\n[TITLE]\n{art['title']}\n[ABSTRACT]\n{art['abstract']}\n[YEAR] {art['year']}\n"

        # REPORT GENERATION
        update_status("Synthesizing evidence-bound report...")
        
        # Adjusted prompt to request Markdown instead of HTML for CLI readability
        report_prompt = f"""
        You are a Senior Bioinformatics Analyst.
        
        ### DATA SOURCE:
        {context}
        
        ### INSTRUCTIONS:
        1. **Evidence-Based:** Use ONLY the provided data. If a mechanism is not in the text, state "Insufficient evidence in retrieved abstracts."
        2. **No Hallucinations:** Do NOT speculate on pathways not mentioned.
        3. **Formatting:** Use Standard Text Formatting (Headings, Bullet points). Do NOT use HTML.
        
        ### REPORT STRUCTURE:
        1. Executive Summary
        (Brief overview of the gene/topic context)
        
        2. Biological Mechanism
        (Detailed molecular interactions based on the data)
        
        3. Key Research Findings
        (Bulleted list of specific outcomes from the abstracts)
        
        4. Study Limitations
        (Sample sizes, specific years, or lack of clinical trials mentioned)
        
        5. References
        (List titles and years)
        
        Query: {query}
        """

        report = ask_llm(report_prompt)
                
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        
    return report

            
