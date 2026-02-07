import os
import re
from datetime import datetime

def save_report(query, report):
    """
    Saves the report to a Markdown file in a 'reports' directory.
    """
    # Create reports directory if it doesn't exist
    if not os.path.exists("reports"):
        os.makedirs("reports")

    # Sanitize filename (replace spaces with underscores, remove weird chars)
    # e.g., "TP53 lung cancer" -> "TP53_lung_cancer"
    clean_query = re.sub(r'[^\w\s-]', '', query).strip().replace(' ', '_')
    
    # Add timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"reports/{clean_query}_{timestamp}.md"

    try:
        with open(filename, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"\n[SUCCESS] Report saved to: {filename}")
    except Exception as e:
        print(f"\n[ERROR] Failed to save report: {e}")