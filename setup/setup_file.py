import sys
import os
from Bio import Entrez
from dotenv import load_dotenv
from llm_config import configure_llm

def perform_setup():    
    # LOAD ENVIRONMENT VARIABLES
    load_dotenv("api_email.env")

    # API SETUP & CONFIGURATION
    API_KEY = os.getenv("GEMINI_API_KEY")
    ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL")

    # Set Entrez Email globally
    if not ENTREZ_EMAIL:
        print("[WARNING] ENTREZ_EMAIL not found in .env file.")
    else:
        Entrez.email = ENTREZ_EMAIL

    # INITIALIZE LLM
    # This creates the 'ask_llm' function used throughout the project.
    # To swap models, simply change the parameters here.
    ask_llm = configure_llm(API_KEY, model_name="gemini-2.5-flash")

    if not ask_llm:
        print("[ERROR] Failed to configure LLM. Exiting.")
        sys.exit()
        
    else:
        return ask_llm

def update_status(msg):
    print(f"[STATUS] {msg}")