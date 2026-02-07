import time
from google import genai

def configure_llm(api_key, model_name="gemini-2.5-flash"):
    """
    Configures the LLM client and returns a generation function.
    This allows swapping the model or provider easily by changing this config.
    """
    try:
        # Initialize client with the provided key
        client = genai.Client(api_key=api_key)
    except Exception as e:
        print(f"API Key Error: {e}")
        return None

    def ask(prompt):
        """Inner function to call the configured LLM with retry logic"""
        for attempt in range(3):
            try:
                response = client.models.generate_content(
                    model=model_name,
                    contents=prompt
                )
                return response.text.strip()
            except Exception as e:
                if "429" in str(e): time.sleep(10)
                elif "503" in str(e): time.sleep(2)
                else: raise e
        raise Exception("API Busy. Try again later.")

    return ask
