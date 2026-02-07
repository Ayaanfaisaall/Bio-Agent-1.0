# Bio-Agent 1.0

LLM-powered, evidence-based bioinformatics data extraction from NCBI Gene and PubMed.

Bio-Agent 1.0 is a Retrieval-Augmented Generation pipeline designed to reduce hallucinations when working with biological and genomic literature. It connects natural language queries to authoritative databases and uses an LLM only for grounded synthesis.

---

## Why This Project Exists

Bioinformatics research faces two common problems.

* Too much literature to read manually.
* AI models that generate confident but unsupported biological claims.

Bio-Agent 1.0 addresses both by enforcing a strict pipeline.

* Data is retrieved first.
* Genes are validated before use.
* Abstracts are checked for relevance.
* The LLM is restricted to the retrieved context only.

Every output is tied to real PubMed records.

---

## Key Capabilities

* Natural language intent parsing.
* Gene symbol normalization and validation using NCBI Gene.
* PubMed retrieval with year and topic filters.
* Soft validation of abstracts against gene aliases and keywords.
* Structured report generation using retrieved evidence only.
* Automated safety auditing to detect overclaims.
* Local report saving for reproducibility.

---

## System Flow

The system follows a linear and auditable pipeline.

1. User submits a natural language query.
2. Intent is extracted into gene, topic, and year.
3. Gene validity and aliases are retrieved from NCBI.
4. PubMed abstracts are fetched and parsed.
5. Abstracts are soft-validated for relevance.
6. A structured report is generated.
7. A safety auditor checks for unsupported claims.
8. The final report is saved locally.

---

## Project Structure

```
bio-agent-1.0/
│
├── setup/
│   ├── __init__.py
│   └── setup_file.py
│
├── llm_config/
│   ├── __init__.py
│   └── llm_config.py
│
├── query_extractor/
│   ├── __init__.py
│   ├── extractor.py
│   └── README.md
│
├── retriever/
│   ├── __init__.py
│   ├── ncbi.py
│   ├── pubmed.py
│   └── README.md
│
├── validator/
│   ├── __init__.py
│   ├── pubmed_validator.py
│   └── README.md
│
├── generator/
│   ├── __init__.py
│   ├── report_generator.py
│   └── README.md
│
├── safety_auditor/
│   ├── __init__.py
│   ├── auditor.py
│   └── README.md
│
├── save_file/
│   ├── __init__.py
│   └── save_report.py
│
├── main.py
├── requirements.txt
├── api_email.env
└── README.md
```

---

## Module Overview

### setup

Handles environment setup.

* Loads API keys and environment variables.
* Initializes required services.

### llm_config

Configures the LLM client.

* API initialization.
* Retry logic and error handling.
* Enforces deterministic, constrained generation.

### query_extractor

Parses user input.

* Extracts gene symbols.
* Detects biological topics.
* Identifies publication years or ranges.

### retriever

Responsible for external data access.

* `ncbi.py` validates genes and retrieves official summaries and aliases.
* `pubmed.py` queries PubMed and parses XML abstracts into structured data.

### validator

Performs soft validation.

* Checks abstracts for gene aliases.
* Confirms topic relevance.
* Flags weak or irrelevant evidence.

### generator

Creates the final report.

* Uses only validated context.
* Produces structured sections such as summaries and mechanisms.

### safety_auditor

Audits the generated output.

* Compares claims against source abstracts.
* Flags unsupported or exaggerated statements.

### save_file

Manages output storage.

* Saves reports locally.
* Uses query-based file naming.

---

## Installation

### Requirements

* Python 3.8 or higher
* Internet access
* NCBI Entrez email
* LLM API access

### Setup

Clone the repository.

```bash
git clone https://github.com/Ayaanfaisaall/Bio-Agent-1.0.git
cd Bio-Agent-1.0
```

Install dependencies.

```bash
pip install -r requirements.txt
```

Configure credentials.

Open `api_email.env` and add:

```env
GEMINI_API_KEY=your_api_key_here
ENTREZ_EMAIL=your_email@example.com
```

---

## Usage

Run the application.

```bash
python main.py
```

Example input.

```
How does TP53 contribute to lung cancer between 2019 and 2020?
```

The system will:

* Validate TP53.
* Retrieve relevant PubMed abstracts.
* Generate a grounded report.
* Perform a safety audit.
* Save the output locally.

---

## Output

Reports are saved automatically.

Example filename:

```
TP53_lung_cancer_2019_2020.md
```

Each report includes:

* Executive summary
* Biological mechanisms
* Evidence-backed findings
* Safety audit result

---

## Limitations

* Uses abstracts only. No full-text mining.
* Soft validation is keyword based.
* Not intended for clinical decision-making.

---

## Disclaimer

This project is for research and educational use only.

Users must verify all biological and medical claims using the original literature.

---

## License

No License.

