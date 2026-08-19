
# 📄 Document RAG Assistant

An end-to-end **Retrieval-Augmented Generation (RAG)** system for querying and analysing research PDFs, built with **Streamlit**.

The system combines **semantic retrieval (FAISS)** and **keyword retrieval (BM25)**, with optional **CrossEncoder reranking**, relevance-based hallucination control, multi-model answer generation, literature review synthesis, and a live **RAG performance and evaluation dashboard**.

---

## 🚀 Demo

Upload cancer research PDFs → ask questions → retrieve relevant passages → generate grounded answers with source attribution and retrieval/performance metrics.

---

## ✨ Features

| Feature | Details |
|---|---|
| **PDF Ingestion** | Multi-file upload with PyPDF text extraction |
| **Chunking** | RecursiveCharacterTextSplitter — 1000 characters, 200 overlap |
| **Embeddings** | `BAAI/bge-base-en-v1.5` via Sentence Transformers |
| **Semantic Retrieval** | FAISS `IndexFlatL2` |
| **Keyword Retrieval** | BM25 via `rank-bm25` |
| **Hybrid Retrieval** | Normalized score fusion with balanced cross-document selection |
| **Reranking** | Optional CrossEncoder `ms-marco-MiniLM-L-6-v2` |
| **Hallucination Control** | Relevance gate before LLM generation |
| **Answer Generation** | Gemini 2.5 Flash or local FLAN-T5 |
| **Literature Review** | Per-document summarization followed by corpus-level synthesis |
| **Evaluation Dashboard** | Retrieval, reranking, generation latency, similarity scores and source distribution |
| **Frontend** | Streamlit |

---

# 🏗️ Architecture

```text
                         User Question
                              │
                              ▼
                    ┌───────────────────┐
                    │  Hybrid Retrieval │
                    │                   │
                    │   FAISS + BM25    │
                    └─────────┬─────────┘
                              │
                              ▼
                    ┌───────────────────┐
                    │  Score Fusion     │
                    │                   │
                    │ Semantic + BM25    │
                    └─────────┬─────────┘
                              │
                              ▼
                    ┌───────────────────┐
                    │ Optional          │
                    │ CrossEncoder      │
                    │ Reranking         │
                    └─────────┬─────────┘
                              │
                              ▼
                    ┌───────────────────┐
                    │ Relevance Gate    │
                    └─────────┬─────────┘
                              │
                     ┌────────┴────────┐
                     │                 │
                  Relevant         Not Relevant
                     │                 │
                     ▼                 ▼
              ┌─────────────┐   ┌──────────────┐
              │ LLM         │   │ No Results   │
              │ Generation  │   │ Response     │
              └──────┬──────┘   └──────┬───────┘
                     │                 │
                     └────────┬────────┘
                              ▼
                    ┌───────────────────┐
                    │ Evaluation        │
                    │ & Metrics         │
                    └───────────────────┘
````

> **Note:** LangGraph-based orchestration is currently under development and is planned as a future extension of the pipeline.

---

# 📁 Project Structure

```text
document-rag-assistant/
│
├── app.py
├── requirements.txt
├── .env
├── .gitignore
├── README.md
│
├── src/
│   ├── config.py
│   ├── pdf_loader.py
│   ├── chunking.py
│   ├── embeddings.py
│   ├── vector_store.py
│   ├── retriever.py
│   ├── reranker.py
│   ├── llm.py
│   ├── pipeline.py
│   └── evaluator.py
│
├── data/
│   ├── uploads/
│   └── index/
│
├── notebooks/
│   └── cancer_rag_dev.ipynb
│
└── tests/
    └── test_pipeline.py
```

---

## Module Responsibilities

| File              | Purpose                                                 |
| ----------------- | ------------------------------------------------------- |
| `app.py`          | Streamlit user interface                                |
| `config.py`       | Central configuration and tunable parameters            |
| `pdf_loader.py`   | Extracts text from PDFs                                 |
| `chunking.py`     | Splits documents into overlapping chunks                |
| `embeddings.py`   | Loads embedding model and generates embeddings          |
| `vector_store.py` | Builds FAISS and BM25 indexes                           |
| `retriever.py`    | Performs hybrid retrieval                               |
| `reranker.py`     | Performs optional CrossEncoder reranking                |
| `llm.py`          | Handles Gemini / FLAN-T5 generation                     |
| `pipeline.py`     | Core RAG and literature-review pipelines                |
| `evaluator.py`    | Calculates and visualizes retrieval/performance metrics |

---

# 🔄 How It Works

## 1. Document Ingestion

Multiple research PDFs can be uploaded through the Streamlit interface.

```text
PDF
 ↓
PyPDF
 ↓
Extracted Text
```

---

## 2. Chunking

Documents are divided into smaller overlapping chunks using:

```text
RecursiveCharacterTextSplitter
```

Configuration:

```text
Chunk size:       1000 characters
Chunk overlap:     200 characters
```

The overlap helps preserve contextual information across chunk boundaries.

---

## 3. Embedding Generation

Each chunk is converted into a dense vector using:

```text
BAAI/bge-base-en-v1.5
```

via Sentence Transformers.

```text
Text Chunk
    ↓
Embedding Model
    ↓
Dense Vector
```

---

# 🔎 Retrieval

The system uses **two retrieval strategies**.

## Semantic Retrieval

FAISS is used to find semantically similar chunks.

```text
Question
   ↓
Embedding
   ↓
FAISS
   ↓
Semantic Candidates
```

Implementation:

```text
FAISS IndexFlatL2
```

---

## Keyword Retrieval

BM25 provides keyword-based retrieval.

```text
Question
   ↓
BM25
   ↓
Keyword Candidates
```

This is particularly useful when queries contain:

* Gene names
* Protein names
* Biomarkers
* Diseases
* Specific terminology
* Technical abbreviations

---

# 🔀 Hybrid Retrieval

Semantic and keyword retrieval scores are normalized and combined.

```text
Hybrid Score =
    0.5 × Semantic Score
  + 0.5 × BM25 Score
```

Configuration:

```python
SEMANTIC_WEIGHT = 0.5
BM25_WEIGHT = 0.5
```

The weights should sum to `1.0`.

The system also applies **balanced cross-document selection** to prevent a single uploaded PDF from dominating the retrieved context.

---

# 🎯 CrossEncoder Reranking

An optional second-stage reranker can be enabled.

```text
Query
  │
  ▼
FAISS + BM25
  │
  ▼
Candidate Pool
  │
  ▼
CrossEncoder
  │
  ▼
Final Top-K Chunks
```

Model:

```text
cross-encoder/ms-marco-MiniLM-L-6-v2
```

Configuration:

```python
ENABLE_RERANKER = True

TOP_K = 5

RERANKER_CANDIDATE_K = 20
```

This allows the system to retrieve a larger candidate pool first and then perform more precise relevance ranking.

---

# 🛡️ Hallucination Control

The system uses a **relevance gate** before calling the LLM.

If retrieved content does not meet the configured relevance threshold:

```text
Question
   ↓
Retrieval
   ↓
Relevant context?
   ├── YES → Generate Answer
   │
   └── NO  → Return "Information not found"
```

The LLM is therefore **not called when the retrieved documents are insufficient**, reducing the risk of unsupported answers.

Configuration:

```python
MAX_DISTANCE = 1.0
HYBRID_MAX_DISTANCE = 0.9
```

---

# 🤖 Answer Generation

Two generation options are supported.

## Gemini

```text
Gemini 2.5 Flash
```

Used through the Gemini API.

## Local Model

```text
google/flan-t5-small
```

This provides a local fallback option without requiring an external LLM API.

---

# 📚 Literature Review Generation

The system can generate a structured literature review from all uploaded research papers.

### Pipeline

```text
PDF 1 ──→ Summary 1 ──┐
PDF 2 ──→ Summary 2 ──┤
PDF 3 ──→ Summary 3 ──┼──→ Final Synthesis
PDF N ──→ Summary N ──┘
```

Each document is first summarized independently.

The summaries are then passed to Gemini to produce a consolidated review.

### Generated Sections

1. Executive Summary
2. Key Findings
3. Important Biomarkers / Concepts
4. Limitations
5. Conclusion

### API Usage

For `N` uploaded PDFs:

```text
N document-summary calls
+
1 final synthesis call
=
N + 1 API calls
```

The number of API calls therefore depends on the number of documents rather than the number of chunks.

---

# 📊 Evaluation & Performance Dashboard

The application includes a live dashboard for monitoring **RAG retrieval and pipeline performance**.

## Metrics

### Retrieval

* Retrieval latency
* Number of retrieved chunks
* Similarity scores
* Indexed chunk count
* Source document distribution
* Semantic vs keyword contribution

### Reranking

* Reranking latency
* Reranker contribution
* Final ranked chunks

### Generation

* Generation latency
* Answer length
* Prompt length
* Context length
* Model used

## Visualizations

The dashboard includes:

### Similarity Score per Retrieved Chunk

Bar chart showing the similarity/relevance score of each retrieved chunk.

### Retrieved Chunks by PDF

Pie chart showing how retrieved chunks are distributed across uploaded documents.

### Retrieval Contribution

Pie chart showing:

* Semantic retrieval
* Keyword retrieval
* Both

### Processing Time

Stacked bar chart showing:

```text
Retrieval → Reranking → Generation
```

These metrics provide visibility into **retrieval behavior, source attribution, and system performance**.

> **Note:** The current dashboard focuses primarily on retrieval analytics and pipeline performance. Ground-truth RAG quality metrics such as Context Precision, Context Recall, Faithfulness, and Answer Correctness can be added using a labeled evaluation dataset.

---

# ⚙️ Configuration

All major parameters are centralized in:

```text
src/config.py
```

## Reranking

```python
ENABLE_RERANKER = True
```

## Retrieval

```python
TOP_K = 5
RERANKER_CANDIDATE_K = 20
```

## Hybrid Weights

```python
SEMANTIC_WEIGHT = 0.5
BM25_WEIGHT = 0.5
```

The weights should sum to:

```text
1.0
```

## Relevance Gate

```python
MAX_DISTANCE = 1.0
HYBRID_MAX_DISTANCE = 0.9
```

## Models

```python
EMBEDDING_MODEL_NAME = "BAAI/bge-base-en-v1.5"

RERANKER_MODEL_NAME = \
    "cross-encoder/ms-marco-MiniLM-L-6-v2"

LLM_MODEL_NAME = \
    "google/flan-t5-small"
```

---

# 🚀 Quickstart

## 1. Clone the Repository

```bash
git clone https://github.com/<your-username>/document-rag-assistant.git

cd document-rag-assistant
```

## 2. Install Dependencies

```bash
pip install -r requirements.txt
```

## 3. Configure Gemini API

Create a `.env` file in the project root:

```text
GEMINI_API_KEY=your_key_here
```

The `.env` file should **never be committed to GitHub**.

## 4. Run the Application

```bash
streamlit run app.py
```


```

---

# 🧪 End-to-End Testing

Place at least one PDF inside:

```text
data/uploads/
```

Then run:

```bash
pytest tests/test_pipeline.py
```

The test executes the main pipeline:

```text
PDF Loading
    ↓
Chunking
    ↓
Embedding
    ↓
FAISS
    ↓
BM25
    ↓
Hybrid Retrieval
    ↓
Answer Generation
    ↓
Evaluation
```

---

# 📦 Requirements

```text
pypdf
langchain-text-splitters
sentence-transformers
faiss-cpu
transformers
sentencepiece
torch
numpy
rank-bm25
plotly
typing_extensions
google-genai
python-dotenv
streamlit
```

Install with:

```bash
pip install -r requirements.txt
```

---

# 🔬 Example Use Cases

The system can be used for research-oriented questions such as:

```text
What are the major biomarkers associated with cancer progression?

Compare the findings of these research papers.

What limitations were reported across the studies?

Which papers discuss targeted therapy?

Summarize the major findings across all uploaded documents.
```

The system answers using the uploaded research documents rather than relying solely on the LLM's pretrained knowledge.

---

# 🔐 Security

Sensitive configuration files and runtime data are excluded from version control.

`.gitignore`:

```text
# Environment
.env

# Runtime data
data/uploads/
data/index/

# Python
__pycache__/
*.pyc
*.pyo

# Virtual environments
.venv/
venv/

# Jupyter
.ipynb_checkpoints/

# OS
.DS_Store
Thumbs.db
```

**Never commit your Gemini API key.**

---

# 🛠️ Future Improvements

Potential extensions include:

* [ ] LangGraph-based workflow orchestration
* [ ] RAGAS-based evaluation
* [ ] Ground-truth evaluation dataset
* [ ] Context Precision
* [ ] Context Recall
* [ ] Faithfulness
* [ ] Answer Relevancy
* [ ] Answer Correctness
* [ ] Automated evaluation over 50–100 benchmark questions
* [ ] Retrieval quality comparison: FAISS vs BM25 vs Hybrid
* [ ] Reranker performance comparison
* [ ] Persistent vector database
* [ ] Metadata filtering
* [ ] PDF page-level citations
* [ ] Conversation memory
* [ ] Streaming responses
* [ ] Experiment tracking
* [ ] Evaluation history across model/configuration changes

---

# 📌 Project Highlights

This project demonstrates an end-to-end RAG architecture rather than a simple LLM chatbot.

The current system combines:

```text
                 ┌──────────────┐
                 │ Research PDF │
                 └──────┬───────┘
                        ↓
                ┌───────────────┐
                │   Chunking    │
                └───────┬───────┘
                        ↓
              ┌───────────────────┐
              │    FAISS + BM25   │
              │   Hybrid Search   │
              └─────────┬─────────┘
                        ↓
                ┌───────────────┐
                │  CrossEncoder  │
                │   Reranking    │
                └───────┬────────┘
                        ↓
                ┌───────────────┐
                │ Relevance Gate│
                └───────┬───────┘
                        ↓
                ┌───────────────┐
                │ Gemini / T5   │
                └───────┬───────┘
                        ↓
                ┌───────────────┐
                │ Grounded      │
                │ Answer        │
                └───────────────┘
```

---

# 📈 Project Status



### Implemented

* ✅ Multi-PDF ingestion
* ✅ PDF text extraction
* ✅ Document chunking
* ✅ BGE embeddings
* ✅ FAISS semantic retrieval
* ✅ BM25 keyword retrieval
* ✅ Hybrid retrieval
* ✅ Balanced cross-document retrieval
* ✅ Optional CrossEncoder reranking
* ✅ Relevance-based hallucination control
* ✅ Gemini / FLAN-T5 generation
* ✅ Literature review generation
* ✅ Streamlit interface
* ✅ Retrieval and performance dashboard
* ✅ End-to-end pipeline testing



