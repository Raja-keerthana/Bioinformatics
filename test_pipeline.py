from src.pdf_loader import extract_text_from_pdf
from src.chunking import build_chunks
from src.embeddings import load_embedding_model, generate_embeddings
from src.vector_store import build_faiss_index, build_metadata_store, save_index
from src.llm import load_llm
from src.pipeline import run_rag_pipeline

# Build the full index 
with open("data/uploads/paper_1.pdf", "rb") as f:
    pdf_result_1 = extract_text_from_pdf(f.read(), "paper_1.pdf")
with open("data/uploads/paper_2.pdf", "rb") as f:
    pdf_result_2 = extract_text_from_pdf(f.read(), "paper_2.pdf")

chunks = build_chunks([pdf_result_1, pdf_result_2])
embed_model = load_embedding_model()
embeddings = generate_embeddings(chunks, embed_model)

index = build_faiss_index(embeddings)
metadata = build_metadata_store(chunks)
save_index(index, metadata)  

llm_pipeline = load_llm()  

# Test 1
result = run_rag_pipeline(
    "What biomarkers are used in breast cancer?",
    embed_model, index, metadata, llm_pipeline
)
print("Q:", result["question"])
print("Sources:", [(s["source"], round(s["score"], 4)) for s in result["sources"]])
print("A:", result["answer"])
print()

# Test 2
result2 = run_rag_pipeline(
    "What is the treatment protocol for pancreatic cancer in elderly patients?",
    embed_model, index, metadata, llm_pipeline
)
print("Q:", result2["question"])
print("Sources:", [(s["source"], round(s["score"], 4)) for s in result2["sources"]])
print("A:", result2["answer"])
assert "not found" in result2["answer"].lower()