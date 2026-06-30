

from src.config import TOP_K, MAX_DISTANCE
from src.retriever import retrieve_relevant_chunks
from src.llm import (
    build_context,
    build_prompt,
    generate_answer,
    generate_answer_gemini,
    build_summary_prompt,
    build_reduce_summaries_prompt,
    build_section_prompt,
)


def run_rag_pipeline(
    question: str,
    model,
    index,
    metadata: list,
    llm_pipeline: dict,
    generator: str = "local",
    k: int = TOP_K,
    max_distance: float = MAX_DISTANCE,
) -> dict:

    
    retrieved_chunks = retrieve_relevant_chunks(question, model, index, metadata, k=k)

    sources = [
        {"source": c["source"], "chunk_index": c["chunk_index"], "score": c["score"],"text": c["text"]}
        for c in retrieved_chunks
    ]

    if not retrieved_chunks:
        
        return {
            "question": question,
            "sources": sources,
            "answer": "The information was not found in the retrieved documents (no documents have been indexed)."
        }

    best_score = min(c["score"] for c in retrieved_chunks)

    if best_score > max_distance:
        return {
            "question": question,
            "sources": sources,
            "answer": "The information was not found in the retrieved documents (no sufficiently relevant content was located)."
        }

    context = build_context(retrieved_chunks)
    prompt = build_prompt(context, question)

    if generator == "gemini":
        answer = generate_answer_gemini(prompt)
    else:
        answer = generate_answer(prompt, llm_pipeline)    
    


    return {"question": question, "sources": sources, "answer": answer}




# Literature Review Generator

#   Stage 1 (map)    : every document's chunks -> one summary per document
#   Stage 2 (reduce) : all document summaries  -> one corpus digest
#   Stage 3 (generate): the digest             -> each of the 8 sections

REDUCE_CHAR_BUDGET = 1500
LITERATURE_REVIEW_SECTIONS = [
    (
        "Executive Summary",
        "Provide a concise overview of what the uploaded research papers collectively discuss."
    ),
    (
        "Key Findings",
        "Summarize the major findings reported across the papers."
    ),
    (
        "Important Biomarkers / Concepts",
        "List the important biomarkers, concepts, methods or techniques mentioned in the papers."
    ),
    (
        "Limitations",
        "Summarize the limitations or research gaps discussed across the papers."
    ),
    (
        "Conclusion",
        "Provide an overall conclusion based only on the uploaded documents."
    ),
]




def _hierarchical_reduce(
    texts: list,
    llm_pipeline: dict,
    build_prompt_fn,
    max_new_tokens: int = 100,
    char_budget: int = REDUCE_CHAR_BUDGET,
) -> str:
    
    if not texts:
        return ""

    combined = "\n".join(texts)

    if len(combined) <= char_budget:
        prompt = build_prompt_fn(combined)
        return generate_answer_gemini(prompt)
        

    # Too large for one call: batch up to char_budget per batch.
    batches = []
    current_batch = []
    current_len = 0
    for text in texts:
        if current_batch and current_len + len(text) > char_budget:
            batches.append(current_batch)
            current_batch = []
            current_len = 0
        current_batch.append(text)
        current_len += len(text)
    if current_batch:
        batches.append(current_batch)

    
    batch_summaries = [
    generate_answer_gemini(
        build_prompt_fn("\n".join(batch))
    )
    for batch in batches
]

    
    return _hierarchical_reduce(batch_summaries, llm_pipeline, build_reduce_summaries_prompt, max_new_tokens, char_budget)


def _summarize_document(chunks_for_doc: list, llm_pipeline: dict) -> str:
    """
    Produces one summary for a single document from ALL of its chunks
    (not just the first few), via the hierarchical reducer -- so
    documents with many chunks are fully covered rather than truncated.
    """
    chunk_texts = [c["text"] for c in chunks_for_doc]
    return _hierarchical_reduce(chunk_texts, llm_pipeline, build_summary_prompt, max_new_tokens=100)


def generate_literature_review(chunk_metadata: list, llm_pipeline: dict) -> dict:
    
    if not chunk_metadata:
        return {
            title: "No documents have been indexed yet."
            for title, _ in LITERATURE_REVIEW_SECTIONS
        }

    # Group chunks by source document.
    chunks_by_source = {}
    for chunk in chunk_metadata:
        chunks_by_source.setdefault(chunk["source"], []).append(chunk)

    # Stage 1 (map): one summary per document.
    document_summaries = [
        f"[{source}]: {_summarize_document(chunks, llm_pipeline)}"
        for source, chunks in chunks_by_source.items()
    ]

    # Stage 2 (reduce): all document summaries -> one corpus digest.
    digest = _hierarchical_reduce(
        document_summaries, llm_pipeline, build_reduce_summaries_prompt, max_new_tokens=150
    )

    # Stage 3 (generate): each section, independently, from the digest.
    review = {}
    for section_title, section_instruction in LITERATURE_REVIEW_SECTIONS:
        prompt = build_section_prompt(section_title, section_instruction, digest)
        review[section_title] = generate_answer_gemini(prompt)

    return review
    
