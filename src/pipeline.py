from src.config import TOP_K, MAX_DISTANCE, HYBRID_MAX_DISTANCE
from src.retriever import retrieve_relevant_chunks
from src.llm import (
    build_context,
    build_prompt,
    generate_answer,
    generate_answer_gemini,
    build_summary_prompt,
    build_final_review_prompt,
)


# ------------------------------------------------------------------
# Question Answering
# ------------------------------------------------------------------

def run_rag_pipeline(
    question: str,
    model,
    index,
    metadata: list,
    llm_pipeline: dict,
    generator: str = "local",
    k: int = TOP_K,
    max_distance: float = MAX_DISTANCE,
    bm25_index=None,
) -> dict:
    
    effective_max_distance = (
        HYBRID_MAX_DISTANCE if bm25_index is not None else max_distance
    )

    retrieved_chunks = retrieve_relevant_chunks(
        question, model, index, metadata, k=k, bm25_index=bm25_index,
    )

    sources = [
        {
            "source": c["source"],
            "chunk_index": c["chunk_index"],
            "score": c["score"],
            "text": c["text"],
        }
        for c in retrieved_chunks
    ]

    if not retrieved_chunks:
        return {
            "question": question,
            "sources": sources,
            "answer": "The information was not found in the retrieved documents.",
        }

    best_score = min(c["score"] for c in retrieved_chunks)

    if best_score > effective_max_distance:
        return {
            "question": question,
            "sources": sources,
            "answer": "The information was not found in the retrieved documents.",
        }

    context = build_context(retrieved_chunks)
    prompt = build_prompt(context, question)

    if generator == "gemini":
        answer = generate_answer_gemini(prompt)
    else:
        answer = generate_answer(prompt, llm_pipeline)

    return {"question": question, "sources": sources, "answer": answer}


# ------------------------------------------------------------------
# Literature Review 
# ------------------------------------------------------------------

_MAX_CHARS_PER_DOCUMENT = 30_000

_SECTION_KEYS = [
    "Executive Summary",
    "Key Findings",
    "Important Biomarkers / Concepts",
    "Limitations",
    "Conclusion",
]


def _prepare_document_text(chunks: list, max_chars: int = _MAX_CHARS_PER_DOCUMENT) -> str:
    
    sampled = chunks[::2] if len(chunks) > 6 else chunks
    text = "\n\n".join(chunk["text"] for chunk in sampled)
    return text[:max_chars]


def _normalize_heading(line: str) -> str:
    
    text = line.strip()
    text = text.lstrip("#").strip()
    text = text.strip("*").strip()
    text = text.rstrip(":").strip()
    return text


def _parse_sections(review_text: str, section_keys: list) -> dict:
   
    sections = {k: "" for k in section_keys}
    current_key = None
    buffer = []

    for line in review_text.splitlines():
        normalized = _normalize_heading(line)
        matched = next(
            (k for k in section_keys if k.lower() == normalized.lower()),
            None,
        )
        if matched:
            if current_key is not None:
                sections[current_key] = "\n".join(buffer).strip()
            current_key = matched
            buffer = []
        elif current_key is not None:
            buffer.append(line)

    if current_key is not None:
        sections[current_key] = "\n".join(buffer).strip()

    if all(not v.strip() for v in sections.values()):
        sections[section_keys[0]] = review_text.strip()

    return sections


def _summarize_document(source_name: str, chunks: list) -> str:
    
    text = _prepare_document_text(chunks)
    prompt = build_summary_prompt(text)
    return generate_answer_gemini(prompt)


def generate_literature_review(
    chunk_metadata: list,
    llm_pipeline: dict,
) -> dict:
    
    if not chunk_metadata:
        return {k: "No documents have been uploaded yet." for k in _SECTION_KEYS}

    chunks_by_source: dict = {}
    for chunk in chunk_metadata:
        chunks_by_source.setdefault(chunk["source"], []).append(chunk)

    document_summaries = [
        f"Document: {source}\n{_summarize_document(source, chunks)}"
        for source, chunks in chunks_by_source.items()
    ]

    review_text = generate_answer_gemini(
        build_final_review_prompt("\n\n".join(document_summaries))
    )

    return _parse_sections(review_text, _SECTION_KEYS)