

from src.config import TOP_K, MAX_DISTANCE
from src.retriever import retrieve_relevant_chunks

from src.llm import (
    build_context,
    build_prompt,
    generate_answer,
    generate_answer_gemini,
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
