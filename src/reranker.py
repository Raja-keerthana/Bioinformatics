from src.config import RERANKER_MODEL_NAME
from sentence_transformers.cross_encoder import CrossEncoder


def load_reranker(model_name: str = RERANKER_MODEL_NAME) -> CrossEncoder:
    
    return CrossEncoder(model_name)


def rerank_chunks(
    question: str,
    chunks: list,
    reranker: CrossEncoder,
    top_k: int,
) -> list:
    
    if not chunks:
        return chunks

    # Build (question, passage) pairs — the CrossEncoder input format.
    pairs = [(question, chunk["text"]) for chunk in chunks]

   
    ce_scores = reranker.predict(pairs)

    # Zip scores with chunks, sort descending, return top_k.
    ranked = sorted(
        zip(ce_scores, chunks),
        key=lambda x: x[0],
        reverse=True,
    )

    return [chunk for _, chunk in ranked[:top_k]]
