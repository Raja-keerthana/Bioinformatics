

import numpy as np

from src.config import TOP_K


def retrieve_relevant_chunks(question: str, model, index, metadata: list, k: int = TOP_K) -> list:
    
    
    # Step 1: convert the question into an embedding, using the same
    # model and same vector space as the document chunks.
    query_embedding = model.encode([question], convert_to_numpy=True)
    query_embedding = np.array(query_embedding, dtype=np.float32)

    # Step 2: search FAISS for the k nearest chunk embeddings.
    distances, indices = index.search(query_embedding, k)

    # Step 3: map the returned positions back to real text + metadata.
    results = []
    for dist, idx in zip(distances[0], indices[0]):
        results.append({
            "text": metadata[idx]["text"],
            "source": metadata[idx]["source"],
            "chunk_index": metadata[idx]["chunk_index"],
            "score": dist
        })

    return results
