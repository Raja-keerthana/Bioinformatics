import numpy as np

from src.config import TOP_K, SEMANTIC_WEIGHT, BM25_WEIGHT


CANDIDATE_POOL_SIZE = 20



def _balanced_select(candidates: list, k: int) -> list:
   
    by_source: dict = {}
    for c in candidates:
        by_source.setdefault(c["source"], []).append(c)

    selected = []
    sources = list(by_source.keys())
    cursor = 0  # how many chunks we have already taken from each source

    while len(selected) < k:
        took_any = False
        for source in sources:
            if len(selected) >= k:
                break
            pool = by_source[source]
            if cursor < len(pool):
                selected.append(pool[cursor])
                took_any = True
        cursor += 1
        if not took_any:  # every source's pool is exhausted
            break

    return selected


def _semantic_retrieve(
    question: str,
    model,
    index,
    metadata: list,
    k: int,
) -> list:
   
    # Embed the query in the same vector space as the document chunks.
    query_embedding = model.encode([question], convert_to_numpy=True)
    query_embedding = np.array(query_embedding, dtype=np.float32)

    
    pool_size = max(CANDIDATE_POOL_SIZE, k)
    distances, indices = index.search(query_embedding, pool_size)

    candidates = []
    for dist, idx in zip(distances[0], indices[0]):
        if idx == -1:
            continue
        candidates.append({
            "text": metadata[idx]["text"],
            "source": metadata[idx]["source"],
            "chunk_index": metadata[idx]["chunk_index"],
            "score": float(dist),
        })
    

    selected = _balanced_select(candidates, k)

    
    selected.sort(key=lambda c: c["score"])

    return selected


def _hybrid_retrieve(
    question: str,
    model,
    index,
    metadata: list,
    bm25_index,
    k: int,
) -> list:
    
    pool_size = max(CANDIDATE_POOL_SIZE, k)

    # ------------------------------------------------------------------
    # Step 1: FAISS semantic search
    # ------------------------------------------------------------------
    query_embedding = model.encode([question], convert_to_numpy=True)
    query_embedding = np.array(query_embedding, dtype=np.float32)
    distances, faiss_indices = index.search(query_embedding, pool_size)

    # Map metadata position → raw L2 distance. Filter -1 padding.
    faiss_scores: dict = {}
    for dist, idx in zip(distances[0], faiss_indices[0]):
        if idx != -1:
            faiss_scores[int(idx)] = float(dist)

    # ------------------------------------------------------------------
    # Step 2: BM25 keyword search
    # ------------------------------------------------------------------
    tokenized_query = question.lower().split()
    bm25_raw = bm25_index.get_scores(tokenized_query)  # shape: (num_chunks,)

    
    bm25_top_positions = np.argsort(bm25_raw)[::-1][:pool_size]

    
    bm25_scores: dict = {}
    for idx in bm25_top_positions:
        score = float(bm25_raw[idx])
        if score > 0.0:
            bm25_scores[int(idx)] = score

    # ------------------------------------------------------------------
    # Step 3: Normalise both score types independently to [0, 1]
    #         where 1.0 = best possible match.
    # ------------------------------------------------------------------
    # FAISS: L2 distance → similarity.  0 distance = perfect match → 1.0.
    faiss_max = max(faiss_scores.values(), default=1.0)
    normalised_semantic = {
        idx: 1.0 - (dist / (faiss_max + 1e-9))
        for idx, dist in faiss_scores.items()
    }

    # BM25: already "higher = better"; divide by max to scale to [0, 1].
    bm25_max = max(bm25_scores.values(), default=1.0)
    normalised_bm25 = {
        idx: score / (bm25_max + 1e-9)
        for idx, score in bm25_scores.items()
    }

    # ------------------------------------------------------------------
    # Step 4 & 5: Merge candidate sets and compute combined score.
    #             Deduplication is automatic: the union of keys from
    #             both dicts gives each unique chunk exactly one entry.
    # ------------------------------------------------------------------
    all_positions = set(normalised_semantic.keys()) | set(normalised_bm25.keys())

    if not all_positions:
        return []

    candidates = []
    for idx in all_positions:
        semantic = normalised_semantic.get(idx, 0.0)
        bm25 = normalised_bm25.get(idx, 0.0)
        combined = SEMANTIC_WEIGHT * semantic + BM25_WEIGHT * bm25

        
        inverted = 1.0 - combined

        candidates.append({
            "text": metadata[idx]["text"],
            "source": metadata[idx]["source"],
            "chunk_index": metadata[idx]["chunk_index"],
            "score": inverted,
            
            "_combined": combined,
        })

    
    candidates.sort(key=lambda c: c["_combined"], reverse=True)

    # ------------------------------------------------------------------
    # Step 6: Balanced round-robin selection across source documents.
    # ------------------------------------------------------------------
    selected = _balanced_select(candidates, k)

    # ------------------------------------------------------------------
    # Step 7: Final sort and clean-up.
    # ------------------------------------------------------------------
    
    selected.sort(key=lambda c: c["score"])

   
    for chunk in selected:
        chunk.pop("_combined", None)

    return selected


# ------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------

def retrieve_relevant_chunks(
    question: str,
    model,
    index,
    metadata: list,
    k: int = TOP_K,
    bm25_index=None,
) -> list:
    
    if bm25_index is not None:
        return _hybrid_retrieve(question, model, index, metadata, bm25_index, k)
    return _semantic_retrieve(question, model, index, metadata, k)