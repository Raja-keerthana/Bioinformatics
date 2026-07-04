import numpy as np

from src.config import TOP_K, SEMANTIC_WEIGHT, BM25_WEIGHT

CANDIDATE_POOL_SIZE = 20


# ------------------------------------------------------------------
# Private helpers
# ------------------------------------------------------------------

def _balanced_select(candidates: list, k: int) -> list:
   
    by_source: dict = {}
    for c in candidates:
        by_source.setdefault(c["source"], []).append(c)

    selected = []
    sources = list(by_source.keys())
    cursor = 0

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
        if not took_any:
            break

    return selected


def _semantic_retrieve(
    question: str,
    model,
    index,
    metadata: list,
    k: int,
) -> list:
   
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
            "retrieval_source": "semantic",
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

    # Step 1: FAISS
    query_embedding = model.encode([question], convert_to_numpy=True)
    query_embedding = np.array(query_embedding, dtype=np.float32)
    distances, faiss_indices = index.search(query_embedding, pool_size)

    faiss_scores: dict = {}
    for dist, idx in zip(distances[0], faiss_indices[0]):
        if idx != -1:
            faiss_scores[int(idx)] = float(dist)

    # Step 2: BM25
    tokenized_query = question.lower().split()
    bm25_raw = bm25_index.get_scores(tokenized_query)
    bm25_top_positions = np.argsort(bm25_raw)[::-1][:pool_size]

    bm25_scores: dict = {}
    for idx in bm25_top_positions:
        score = float(bm25_raw[idx])
        if score > 0.0:
            bm25_scores[int(idx)] = score

    # Step 3: Normalise
    faiss_max = max(faiss_scores.values(), default=1.0)
    normalised_semantic = {
        idx: 1.0 - (dist / (faiss_max + 1e-9))
        for idx, dist in faiss_scores.items()
    }

    bm25_max = max(bm25_scores.values(), default=1.0)
    normalised_bm25 = {
        idx: score / (bm25_max + 1e-9)
        for idx, score in bm25_scores.items()
    }

    # Step 4 & 5: Merge, assign retrieval_source per chunk
    all_positions = set(normalised_semantic.keys()) | set(normalised_bm25.keys())

    if not all_positions:
        return []

    candidates = []
    for idx in all_positions:
        in_faiss = idx in normalised_semantic
        in_bm25 = idx in normalised_bm25

        # Determine which search(es) found this chunk
        if in_faiss and in_bm25:
            retrieval_source = "both"
        elif in_faiss:
            retrieval_source = "semantic"
        else:
            retrieval_source = "keyword"

        semantic = normalised_semantic.get(idx, 0.0)
        bm25 = normalised_bm25.get(idx, 0.0)
        combined = SEMANTIC_WEIGHT * semantic + BM25_WEIGHT * bm25

        candidates.append({
            "text": metadata[idx]["text"],
            "source": metadata[idx]["source"],
            "chunk_index": metadata[idx]["chunk_index"],
            "score": 1.0 - combined,      # inverted: lower = better
            "retrieval_source": retrieval_source,
            "_combined": combined,          # internal sort key only
        })

    # Sort best-first so _balanced_select picks the best from each source
    candidates.sort(key=lambda c: c["_combined"], reverse=True)

    # Step 6: Balanced selection
    selected = _balanced_select(candidates, k)

    # Step 7: Final sort and remove internal key
    selected.sort(key=lambda c: c["score"])
    for chunk in selected:
        chunk.pop("_combined", None)

    return selected



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