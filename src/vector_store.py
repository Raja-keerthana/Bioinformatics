import json
from pathlib import Path

import faiss
import numpy as np
from rank_bm25 import BM25Okapi

from src.config import INDEX_DIR

DEFAULT_INDEX_FILENAME = "faiss_index.bin"
DEFAULT_METADATA_FILENAME = "chunk_metadata.json"


def build_faiss_index(embeddings: np.ndarray) -> faiss.Index:
    """Build a FAISS index from embeddings."""
    dimension = embeddings.shape[1]
    index = faiss.IndexFlatL2(dimension)

    embeddings = np.ascontiguousarray(embeddings, dtype=np.float32)
    index.add(embeddings)

    return index


def build_metadata_store(chunks: list) -> list:
    """Store chunk information aligned with the FAISS index."""
    return [
        {
            "text": chunk["text"],
            "source": chunk["source"],
            "chunk_index": chunk["chunk_index"],
        }
        for chunk in chunks
    ]


def build_bm25_index(chunks: list) -> BM25Okapi:
    """Build a BM25 index for keyword search."""
    tokenized_corpus = [
        chunk["text"].lower().split()
        for chunk in chunks
    ]
    return BM25Okapi(tokenized_corpus)


def save_index(
    index: faiss.Index,
    metadata: list,
    index_dir: Path = INDEX_DIR,
    index_filename: str = DEFAULT_INDEX_FILENAME,
    metadata_filename: str = DEFAULT_METADATA_FILENAME,
) -> None:
    """Save the FAISS index and metadata to disk."""
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)

    faiss.write_index(index, str(index_dir / index_filename))

    with open(index_dir / metadata_filename, "w") as f:
        json.dump(metadata, f)


def load_index(
    index_dir: Path = INDEX_DIR,
    index_filename: str = DEFAULT_INDEX_FILENAME,
    metadata_filename: str = DEFAULT_METADATA_FILENAME,
) -> tuple:
    """Load a saved FAISS index and metadata."""
    index_dir = Path(index_dir)

    index_path = index_dir / index_filename
    metadata_path = index_dir / metadata_filename

    if not index_path.exists():
        raise FileNotFoundError(f"No saved FAISS index found at {index_path}")

    if not metadata_path.exists():
        raise FileNotFoundError(f"No saved metadata found at {metadata_path}")

    index = faiss.read_index(str(index_path))

    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    return index, metadata