import json
from pathlib import Path

import faiss
import numpy as np

from src.config import INDEX_DIR

DEFAULT_INDEX_FILENAME = "faiss_index.bin"
DEFAULT_METADATA_FILENAME = "chunk_metadata.json"


def build_faiss_index(embeddings: np.ndarray) -> faiss.Index:
    
    dimension = embeddings.shape[1]
    index = faiss.IndexFlatL2(dimension)
    embeddings = np.ascontiguousarray(embeddings, dtype=np.float32)
    index.add(embeddings)
    return index


def build_metadata_store(chunks: list) -> list:
    
    metadata = []
    for chunk in chunks:
        metadata.append({
            "text": chunk["text"],
            "source": chunk["source"],
            "chunk_index": chunk["chunk_index"]
        })
    return metadata


def save_index(
    index: faiss.Index,
    metadata: list,
    index_dir: Path = INDEX_DIR,
    index_filename: str = DEFAULT_INDEX_FILENAME,
    metadata_filename: str = DEFAULT_METADATA_FILENAME,
) -> None:
    
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)

    index_path = index_dir / index_filename
    metadata_path = index_dir / metadata_filename

    faiss.write_index(index, str(index_path))

    with open(metadata_path, "w") as f:
        json.dump(metadata, f)


def load_index(
    index_dir: Path = INDEX_DIR,
    index_filename: str = DEFAULT_INDEX_FILENAME,
    metadata_filename: str = DEFAULT_METADATA_FILENAME,
) -> tuple:
    
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
