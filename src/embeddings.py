"""
Embedding generation.

Moved from Section 3 of the original Colab notebook. Loads a
sentence-transformers model and converts chunk text into numeric
vectors (embeddings) that FAISS can later index and search.

This module has no dependency on Colab or any UI layer — it takes the
`chunks` list produced by chunking.py and returns a NumPy array of
embeddings, nothing else.
"""

from sentence_transformers import SentenceTransformer
import numpy as np

from src.config import EMBEDDING_MODEL_NAME


def load_embedding_model(model_name: str = EMBEDDING_MODEL_NAME) -> SentenceTransformer:
    """
    Loads a pretrained sentence embedding model from sentence-transformers.
    Downloads the model weights on first use; cached locally afterward.

    Parameters
    ----------
    model_name : str
        Defaults to src/config.py's EMBEDDING_MODEL_NAME ("all-MiniLM-L6-v2").

    Returns
    -------
    A loaded SentenceTransformer model, ready for .encode() calls.
    """
    model = SentenceTransformer(model_name)
    return model


def generate_embeddings(chunks: list, model: SentenceTransformer, show_progress_bar: bool = True) -> np.ndarray:
    """
    Converts every chunk's text into a numeric vector (embedding).

    Parameters
    ----------
    chunks : list
        The list of dicts produced by chunking.build_chunks. Each dict
        has 'text', 'source', 'chunk_index'.
    model : SentenceTransformer
        A model returned by load_embedding_model.
    show_progress_bar : bool
        Whether sentence-transformers prints its own progress bar
        during encoding. Defaults to True to match the notebook's
        original behavior; can be set to False in non-interactive
        environments (e.g. inside a Streamlit app).

    Returns
    -------
    NumPy array of shape (num_chunks, embedding_dimension). Row i
    corresponds to chunks[i] — this alignment is what later lets FAISS
    map a matched vector back to its original chunk.
    """
    texts = [chunk["text"] for chunk in chunks]
    embeddings = model.encode(
        texts,
        show_progress_bar=show_progress_bar,
        convert_to_numpy=True
    )
    return embeddings
