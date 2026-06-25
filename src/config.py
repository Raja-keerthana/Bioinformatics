from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

# Paths
DATA_DIR = BASE_DIR / "data"
UPLOADS_DIR = DATA_DIR / "uploads"
INDEX_DIR = DATA_DIR / "index"

# Chunking
CHUNK_SIZE = 1000
CHUNK_OVERLAP = 200

# Embeddings

EMBEDDING_MODEL_NAME = "BAAI/bge-base-en-v1.5"

# Retrieval
TOP_K = 5

# Generation
LLM_MODEL_NAME = "google/flan-t5-small"

MAX_NEW_TOKENS = 200
MAX_INPUT_TOKENS = 512
MAX_CHARS_PER_CHUNK = 600

# Relevance filtering
MAX_DISTANCE = 1.0