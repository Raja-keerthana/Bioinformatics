

from langchain_text_splitters import RecursiveCharacterTextSplitter

from src.config import CHUNK_SIZE, CHUNK_OVERLAP


def build_chunks(pdf_results: list, chunk_size: int = CHUNK_SIZE, chunk_overlap: int = CHUNK_OVERLAP) -> list:
    
    splitter = RecursiveCharacterTextSplitter(
        chunk_size=chunk_size,
        chunk_overlap=chunk_overlap,
        separators=["\n\n", "\n", ". ", " ", ""]
    )

    all_chunks = []

    for result in pdf_results:
        if result["error"]:
            continue  # skip PDFs that failed to load

        text_pieces = splitter.split_text(result["text"])

        for i, piece in enumerate(text_pieces):
            all_chunks.append({
                "text": piece,
                "source": result["filename"],
                "chunk_index": i
            })

    return all_chunks
