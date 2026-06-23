import io

from pypdf import PdfReader


def extract_text_from_pdf(file_bytes: bytes, filename: str) -> dict:
    
    try:
        pdf_stream = io.BytesIO(file_bytes)
        reader = PdfReader(pdf_stream)
        num_pages = len(reader.pages)

        full_text = ""
        for page in reader.pages:
            page_text = page.extract_text()
            if page_text:  # some pages (e.g. scanned images) return None
                full_text += page_text + "\n"

        if not full_text.strip():
            return {
                "filename": filename,
                "num_pages": num_pages,
                "text": "",
                "error": "No extractable text found (this PDF may be scanned images, not real text)."
            }

        return {
            "filename": filename,
            "num_pages": num_pages,
            "text": full_text,
            "error": None
        }

    except Exception as e:
        return {
            "filename": filename,
            "num_pages": 0,
            "text": "",
            "error": f"Could not read this PDF: {e}"
        }


def extract_text_from_pdfs(uploaded_files: dict) -> list:
    
    results = []
    for filename, file_bytes in uploaded_files.items():
        result = extract_text_from_pdf(file_bytes, filename)
        results.append(result)
    return results
