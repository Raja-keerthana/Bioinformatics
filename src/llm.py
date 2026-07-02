import os

from dotenv import load_dotenv
from google import genai
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch

from src.config import (
    LLM_MODEL_NAME,
    MAX_NEW_TOKENS,
    MAX_INPUT_TOKENS,
    MAX_CHARS_PER_CHUNK,
)

load_dotenv()

GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
client = genai.Client(api_key=GEMINI_API_KEY)


# ------------------------------------------------------------------
# Local FLAN-T5
# ------------------------------------------------------------------

def load_llm(model_name: str = LLM_MODEL_NAME) -> dict:
    """Loads the local FLAN-T5 tokenizer and model."""
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
    return {"tokenizer": tokenizer, "model": model}


def generate_answer(
    prompt: str,
    llm_pipeline: dict,
    max_new_tokens: int = MAX_NEW_TOKENS,
    max_input_tokens: int = MAX_INPUT_TOKENS,
) -> str:
    """Generates an answer using the local FLAN-T5 model."""
    tokenizer = llm_pipeline["tokenizer"]
    model = llm_pipeline["model"]

    inputs = tokenizer(
        prompt,
        return_tensors="pt",
        truncation=True,
        max_length=max_input_tokens,
    )

    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
        )

    return tokenizer.decode(outputs[0], skip_special_tokens=True)


# ------------------------------------------------------------------
# Gemini
# ------------------------------------------------------------------

def generate_answer_gemini(prompt: str) -> str:
    """Generates an answer using the Gemini API."""
    try:
        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=prompt,
        )
        return response.text.strip()
    except Exception as e:
        return f"Gemini Error: {e}"


# ------------------------------------------------------------------
# QA prompt builders 
# ------------------------------------------------------------------

def build_context(
    retrieved_chunks: list,
    max_chars_per_chunk: int = MAX_CHARS_PER_CHUNK,
) -> str:
    """Builds a labeled context block from retrieved chunks."""
    context_parts = []
    for chunk in retrieved_chunks:
        label = f"[Source: {chunk['source']}, chunk {chunk['chunk_index']}]"
        trimmed = chunk["text"][:max_chars_per_chunk]
        context_parts.append(f"{label}\n{trimmed}")
    return "\n\n".join(context_parts)


def build_prompt(context: str, question: str) -> str:
    """Builds the QA prompt for local or Gemini generation."""
    return f"""
Answer the question using ONLY the context below.

If the answer is not present, clearly say it was not found.

Question:
{question}

Context:
{context}
"""


# ------------------------------------------------------------------
# Literature Review prompt builders
# ------------------------------------------------------------------

def build_summary_prompt(text: str) -> str:
    """
    Prompt for summarizing one document's representative text.
    Called exactly once per uploaded document.
    """
    return f"""Summarize the following research paper extract.

Write 5-8 concise bullet points covering:
- Main research objective
- Key findings and results
- Methods and techniques used
- Important biomarkers or concepts mentioned
- Limitations, if any

Use ONLY the provided text. Do not add outside information.

Paper:

{text}"""


def build_final_review_prompt(document_summaries: str) -> str:
    
    return f"""You are writing a structured academic literature review.

Use ONLY the document summaries provided below. Do not add any information not present in those summaries.

Write your response using EXACTLY these five section headings. Each heading must appear on its own line, preceded by ## with no other text on that line:

## Executive Summary
## Key Findings
## Important Biomarkers / Concepts
## Limitations
## Conclusion

Write 3-5 sentences per section. If a section cannot be answered from the summaries, write: "Not enough information available."

Document Summaries:

{document_summaries}"""