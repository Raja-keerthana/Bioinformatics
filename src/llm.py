import os
from dotenv import load_dotenv
from google import genai
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch

from src.config import LLM_MODEL_NAME, MAX_NEW_TOKENS, MAX_INPUT_TOKENS, MAX_CHARS_PER_CHUNK
load_dotenv()

GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")

def load_llm(model_name: str = LLM_MODEL_NAME) -> dict:
    
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
    return {"tokenizer": tokenizer, "model": model}


def build_context(retrieved_chunks: list, max_chars_per_chunk: int = MAX_CHARS_PER_CHUNK) -> str:
   
    context_parts = []
    for chunk in retrieved_chunks:
        label = f"[Source: {chunk['source']}, chunk {chunk['chunk_index']}]"
        trimmed_text = chunk["text"][:max_chars_per_chunk]
        context_parts.append(f"{label}\n{trimmed_text}")
    return "\n\n".join(context_parts)


def build_prompt(context: str, question: str) -> str:
    
    prompt = f"""Answer the question using only the context below. If the answer is not present in the context, clearly state that the information was not found in the retrieved documents.

Question: {question}

Context: {context}"""
    return prompt


def generate_answer(
    prompt: str,
    llm_pipeline: dict,
    max_new_tokens: int = MAX_NEW_TOKENS,
    max_input_tokens: int = MAX_INPUT_TOKENS,
    verbose: bool = False,
) -> str:
   
    try:
        tokenizer = llm_pipeline["tokenizer"]
        model = llm_pipeline["model"]

        if verbose:
            print("  -> tokenizing prompt...")
        inputs = tokenizer(prompt, return_tensors="pt", truncation=True, max_length=max_input_tokens)
        if verbose:
            print("  -> tokenized. input length:", inputs["input_ids"].shape[1])

        if verbose:
            print("  -> generating (may take a moment on CPU)...")
        with torch.no_grad():
            output_ids = model.generate(**inputs, max_new_tokens=max_new_tokens)
        if verbose:
            print("  -> generation finished.")

        answer = tokenizer.decode(output_ids[0], skip_special_tokens=True)
        return answer
    except Exception as e:
        return f"Error while generating answer: {e}"

def generate_answer_gemini(prompt: str) -> str:
    try:
        client = genai.Client(api_key=GEMINI_API_KEY)

        response = client.models.generate_content(
            model="gemini-2.5-flash",
            contents=prompt,
        )

        return response.text.strip()

    except Exception as e:
        return f"Gemini Error: {e}"

def build_summary_prompt(text: str) -> str:
    return f"""Summarize the key points of the following research text in 2-3 sentences. Stay strictly within the information given; do not add outside knowledge.

Text: {text}"""


def build_reduce_summaries_prompt(summaries: str) -> str:
    return f"""Combine the following summaries into a single, more concise summary that preserves the most important point from each one. Do not add information that isn't present in the summaries.

Summaries:
{summaries}"""


def build_section_prompt(
    section_title: str,
    section_instruction: str,
    digest: str,
) -> str:
    return f"""Write the "{section_title}" section of a literature review, using only the information in the digest below.

{section_instruction}

If the digest does not contain relevant information for this section, state that clearly instead of guessing.

Digest of the source documents:
{digest}"""