"""
LLM loading, prompt construction, and answer generation.

Moved from Section 6 of the original Colab notebook. Covers everything
needed to turn a question plus retrieved context into a generated
answer: loading the model, building the context block, filling the
prompt template, and generating + decoding the response.

The orchestration that ties retrieval and generation together
(run_rag_pipeline, including the relevance gate from "Hallucination
Control") is deliberately NOT here — it belongs in pipeline.py, built
in the next step.

This module has no dependency on Colab or any UI layer.
"""

from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch

from src.config import LLM_MODEL_NAME, MAX_NEW_TOKENS, MAX_INPUT_TOKENS, MAX_CHARS_PER_CHUNK


def load_llm(model_name: str = LLM_MODEL_NAME) -> dict:
    """
    Loads a tokenizer + text-to-text generation model directly via
    AutoTokenizer / AutoModelForSeq2SeqLM. This bypasses transformers'
    pipeline() task-registry lookup, which proved unreliable across
    transformers versions during development (see the
    "Unknown task text2text-generation" issue from earlier).

    Parameters
    ----------
    model_name : str
        Defaults to src/config.py's LLM_MODEL_NAME ("google/flan-t5-base").

    Returns
    -------
    dict: {"tokenizer": ..., "model": ...}
    """
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
    return {"tokenizer": tokenizer, "model": model}


def build_context(retrieved_chunks: list, max_chars_per_chunk: int = MAX_CHARS_PER_CHUNK) -> str:
    """
    Combines retrieved chunks into a single labeled context block.

    Parameters
    ----------
    retrieved_chunks : list
        The list of dicts from retriever.retrieve_relevant_chunks.
        Each dict has 'text', 'source', 'chunk_index', 'score'.
    max_chars_per_chunk : int
        Each chunk's text is trimmed to this length so more chunks fit
        inside the LLM's fixed input token budget. Defaults to
        src/config.py's MAX_CHARS_PER_CHUNK.

    Returns
    -------
    A single string: one labeled passage per chunk, separated by blank
    lines.
    """
    context_parts = []
    for chunk in retrieved_chunks:
        label = f"[Source: {chunk['source']}, chunk {chunk['chunk_index']}]"
        trimmed_text = chunk["text"][:max_chars_per_chunk]
        context_parts.append(f"{label}\n{trimmed_text}")
    return "\n\n".join(context_parts)


def build_prompt(context: str, question: str) -> str:
    """
    Fills the RAG prompt template.

    Instructions and the question come FIRST, context LAST. This
    ordering matters: truncation (in generate_answer) cuts content
    from the end when the combined text exceeds the model's input
    limit — putting context last means only context gets trimmed if
    needed, never the instructions or the question itself.

    Parameters
    ----------
    context : str
        The string returned by build_context.
    question : str
        The user's natural language question.

    Returns
    -------
    The complete prompt string, ready to tokenize.
    """
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
    """
    Tokenizes the prompt, generates a response, and decodes it back
    into text.

    Parameters
    ----------
    prompt : str
        The string returned by build_prompt.
    llm_pipeline : dict
        The dict returned by load_llm: {"tokenizer": ..., "model": ...}.
    max_new_tokens : int
        Cap on generated answer length. Defaults to src/config.py's
        MAX_NEW_TOKENS.
    max_input_tokens : int
        The model's hard input length limit; longer prompts are
        truncated from the end. Defaults to src/config.py's
        MAX_INPUT_TOKENS.
    verbose : bool
        If True, prints progress at each stage (tokenizing,
        generating, done) — useful when debugging a slow or stuck
        generation call. Defaults to False: a library function
        shouldn't print by default, unlike the debug version used
        while troubleshooting earlier.

    Returns
    -------
    The generated answer string, or an "Error while generating
    answer: ..." string if generation failed for any reason — this
    function never raises, so one bad call can't crash a calling loop.
    """
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
