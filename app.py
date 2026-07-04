import streamlit as st

from src.config import EMBEDDING_MODEL_NAME, LLM_MODEL_NAME
from src.pdf_loader import extract_text_from_pdfs
from src.chunking import build_chunks
from src.embeddings import load_embedding_model, generate_embeddings
from src.vector_store import build_faiss_index, build_metadata_store, build_bm25_index
from src.llm import load_llm
from src.pipeline import run_rag_pipeline, generate_literature_review
from src.evaluator import (                      # ← new import for dashboard
    compute_metrics,
    build_similarity_bar_chart,
    build_source_pie_chart,
    build_contribution_pie_chart,
    build_timing_chart,
)


st.set_page_config(
    page_title="Document RAG Assistant",
    page_icon="📄",
    layout="centered",
)


# ------------------------------------------------------------------
# Session state
# ------------------------------------------------------------------

def init_session_state():
    if "embedding_model" not in st.session_state:
        with st.spinner("Loading embedding model..."):
            st.session_state.embedding_model = load_embedding_model()

    if "llm_pipeline" not in st.session_state:
        with st.spinner("Loading language model..."):
            st.session_state.llm_pipeline = load_llm()

    defaults = {
        "faiss_index": None,
        "chunk_metadata": None,
        "bm25_index": None,
        "processed_files": [],
        "num_chunks": 0,
        "status": "Not started",
        "last_question": None,
        "last_result": None,
        "literature_review": None,
    }

    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


# ------------------------------------------------------------------
# Backend orchestration
# ------------------------------------------------------------------

def process_uploaded_files(uploaded_files):
    uploaded_dict = {f.name: f.getvalue() for f in uploaded_files}

    with st.status("Processing documents...", expanded=True) as status:
        st.write("📥 Reading and extracting text from PDFs...")
        pdf_results = extract_text_from_pdfs(uploaded_dict)

        st.write("✂️ Creating chunks...")
        chunks = build_chunks(pdf_results)

        st.write("🧮 Generating embeddings...")
        embeddings = generate_embeddings(
            chunks,
            st.session_state.embedding_model,
            show_progress_bar=False,
        )

        st.write("📚 Building FAISS index...")
        faiss_index = build_faiss_index(embeddings)
        chunk_metadata = build_metadata_store(chunks)

        st.write("🔍 Building BM25 index...")
        bm25_index = build_bm25_index(chunks)

        status.update(
            label="✅ Ready for questions",
            state="complete",
            expanded=False,
        )

    return pdf_results, chunks, faiss_index, chunk_metadata, bm25_index


# ------------------------------------------------------------------
# Sidebar
# ------------------------------------------------------------------

def render_sidebar():
    with st.sidebar:
        st.header("Project status")

        st.subheader("Documents uploaded")
        if st.session_state.processed_files:
            for filename in st.session_state.processed_files:
                st.write(f"📄 {filename}")
        else:
            st.caption("No documents processed yet.")

        st.metric("Chunks created", st.session_state.num_chunks)

        st.divider()

        st.subheader("Models in use")
        st.caption("Embedding model")
        st.code(EMBEDDING_MODEL_NAME, language=None)
        st.caption("Available Answer Models")
        st.code(f"Local: {LLM_MODEL_NAME}", language=None)
        st.code("API: Gemini 2.5 Flash", language=None)

        st.divider()

        st.subheader("Processing status")
        if st.session_state.faiss_index is not None:
            st.success(st.session_state.status)
        else:
            st.info(st.session_state.status)


# ------------------------------------------------------------------
# Main page sections
# ------------------------------------------------------------------

def render_upload_section():
    st.header("1. Upload documents")

    uploaded_files = st.file_uploader(
        "Upload one or more PDF files",
        type=["pdf"],
        accept_multiple_files=True,
        help="Upload the research PDFs you'd like to ask questions about.",
    )

    if uploaded_files and st.button("Process documents", type="primary"):
        pdf_results, chunks, faiss_index, chunk_metadata, bm25_index = (
            process_uploaded_files(uploaded_files)
        )

        succeeded = [r["filename"] for r in pdf_results if not r["error"]]
        failed = [(r["filename"], r["error"]) for r in pdf_results if r["error"]]

        st.session_state.faiss_index = faiss_index
        st.session_state.chunk_metadata = chunk_metadata
        st.session_state.bm25_index = bm25_index
        st.session_state.processed_files = succeeded
        st.session_state.num_chunks = len(chunks)
        st.session_state.status = (
            f"{len(succeeded)} document(s) indexed · {len(chunks)} chunks"
        )

        if succeeded:
            st.success(
                f"Processed {len(succeeded)} document(s) into {len(chunks)} chunks. "
                f"You can now ask questions below."
            )

        for filename, error in failed:
            st.warning(f"Could not process **{filename}**: {error}")

    if st.session_state.processed_files:
        st.caption(
            "Currently indexed: " + ", ".join(st.session_state.processed_files)
        )


def render_conversation(question, result):
    st.divider()

    with st.chat_message("user"):
        st.write(question)

    with st.chat_message("assistant"):
        st.write(result["answer"])

    st.subheader("Retrieved sources")

    if not result["sources"]:
        st.caption("No sources were retrieved for this question.")
        return

    for i, s in enumerate(result["sources"], start=1):
        with st.expander(f"Source {i}: {s['source']} · chunk {s['chunk_index']}"):
            col1, col2 = st.columns(2)
            col1.metric("Chunk index", s["chunk_index"])
            col2.metric("Similarity score", f"{s['score']:.4f}")

            st.markdown("**Preview**")
            preview = s["text"][:500]
            if len(s["text"]) > 500:
                preview += "..."
            st.text(preview)

            st.caption("Lower similarity score means a closer match.")


def render_qa_section():
    st.header("2. Ask a question")

    generator = st.radio(
        "Answer Model",
        ["Local FLAN-T5", "Gemini 2.5 Flash"],
        horizontal=True,
    )

    question = st.text_input(
        "Your question",
        placeholder="e.g. What biomarkers are used in breast cancer?",
        label_visibility="collapsed",
    )

    ask_clicked = st.button("Ask", type="primary")

    if ask_clicked:
        if st.session_state.faiss_index is None:
            st.error(
                "Please upload and process at least one PDF before asking a question."
            )
        elif not question.strip():
            st.error("Please enter a question.")
        else:
            with st.spinner("Searching documents and generating answer..."):
                result = run_rag_pipeline(
                    question=question,
                    model=st.session_state.embedding_model,
                    index=st.session_state.faiss_index,
                    metadata=st.session_state.chunk_metadata,
                    llm_pipeline=st.session_state.llm_pipeline,
                    generator="gemini" if generator == "Gemini 2.5 Flash" else "local",
                    bm25_index=st.session_state.bm25_index,
                )

            st.session_state.last_question = question
            st.session_state.last_result = result

    if st.session_state.last_result:
        render_conversation(
            st.session_state.last_question,
            st.session_state.last_result,
        )


def render_literature_review_section():
    st.header("3. Generate Literature Review")

    if st.button("Generate Literature Review", type="primary"):

        if st.session_state.chunk_metadata is None:
            st.error("Please upload and process documents first.")
            return

        with st.spinner("Generating literature review..."):
            review = generate_literature_review(
                chunk_metadata=st.session_state.chunk_metadata,
                llm_pipeline=st.session_state.llm_pipeline,
            )

        st.session_state.literature_review = review

    if st.session_state.literature_review:
        st.subheader("Literature Review")
        for title, content in st.session_state.literature_review.items():
            with st.expander(title, expanded=True):
                st.write(content)


def render_evaluation_dashboard():
    
    st.header("4. Evaluation Dashboard")

    if not st.session_state.last_result:
        st.info("Ask a question in Section 2 to populate the dashboard.")
        return

    result = st.session_state.last_result
    metrics = compute_metrics(
        result,
        num_pdfs=len(st.session_state.processed_files),
        total_chunks_indexed=st.session_state.num_chunks,
    )

    # ------------------------------------------------------------------
    # 4A: Retrieval Metrics
    # ------------------------------------------------------------------
    st.subheader("Retrieval Metrics")

    col1, col2, col3 = st.columns(3)
    col1.metric(
        "Retrieval Time",
        f"{metrics['retrieval_time_ms']:.1f} ms",
    )
    col2.metric(
        "Generation Time",
        f"{metrics['generation_time_ms']:.1f} ms",
    )
    col3.metric(
        "Total Response Time",
        f"{metrics['total_time_ms']:.1f} ms",
    )

    col4, col5, col6 = st.columns(3)
    col4.metric("Uploaded PDFs", metrics["num_pdfs"])
    col5.metric("Indexed Chunks", metrics["total_chunks_indexed"])
    col6.metric("Retrieved Chunks", metrics["num_retrieved_chunks"])

    st.divider()

    # ------------------------------------------------------------------
    # 4B: Answer Information
    # ------------------------------------------------------------------
    st.subheader("Answer Information")

    model_label = (
        "Gemini 2.5 Flash"
        if metrics["generator"] == "gemini"
        else "Local FLAN-T5"
    )

    col7, col8, col9, col10 = st.columns(4)
    col7.metric("Model Used", model_label)
    col8.metric("Answer Length", f"{metrics['answer_length']} chars")
    col9.metric("Prompt Length", f"{metrics['prompt_length']} chars")
    col10.metric("Context Length", f"{metrics['context_length']} chars")

    st.divider()

    # ------------------------------------------------------------------
    # 4C: Visualizations — row 1
    # ------------------------------------------------------------------
    st.subheader("Visualizations")

    col_left, col_right = st.columns(2)

    with col_left:
        st.plotly_chart(
            build_similarity_bar_chart(result["sources"]),
            use_container_width=True,
        )

    with col_right:
        st.plotly_chart(
            build_source_pie_chart(metrics["source_counts"]),
            use_container_width=True,
        )

    # ------------------------------------------------------------------
    # 4D: Visualizations — row 2
    # ------------------------------------------------------------------
    col_left2, col_right2 = st.columns(2)

    with col_left2:
        st.plotly_chart(
            build_contribution_pie_chart(metrics["contribution_counts"]),
            use_container_width=True,
        )

    with col_right2:
        st.plotly_chart(
            build_timing_chart(
                metrics["retrieval_time_ms"],
                metrics["generation_time_ms"],
            ),
            use_container_width=True,
        )


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

def main():
    init_session_state()
    render_sidebar()

    st.title("📄 Document RAG Assistant")
    st.caption("Ask questions grounded in your own uploaded PDF documents.")

    render_upload_section()
    st.divider()
    render_qa_section()
    st.divider()
    render_literature_review_section()
    st.divider()
    render_evaluation_dashboard()    


if __name__ == "__main__":
    main()