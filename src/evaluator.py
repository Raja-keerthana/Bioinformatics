
from collections import Counter

import plotly.express as px
import plotly.graph_objects as go


# ------------------------------------------------------------------
# Metric computation
# ------------------------------------------------------------------

def compute_metrics(
    result: dict,
    num_pdfs: int,
    total_chunks_indexed: int,
) -> dict:
    
    sources = result.get("sources", [])
    answer = result.get("answer", "")
    meta = result.get("eval_meta", {})

    source_counts = Counter(s["source"] for s in sources)
    contribution_counts = Counter(
        s.get("retrieval_source", "semantic") for s in sources
    )

    return {
        # Timing (ms)
        "retrieval_time_ms": meta.get("retrieval_time_ms", 0.0),
        "generation_time_ms": meta.get("generation_time_ms", 0.0),
        "total_time_ms": meta.get("total_time_ms", 0.0),

        # Index / retrieval counts
        "num_pdfs": num_pdfs,
        "total_chunks_indexed": total_chunks_indexed,
        "num_retrieved_chunks": len(sources),

        # Answer information
        "generator": meta.get("generator", "local"),
        "prompt_length": meta.get("prompt_length", 0),
        "context_length": meta.get("context_length", 0),
        "answer_length": len(answer),

        # Pre-aggregated for chart builders
        "sources": sources,
        "source_counts": dict(source_counts),
        "contribution_counts": dict(contribution_counts),
    }


# ------------------------------------------------------------------
# Plotly figure builders
# ------------------------------------------------------------------

_CONTRIBUTION_COLOR_MAP = {
    "Semantic (FAISS)": "#4C78A8",
    "Keyword (BM25)":   "#F58518",
    "Both":             "#54A24B",
}

_SOURCE_COLOR_MAP = {
    "semantic": "#4C78A8",
    "keyword":  "#F58518",
    "both":     "#54A24B",
}


def build_similarity_bar_chart(sources: list) -> go.Figure:
    
    if not sources:
        fig = go.Figure()
        fig.update_layout(title="No chunks retrieved", height=380)
        return fig

    labels = [
        f"{s['source']}\nchunk {s['chunk_index']}" for s in sources
    ]
    scores = [round(s["score"], 4) for s in sources]
    colors = [
        _SOURCE_COLOR_MAP.get(s.get("retrieval_source", "semantic"), "#4C78A8")
        for s in sources
    ]

    fig = go.Figure(
        go.Bar(
            x=labels,
            y=scores,
            marker_color=colors,
            text=[str(sc) for sc in scores],
            textposition="outside",
            hovertemplate="<b>%{x}</b><br>Score: %{y:.4f}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Similarity Score per Retrieved Chunk",
        xaxis_title="Chunk",
        yaxis_title="Score  (lower = closer match)",
        height=380,
        margin=dict(t=55, b=10, l=10, r=10),
        showlegend=False,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


def build_source_pie_chart(source_counts: dict) -> go.Figure:
    """
    Pie chart showing how retrieved chunks are distributed across the
    uploaded PDF files.
    """
    if not source_counts:
        fig = go.Figure()
        fig.update_layout(title="No data", height=380)
        return fig

    fig = px.pie(
        names=list(source_counts.keys()),
        values=list(source_counts.values()),
        title="Retrieved Chunks by Source Document",
        color_discrete_sequence=px.colors.qualitative.Set2,
        hole=0.3,
    )
    fig.update_layout(
        height=380,
        margin=dict(t=55, b=10, l=10, r=10),
        paper_bgcolor="rgba(0,0,0,0)",
    )
    fig.update_traces(textposition="inside", textinfo="percent+label")
    return fig


def build_contribution_pie_chart(contribution_counts: dict) -> go.Figure:
    
    if not contribution_counts:
        fig = go.Figure()
        fig.update_layout(title="No data", height=380)
        return fig

    label_map = {
        "semantic": "Semantic (FAISS)",
        "keyword":  "Keyword (BM25)",
        "both":     "Both",
    }
    labels = [label_map.get(k, k) for k in contribution_counts.keys()]
    values = list(contribution_counts.values())
    colors = [_CONTRIBUTION_COLOR_MAP.get(l, "#888888") for l in labels]

    fig = px.pie(
        names=labels,
        values=values,
        title="Retrieval Contribution: Semantic vs Keyword",
        hole=0.3,
    )
    fig.update_traces(
        marker_colors=colors,
        textposition="inside",
        textinfo="percent+label",
    )
    fig.update_layout(
        height=380,
        margin=dict(t=55, b=10, l=10, r=10),
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


def build_timing_chart(
    retrieval_time_ms: float,
    generation_time_ms: float,
) -> go.Figure:
    
    fig = go.Figure()

    fig.add_trace(go.Bar(
        name="Retrieval",
        x=[retrieval_time_ms],
        y=["Pipeline"],
        orientation="h",
        marker_color="#4C78A8",
        text=[f"Retrieval  {retrieval_time_ms:.1f} ms"],
        textposition="inside",
        insidetextanchor="middle",
        hovertemplate="Retrieval: %{x:.1f} ms<extra></extra>",
    ))

    fig.add_trace(go.Bar(
        name="Generation",
        x=[generation_time_ms],
        y=["Pipeline"],
        orientation="h",
        marker_color="#F58518",
        text=[f"Generation  {generation_time_ms:.1f} ms"],
        textposition="inside",
        insidetextanchor="middle",
        hovertemplate="Generation: %{x:.1f} ms<extra></extra>",
    ))

    total = retrieval_time_ms + generation_time_ms
    fig.update_layout(
        title=f"Processing Time Breakdown  (total {total:.1f} ms)",
        barmode="stack",
        xaxis_title="Time (ms)",
        height=230,
        margin=dict(t=55, b=10, l=10, r=10),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.08,
            xanchor="right",
            x=1,
        ),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig
