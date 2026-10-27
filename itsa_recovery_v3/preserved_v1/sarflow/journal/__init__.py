"""Publication styling; legacy portfolio is loaded only when explicitly called."""
from .style import apply_publication_style


def render_portfolio(*args, **kwargs):
    from .portfolio import render_portfolio as render
    return render(*args, **kwargs)
