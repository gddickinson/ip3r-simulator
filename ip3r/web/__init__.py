"""Web interface: FastAPI app + static front end (also exported to GitHub Pages).

Contents: ``app`` (routes under /api, ``STATIC_ENDPOINTS``), ``codec`` (int16/base64
coordinate encoding), ``export_structures`` and ``export_dynamics`` (precompute the
viewer files, the ANM, the 8TKG -> 8TKF transition and the state panel with the
repo's own code), ``gating_api`` (gating constants + reference bells),
``precomputed/`` (their committed JSON) and ``static/`` (the ES-module front end).
"""
