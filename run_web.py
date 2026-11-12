"""Start the IP3R Structural Simulator web app on http://127.0.0.1:8194 (local only).

Run: python run_web.py   (needs requirements-web.txt; the desktop app is `python -m ip3r`)
"""
import uvicorn

if __name__ == "__main__":
    uvicorn.run("ip3r.web.app:app", host="127.0.0.1", port=8194)
