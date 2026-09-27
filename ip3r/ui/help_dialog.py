"""The non-modal guide window (Help → Guide, F1; Round 7.18, after PIEZO1's).

A topic list on the left and the topic on the right; the last entry is the
shortcut table, built from :data:`~ip3r.ui.help_content.SHORTCUTS`. Shipped
documents open with the system's viewer (:func:`open_document`).
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import QDialog, QHBoxLayout, QListWidget, QTextBrowser

from ..config import PROJECT_ROOT
from .help_content import SHORTCUTS, TOPICS

__all__ = ["HelpDialog", "open_document", "SHORTCUT_TOPIC"]

SHORTCUT_TOPIC = "Keyboard shortcuts"


def shortcut_html() -> str:
    rows = "".join(f"<tr><td><b>{k}</b>&nbsp;&nbsp;</td><td>{v}</td></tr>"
                   for k, v in SHORTCUTS)
    return f"<table>{rows}</table>"


def open_document(relative: str) -> bool:
    """Open a shipped document; False when it is not there."""
    path = PROJECT_ROOT / relative
    return path.exists() and QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))


class HelpDialog(QDialog):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("IP3R Structural Simulator — guide")
        self.setModal(False)
        self.resize(820, 560)
        lay = QHBoxLayout(self)
        self.topics = QListWidget()
        self.topics.addItems([*TOPICS, SHORTCUT_TOPIC])
        self.topics.setMaximumWidth(230)
        self.body = QTextBrowser()
        self.body.setOpenExternalLinks(True)
        lay.addWidget(self.topics)
        lay.addWidget(self.body, 1)
        self.topics.currentTextChanged.connect(self.show_topic)
        self.topics.setCurrentRow(0)

    def show_topic(self, name: str) -> None:
        html = shortcut_html() if name == SHORTCUT_TOPIC else TOPICS.get(name, "")
        self.body.setHtml(f"<h2>{name}</h2>{html}")
        items = self.topics.findItems(name, Qt.MatchFlag.MatchExactly)
        if items and self.topics.currentItem() is not items[0]:
            self.topics.setCurrentItem(items[0])
