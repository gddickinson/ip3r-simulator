"""A small Newick reader and the clade questions the origin checks ask.

Independent of ``ip3r_genes/scripts/s7_lib.py`` on purpose: the sister-pair
check parses the committed tree with this code and asks its own question of
it. Internal node labels (IQ-TREE's ``aLRT/UFBoot``) are kept as ``label``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Node", "parse", "leaves", "mrca", "smallest_clade_containing"]


@dataclass
class Node:
    label: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = field(default=None, repr=False)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def _unquote(s: str) -> str:
    s = s.strip()
    return s[1:-1] if len(s) > 1 and s[0] == s[-1] == "'" else s


def parse(text: str) -> Node:
    """Parse one Newick tree (terminated by ``;``)."""
    s = text.strip()
    if not s.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    pos = 0

    def read_label() -> str:
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            end = s.index("'", pos + 1)
            lab = s[pos:end + 1]
            pos = end + 1
            return _unquote(lab)
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        return s[start:pos].strip()

    def read_node(parent: Node | None) -> Node:
        nonlocal pos
        node = Node(parent=parent)
        if s[pos] == "(":
            pos += 1
            while True:
                node.children.append(read_node(node))
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"unexpected {s[pos]!r} at {pos}")
        node.label = read_label()
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",);":
                pos += 1
            node.length = float(s[start:pos])
        return node

    root = read_node(None)
    if s[pos] != ";":
        raise ValueError(f"trailing text at {pos}")
    return root


def leaves(node: Node) -> list[Node]:
    out, stack = [], [node]
    while stack:
        n = stack.pop()
        if n.is_leaf:
            out.append(n)
        else:
            stack.extend(n.children)
    return out


def _ancestors(n: Node) -> list[Node]:
    out = []
    while n is not None:
        out.append(n)
        n = n.parent
    return out


def mrca(root: Node, labels: set[str]) -> Node:
    """Most recent common ancestor of the leaves carrying ``labels``."""
    tips = [lf for lf in leaves(root) if lf.label in labels]
    missing = labels - {t.label for t in tips}
    if missing:
        raise KeyError(f"{len(missing)} labels not in tree, e.g. {sorted(missing)[:2]}")
    common = _ancestors(tips[0])
    for t in tips[1:]:
        anc = set(map(id, _ancestors(t)))
        common = [a for a in common if id(a) in anc]
    return common[0]


def smallest_clade_containing(root: Node, labels: set[str]) -> set[str]:
    """All leaf labels under the MRCA of ``labels`` (rooted reading)."""
    return {lf.label for lf in leaves(mrca(root, labels))}
