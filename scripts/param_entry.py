"""The one constructor for a parameter-table entry (shared by the tables)."""


def entry(key, name, value, unit, kind, category, citation, description,
          source_note="", minimum=None, maximum=None):
    return {"key": key, "name": name, "value": value, "unit": unit,
            "kind": kind, "category": category, "citation": citation,
            "description": description, "source_note": source_note,
            "minimum": minimum, "maximum": maximum}
