#!/usr/bin/env python3
"""Assemble docs/manuscript/MANUSCRIPT_DRAFT.md from the section drafts (no text is written here, only moved).

Order: title, abstract, introduction, results 1-15, discussion, methods (evaluation framework, model, end
analyses), figure legends (1-6, supplementary), references.  Author notes (Chinese) are left out.
"""
import os
import re

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "docs", "manuscript")


def read(name):
    return open(os.path.join(D, name)).read()


def block(text, heading):
    """Body under `heading` (exact line) up to the next heading of the same or higher level."""
    lines = text.split("\n")
    i = lines.index(heading)
    lvl = len(heading) - len(heading.lstrip("#"))
    out = []
    for line in lines[i + 1:]:
        m = re.match(r"^(#+) ", line)
        if m and len(m.group(1)) <= lvl:
            break
        out.append(line)
    while out and out[-1].strip() in ("", "---"):
        out.pop()
    return "\n".join(out).strip("\n")


def demote(text, by=1):
    return re.sub(r"^(#+) ", lambda m: "#" * (len(m.group(1)) + by) + " ", text, flags=re.M)


def main():
    ab, intro, bm, lo, en = (read(f) for f in ("DRAFT_abstract.md", "DRAFT_intro_discussion.md", "DRAFT_results_benchmark.md",
                                                "DRAFT_results_locus.md", "DRAFT_results_ends.md"))
    leg = read("figures/LEGENDS.md")
    title = " ".join(block(intro, "## Working title (proposed, replaces the October title)").split("\n\n")[0].strip("*\n ").split())
    legends = {}
    for m in re.finditer(r"^\*\*((?:Supplementary )?Figure (\w+)) \|.*?(?=^\*\*(?:Supplementary )?Figure |^# |\Z)", leg, flags=re.M | re.S):
        legends[m.group(1)] = m.group(0).strip()
    order = [f"Figure {i}" for i in range(1, 7)] + ["Supplementary Figure 1"]
    refs = block(en, "## References cited in sections 12–15") + "\n" + block(intro, "## References to add (beyond the results drafts)")
    parts = [
        f"# {title}",
        "> Assembled by `scripts/manuscript/assemble.py` from the section drafts; edit the drafts, not this file. "
        "Bracketed file names are provenance pointers for the authors and are removed at submission.",
        "## Abstract", block(ab, "## Abstract (~190 words)"),
        "## Introduction", block(intro, "## Introduction"),
        "## Results",
        block(bm, "## Results"), block(lo, "## Results"), block(en, "## Results"),
        "## Discussion", block(intro, "## Discussion"),
        "## Methods",
        "### Evaluation framework", demote(block(bm, "## Methods (evaluation framework)"), 2),
        "### Annotator, held-out reads and unified annotation", demote(block(lo, "## Methods (additions)"), 2),
        "### End analyses", demote(block(en, "## Methods (additions)"), 2),
        "## Figure legends", "\n\n".join(legends[k] for k in order if k in legends),
        "## References (to be completed and formatted)", refs,
    ]
    txt = "\n\n".join(p for p in parts if p is not None) + "\n"
    open(os.path.join(D, "MANUSCRIPT_DRAFT.md"), "w").write(txt)
    missing = [k for k in order if k not in legends]
    words = len(re.sub(r"\[`[^\]]*\]", "", txt).split())
    print({"words": words, "missing_legends": missing})


if __name__ == "__main__":
    main()
