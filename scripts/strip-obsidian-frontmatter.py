#!/usr/bin/env python3
"""Git `clean` filter — strip Obsidian-injected metadata from staged content.

This repository lives inside an Obsidian vault (`Documents/Personal Projects`), whose
maintenance script `<vault>/.obsidian/add_properties.py` stamps a YAML frontmatter block
(`project:` / `type:` / `authorship:` / `tags:`) into `CLAUDE.md`, `README.md` and
`docs/*.md`, and may append a trailing "## Related notes (Obsidian)" wiki-link section.
It runs from a Claude Code PostToolUse hook, *detached*, so a tracked file can change on
disk seconds after any edit, with no visible output. It also rewrites files as CRLF.

Those files are git-tracked, so the additions show up as spurious diffs and would put
vault metadata into the project's public history.

This filter lets the working tree KEEP the metadata (Obsidian needs it to index the
notes) while git stores the clean version. Wiring, via `.gitattributes`
(`filter=obsidian-strip`) plus per-clone local git config:

    git config filter.obsidian-strip.clean "python scripts/strip-obsidian-frontmatter.py"
    git config filter.obsidian-strip.smudge cat

`sh scripts/setup-clone.sh` sets both, along with `core.hooksPath`.

Design guarantees
-----------------
* **Idempotent** — running it on already-clean content is a no-op.
* **Conservative** — removes a LEADING frontmatter block only when it carries an Obsidian
  marker (a `project:` / `type:` / `authorship:` key, or a `project/` / `type/` tag), and
  a TRAILING block only when introduced by the literal "## Related notes (Obsidian)"
  heading. Any other frontmatter is legitimate content and is left alone.
* **Fail-safe** — on ANY error it passes the input through unchanged. A clean filter that
  raises would otherwise truncate the blob and silently commit an empty file.
* **EOL-normalizing** — the stamper writes CRLF; the repo stores LF. Without this the
  files stay perpetually "modified" against the index even once the frontmatter is gone.

Reads stdin, writes stdout. Never import this; git invokes it as a subprocess.
"""
import re
import sys

FOOTER_MARKER = "## Related notes (Obsidian)"

#: Keys / tag prefixes that identify the block as vault metadata rather than real content.
_OBSIDIAN_KEY = re.compile(r"(?m)^\s*(project|type|authorship)\s*:")


def strip_frontmatter(text: str) -> str:
    m = re.match(r"^---\r?\n(.*?)\r?\n---[ \t]*\r?\n", text, re.DOTALL)
    if not m:
        return text
    block = m.group(1)
    is_obsidian = (
        _OBSIDIAN_KEY.search(block) is not None
        or "project/" in block
        or "type/" in block
    )
    if not is_obsidian:
        return text  # some other, legitimate frontmatter — leave it alone
    rest = text[m.end():]
    # The stamper writes "---\n" + "\n" + body; drop the one blank separator line.
    return re.sub(r"^\r?\n", "", rest, count=1)


def strip_footer(text: str) -> str:
    idx = text.rfind(FOOTER_MARKER)
    if idx == -1:
        return text
    head = text[:idx].rstrip()
    # Drop the horizontal rule the linker writes above the heading.
    if head.endswith("---"):
        head = head[:-3].rstrip()
    return head + "\n"


def strip(text: str) -> str:
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return strip_footer(strip_frontmatter(text))


def main() -> None:
    data = sys.stdin.buffer.read()
    try:
        out = strip(data.decode("utf-8")).encode("utf-8")
    except Exception:
        out = data  # never fail closed — pass through unchanged
    sys.stdout.buffer.write(out)


if __name__ == "__main__":
    main()
