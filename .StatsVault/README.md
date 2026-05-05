# StatsVault Workspace

This folder is a manuscript-local cache backed by the main StatsVault library.

Use it to keep all literature-backed writing grounded in local copied resources.

Key files:
- `manifest.json`: manuscript-local state, including prepared sources and claim records
- `query-cache/`: recent search results — slim `index.json` plus one CSV per search
- `resources/`: copied library-backed PDFs, canonical summaries, BibTeX snippets, and metadata
- `contributions/`: project-local contributed PDFs, contributed summaries, provisional metadata, and `contributions.bib`
- `references.bib`: manuscript-local bibliography built only from citation-ready sources
- `tex/statsvault-macros.tex`: TeX support for unresolved `\SV{...}` markers
- `pending-markers.json`: durable tracker of unresolved `\SV{...}` manuscript gaps
- `provenance.jsonl`: append-only record of contributed summaries, grounded writes, marker replacements, and bibliography exports, including task/section/claim context
- `project-memory.md`: persistent working memory for the agent
- `research-log.md`: human-readable log of searches and materialization steps
- `ACTION_REQUIRED.md`: human-required and agent-required source-preparation tasks
