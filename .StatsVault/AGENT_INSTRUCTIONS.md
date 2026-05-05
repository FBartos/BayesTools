# StatsVault Agent Instructions

When asked to write literature-backed text with StatsVault, follow this workflow:

1. Read `.StatsVault/project-memory.md`, `.StatsVault/key-papers.md`, `.StatsVault/ACTION_REQUIRED.md`, and `.StatsVault/research-log.md`.
2. Start with `prepare_writing_context(project_root=..., task=..., query=..., claim=...)` when you are working toward a specific proposition. Use `claim=` whenever possible.
3. If needed, search directly with `search_papers(..., project_root=..., task=...)`, `discover_claim_support(...)`, or `find_author_works(...)`.
4. If discovery is returning several strong candidates, foundational leads, or contrastive leads, keep materializing and preparing them. Do not stop after the first acceptable citation if the search vein is still productive.
5. Run `report_task_coverage(...)` before drafting or stopping research so blind spots, missing foundational work, and missing comparison papers are visible.
6. Materialize every paper you plan to rely on with `materialize_papers(...)`.
7. If a paper has a PDF but no canonical summary, use `prepare_contribution(...)` to generate the extraction command — this runs in a separate terminal and uploads both the PDF and summary to the library, and writes project-local copies into `.StatsVault/contributions/` when `project_root` is supplied.
   If the user asks for a local summary script/command for a PDF in `.StatsVault/unprocessed-input/`, do not improvise. Call `prepare_contribution(project_root=..., pdf_path="filename.pdf", extraction_cli="copilot")` when the user says GitHub Copilot; otherwise omit `extraction_cli`. The tool always returns the operator command and that command always contributes upstream. Copilot defaults to `reasoning_effort="xhigh"`.
8. Draft only from local copied files in `.StatsVault/resources/` or `.StatsVault/contributions/` and cite only keys present in `.StatsVault/references.bib`.
9. Save literature-backed `.tex` fragments with `write_grounded_tex_fragment(...)` or use the one-shot `prepare_and_write_grounded_tex_fragment(...)` path. Include `section` and `claim` when you know them.
10. When a source is not citation-ready, insert `\SV{...}` and surface `.StatsVault/ACTION_REQUIRED.md` instead of inventing support.
11. Use `resolve_pending_markers(...)` to turn unresolved `\SV{...}` markers back into concrete search/materialization steps.
12. Export `.StatsVault/references.bib` to a manuscript-root `.bib` file with `export_project_bib(...)` when the draft is ready to compile.
13. Replace resolved markers with grounded text via `replace_resolved_markers(...)`.
14. After substantial research, record what changed with `append_project_note(...)`.
15. Run `validate_manuscript_sources(...)` before considering the draft complete.
16. Re-check `.StatsVault/pending-markers.json` until no unresolved manuscript markers remain.

Hard rules:
- Do not cite papers that have not been materialized into `.StatsVault/resources/` or discovered from `.StatsVault/contributions/`.
- Do not cite PDF-only papers until `prepare_contribution(...)` has produced a project-local contribution summary or the library later gains a canonical summary.
- Do not invent BibTeX keys; use only keys in `.StatsVault/references.bib`.
- Do not treat search hits as sufficient evidence by themselves. Materialize first.
- Do not stop research after one usable citation if `report_task_coverage(...)` or the search results still show obvious high-value follow-ups.
