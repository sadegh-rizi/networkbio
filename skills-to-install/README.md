# skills-to-install/

Agent skills for this project. Claude Code reads skills from
`.claude/skills/<name>/SKILL.md`. Agents in this workflow do not write into
`.claude/`, so the skills are staged here and copied by the user.

From the repository root, in WSL:

```bash
mkdir -p .claude/skills
cp -r skills-to-install/* .claude/skills/
rm .claude/skills/README.md
```

| Skill | Use it for |
|---|---|
| `careful-planning` | the plan gate before any implementation |
| `network-inference` | building or interpreting a CORNETO / PKN inference |
| `statistics-review` | tests, small-sample limits, nulls, enrichment |
| `python-analysis` | reusable Python logic and parsing |
| `publication-plot` | figures, including network figures |
| `research-code-style` | keeping code readable at research level |
| `ai-provenance` | recording prompts and stamping outputs with a commit |
