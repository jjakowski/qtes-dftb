# Restart prompt

Paste this into a new agent session (Claude Code, Codex, ...) started in
`/home/j2c/Development/QTES/qtes-dftb`:

```text
Read AGENTS.md, then agent_state.md, then the last entries of
agent_progress.md. Run `git status` and `git log --oneline -10`. Activate
the rmg-dev environment and check that the code builds as described in
AGENTS.md (do not run long tests).

Then summarize in your own words: the current state, the current milestone
and next step, the open questions, and anything uncommitted or inconsistent
you found. Do not change any files until I approve a plan.
```

## To continue a specific step

After the summary, give the step, e.g.:

```text
Proceed with milestone 1, step 1 as described in agent_state.md. Explain
the exact commands first and wait for my approval. Afterwards update
agent_state.md and agent_progress.md in the same commit as the work.
```

## Notes

- The agent's own memory from earlier sessions may not be available. All
  needed information must be in the files above, `docs/KNOWN_ISSUES.md`,
  `docs/DFTBPLUS_AUDIT.md`, and the git history.
- Starting Codex, as in the DFTB+ project:
  `source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev && cd /home/j2c/Development/QTES/qtes-dftb && codex .`
