# Compute environment

Last verified: 2026-09-29 (partly inherited from the melanoma repository;
re-verify for this project).

- The user works in Ubuntu under WSL2 with PyCharm. In the sibling melanoma
  repository WSL has 7.7 GB RAM; this project has not measured it.
- The Claude Science app runs natively on Windows. Its sandbox has its own
  Windows-side Python and R, not the project venv, and `git` does not run in
  it. The agent therefore writes files and gives commands; the user runs
  them in WSL and runs git.
- Repository path: `D:\#University\#MS-SystemsBio\Master_Thesis\repos\networkbio`
  on Windows, `/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio`
  in WSL. Always quote it.
- WSL distribution: Ubuntu 24.04 (Python 3.12.3). Kernel string from the
  toy run: `6.18.40.1-microsoft-standard-WSL2`.
- The project venv exists (2026-10-06): `.venv`, built from
  `requirements.lock.txt` by `scripts/setup/create_venv.sh`; versions in
  `logs/pip_freeze_2026-10-06.txt`.
- MILP solver: HiGHS (highspy 1.15.1) passed the toy problem; SCIP is also
  installed. Record time limit and optimality gap with each real result.
- The Cowork desktop agent cannot type into the WSL terminal (Windows
  Terminal is view/click only for it) and its own Linux VM is not WSL. The
  user runs the pipeline commands; the agent reads logs and results.
- The repository is on the Windows filesystem (`/mnt/d`), which is slow for
  many small files; expect sluggish imports.
