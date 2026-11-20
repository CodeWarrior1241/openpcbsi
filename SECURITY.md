# Security

## Reporting a problem

Please report security problems privately through GitHub's private vulnerability reporting:
[Security → Report a vulnerability](https://github.com/CodeWarrior1241/openpcbsi/security/advisories/new).
If that is not available to you, open an issue that asks for a private contact, without the details.
Include the command, the openpcbsi commit (`git describe --always --dirty`) and, where you can, a small
input that shows the problem.

Only the latest commit on `main` is supported.

## What openpcbsi treats as untrusted

openpcbsi is a design-review tool that reads files other people send. These inputs are treated as
untrusted, and a problem with any of them is a security problem:

- board files (`.kicad_pcb`), their project files (`.kicad_pro`, `.kicad_dru`, `.kicad_sch`) and
  everything named inside them: nets, parts, layers, rails, text;
- board configurations (`-c CONFIG`), waivers, region specs and fabricator profiles (`fab-check`
  fetches the pages a profile cites, over `https://` only);
- folders received from someone else: a sign-off record with its board and evaluation
  (`verify-signoff`, `--rerun`), run folders given to `rediff`, files given to `dfm --gerbers`;
- the reports openpcbsi writes (HTML, JSON), which are opened in a browser or passed on.

What openpcbsi does with them:

- names from inputs never become paths outside the run folder (`si.safe_name`, `si.inside`) and are
  escaped in HTML;
- paths in written reports are made relative, so a report passed on does not show the author's
  folders (`paths.portable`);
- no interpreter is started that a received folder could replace: `verify-signoff --rerun` evaluates
  in-process, and the openEMS and KiPI interpreters, which must be named explicitly
  (`OPENPCBSI_OPENEMS_PYTHON`, `OPENPCBSI_KIPY_PYTHON`), run isolated (`python -I`);
- a folder is emptied only when openpcbsi made it (`dfm`'s plot folder carries a mark file);
- the file parsers and KiPI are pinned by commit and by hash, checked against their published
  releases, and every file in their installed sites is checked against its package's RECORD
  (`deps/pins.py`, `deps/setup_*.sh --check`).

A sign-off record is authenticated only by its SSH signature, checked with
`verify-signoff --allowed-signers`; without the signature a record proves only that it is consistent
with its own evaluation, not who wrote it.

KiCad itself (`kicad-cli`), openEMS and KiPI are outside this scope: report their problems upstream.
