# Security Policy

## Scope of this document

`bridgeland_stability` is a pure-mathematical-computation library. Its core
performs local, exact arithmetic (`fractions.Fraction`) and has **no network
I/O, no subprocess execution, and no deserialization of untrusted data** — the
core imports only the Python standard library. As a result the meaningful
security surface is narrow and quite different from a web service or a
data-parsing library.

**Mathematical correctness bugs** — a wrong result, an incorrect formula, or a
discrepancy with the cited literature — are the primary failure mode of a
library like this. Those are **not** security vulnerabilities; please report
them as regular GitHub issues. The project maintains
[`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) precisely for this and takes
correctness seriously.

This policy covers genuine security vulnerabilities: potential for arbitrary
code execution, privilege escalation, or other unsafe behavior that could harm
someone who installs or runs the package.

## Supported versions

This project is at an early stage; only the latest published release is
actively maintained.

| Version | Supported |
|---------|-----------|
| 0.1.x (latest) | ✅ |
| < 0.1.0 | ❌ |

## Reporting a vulnerability

Please **do not** open a public issue for a security vulnerability — that could
expose it before a fix exists.

- **Preferred — GitHub Private Vulnerability Reporting.** Use the
  **Security → Advisories → Report a vulnerability** button on the repository.
  This keeps the report confidential and allows a fix to be drafted privately.
- **Fallback — email.** If you cannot use GitHub's private reporting, email
  **me@chrisdare.net** with the subject `[bridgeland-stability] Security`. PGP
  is not required; please describe the issue and include a minimal reproducer if
  you can.

## Out of scope

The following are **not** treated as security vulnerabilities under this policy:

- **Mathematical or algorithmic errors** (wrong discriminant values, incorrect
  wall formulas, disagreement with the literature). Report as regular issues;
  see [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md).
- **Exceptions on malformed or out-of-range input** (e.g. `ValueError`,
  `ZeroDivisionError`). The library is not designed to process untrusted
  external data; callers supply their own mathematical parameters.
- **Resource consumption** reachable only by a caller who can already execute
  arbitrary Python. If you can run Python in the process, you can exhaust
  resources far more simply than through this library.
- **Vulnerabilities in the optional visualization dependencies** (plotly,
  matplotlib). Those are third-party projects with their own advisories — see
  below.

## Dependencies and attack surface

The core (`pip install bridgeland-stability`) has **zero runtime dependencies**;
it uses only `fractions`, `math`, `dataclasses`, and `typing`. No network
requests, subprocess calls, `eval`/`exec`, or deserialization occur anywhere in
the core.

The optional `[viz]` extra (`pip install "bridgeland-stability[viz]"`) adds
**plotly ≥ 5.18** and **matplotlib ≥ 3.7**, which have their own security
advisories. If you install the visualization extras, track upstream advisories
for those packages independently (e.g. with `pip-audit` or Dependabot). Note
that Plotly's interactive HTML output references CDN-hosted assets, which the
**browser** fetches when you open the HTML — the Python process itself makes no
network requests.

## Response expectations

This is a solo project maintained in spare time. There is no security team and
no formal SLA. On a best-effort basis:

- Reports submitted via GitHub Private Vulnerability Reporting or email will be
  acknowledged within about **7 days**.
- If confirmed, a fix is prioritized by severity and maintainer availability.
- You will be kept informed of progress and credited in the release notes and
  advisory, unless you prefer to remain anonymous.

## Coordinated disclosure

The maintainer will:

1. Confirm receipt and whether the report is accepted as a security issue.
2. Work with the reporter on a coordinated fix before any public disclosure.
3. Publish a fix and a GitHub Security Advisory at or before public disclosure
   of the details.
4. Credit the reporter in the advisory (unless anonymity is requested).

There is no fixed embargo window; the timeline is agreed with the reporter
case by case, with a good-faith effort toward timely remediation.
