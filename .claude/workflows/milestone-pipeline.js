export const meta = {
  name: 'milestone-pipeline',
  description: 'Run one ROADMAP.md milestone through Research -> Implement -> Adversarial -> Improve, editing the real repo, with the full test suite as the gate.',
  phases: [
    { title: 'Research', detail: 'read the milestone spec + goal + touched code; produce an exact implementation plan (signatures, tests with literature-anchored values, invariants, risks)' },
    { title: 'Implement', detail: 'write the code + tests per the plan; run the full suite; report green/diff' },
    { title: 'Adversarial', detail: 'independently review vs acceptance criteria + inviolable invariants; re-run the suite; confirm no pinned value changed' },
    { title: 'Improve', detail: 'apply must-fix findings; re-run the suite; report final state' },
  ],
}

// Milestone is hardcoded here (args does not propagate reliably via scriptPath).
// Edit MS before each run.  args overrides are still honored if present.
const REPO = (args && args.repo) ? args.repo : 'C:/Users/cedar/Documents/Personal Projects/Source Code/stability-mflds'
const EPIC = (args && args.epic) ? args.epic : 'E1'
const MS = (args && args.milestone) ? args.milestone : 'E11-M3'

const TOOLING = [
  'ENVIRONMENT & TOOLING RULES (this host):',
  ' - The repo root is: ' + REPO + '  (note the SPACE in "Source Code" — always quote the path).',
  ' - Make ALL file changes with the Write/Edit tools (they hit the real filesystem reliably). Do NOT use shell rm/mkdir/cp/mv',
  '   for source files — a prior stale-path shell mishap made this the safe rule.',
  ' - Run the test suite with PowerShell, from the repo root, EXACTLY:',
  '       Set-Location "' + REPO + '"; python -m pytest -q',
  '   (Read the last ~15 lines: the summary line with pass/fail counts.) Use PowerShell (not Bash) for filesystem/test ops here.',
].join('\n')

const INVARIANTS = [
  'INVIOLABLE INVARIANTS (from CLAUDE.md — violating any is a hard failure):',
  ' 1. EXACT ARITHMETIC: every invariant stays fractions.Fraction. Float ONLY at the display/geometry layer',
  '    (Wall.radius, central_charge, bg_boundary_curve grid, alpha_crit sqrt, viz). Never introduce float into the core.',
  ' 2. Coskun-Huizenga discriminant convention is FIXED: Delta = 1/2 mu^2 - ch2/(rd), d=H^2. discriminant_brief (=2 Delta)',
  '    is comparison-only; never the default.',
  ' 3. ZERO RUNTIME DEPENDENCIES: "import bridgeland_stability" MUST succeed with only the Python stdlib. matplotlib/plotly',
  '    (and any CAS) stay lazily imported inside viz/ (or oracle/) — never at module top level. Any new core module must be',
  '    stdlib-only at import time.',
  ' 4. FROZEN DATACLASSES: create new instances; do not mutate. New optional fields must have defaults (backward-compatible).',
  ' 5. THE 42 PINNED TESTS ARE GROUND TRUTH. If code produces a different answer, the CODE is wrong. NEVER edit or weaken a',
  '    pinned test to make it pass. This milestone is in Epic 1 ("Harden and document the existing code") which introduces',
  '    NO NEW MATHEMATICS: NO computed value may change — only tests, citations, provenance metadata, and docs are added.',
  ' 6. Every asserted test value must be literature-anchored or exactly re-derived; two-way verification (exact Fraction',
  '    recompute + primary-source citation) for anything mathematical. docs/CORRECTIONS.md is the template.',
  'The whole suite must stay green (currently 42 passing); this milestone should ADD tests (strictly more than 42) and keep',
  'all prior tests passing unchanged.',
].join('\n')

const CONTEXT = [
  'You are working in the bridgeland_stability repo at: ' + REPO,
  'It is a pure-Python, zero-dependency, exact-Fraction toolkit for Bridgeland stability. The roadmap is docs/ROADMAP.md;',
  'the goal set is docs/GOALS.md; corrections ledger is docs/CORRECTIONS.md; agent guidance is CLAUDE.md.',
  'You are implementing milestone ' + MS + ' of epic ' + EPIC + '. READ its full section from docs/ROADMAP.md (find the header',
  '"### ' + MS + '") and the goal it maps to in docs/GOALS.md, plus the actual source files it touches.',
  '', TOOLING, '', INVARIANTS,
].join('\n')

// ---------------- Research ----------------
const RESEARCH_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    plan_markdown: { type: 'string', description: 'Concrete implementation plan: exact functions/signatures, exact file edits, and the precise test cases with literature-anchored asserted values.' },
    files_to_modify: { type: 'array', items: { type: 'string' } },
    new_test_files: { type: 'array', items: { type: 'string' } },
    test_cases: { type: 'array', items: { type: 'object', additionalProperties: false, properties: { name: { type: 'string' }, file: { type: 'string' }, asserts: { type: 'string' } }, required: ['name', 'file', 'asserts'] } },
    acceptance_criteria: { type: 'array', items: { type: 'string' } },
    invariant_risks: { type: 'array', items: { type: 'string' }, description: 'Where this milestone could accidentally violate an invariant, and how to avoid it.' },
    open_questions: { type: 'array', items: { type: 'string' } },
  },
  required: ['plan_markdown', 'files_to_modify', 'new_test_files', 'test_cases', 'acceptance_criteria', 'invariant_risks', 'open_questions'],
}

phase('Research')
const research = await agent([
  'STAGE 1 — RESEARCH (read-only planning). Do NOT edit any file.',
  '', CONTEXT, '',
  'Read the ' + MS + ' section of docs/ROADMAP.md, the mapped goal in docs/GOALS.md, and every source/test file this milestone',
  'touches. Produce a precise, contributor-ready implementation plan: exact function names/signatures, exact edits (which strings',
  'change and to what), the new test file(s) and every test case with its EXACT asserted value (literature-anchored or exactly',
  're-derived — recompute by hand and show the arithmetic for any numeric assertion), the acceptance criteria copied from the',
  'roadmap, and the specific places this milestone risks violating an invariant (esp. "no computed value changes",',
  '"zero-dependency import", "frozen dataclass defaults"). List any genuine open question. Ground every claim in the real files.',
].join('\n'), { label: MS + ':research', phase: 'Research', schema: RESEARCH_SCHEMA, agentType: 'general-purpose' })

// ---------------- Implement ----------------
const IMPL_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    summary: { type: 'string' },
    files_changed: { type: 'array', items: { type: 'string' } },
    tests_added: { type: 'array', items: { type: 'string' } },
    pytest_tail: { type: 'string', description: 'The last ~15 lines of the pytest -q output (the summary line with counts).' },
    all_green: { type: 'boolean' },
    total_tests: { type: 'number' },
    notes: { type: 'string', description: 'Anything the reviewer should scrutinize.' },
  },
  required: ['summary', 'files_changed', 'tests_added', 'pytest_tail', 'all_green', 'total_tests', 'notes'],
}

phase('Implement')
const impl = await agent([
  'STAGE 2 — IMPLEMENT. Apply the plan to the REAL repo with Edit/Write.',
  '', CONTEXT, '',
  'Implement milestone ' + MS + ' exactly per the RESEARCH PLAN below. Edit the real files under ' + REPO + ' using Write/Edit',
  '(never shell fs commands). Add the new test file(s) with the exact literature-anchored assertions. Then RUN THE FULL SUITE',
  '(Set-Location "' + REPO + '"; python -m pytest -q) and confirm it is green with strictly MORE than 42 tests. If any prior',
  '(pinned) test fails, your code is wrong — FIX THE CODE, never the pinned test. Do not add any runtime dependency; keep new',
  'core modules stdlib-only at import time; keep frozen dataclasses backward-compatible (new fields need defaults). Report the',
  'pytest summary line verbatim, the files changed, and anything the adversarial reviewer should scrutinize. If you could not get',
  'green, say so honestly in notes and set all_green=false.',
  '',
  '===== RESEARCH PLAN =====',
  (research && research.plan_markdown) ? research.plan_markdown : '[research plan missing — read the ' + MS + ' section of docs/ROADMAP.md and proceed from it]',
  '',
  'ACCEPTANCE CRITERIA: ' + ((research && research.acceptance_criteria) ? research.acceptance_criteria.join(' | ') : '(see ROADMAP)'),
].join('\n'), { label: MS + ':implement', phase: 'Implement', schema: IMPL_SCHEMA, agentType: 'general-purpose' })

// ---------------- Adversarial ----------------
const ADV_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    verdict: { type: 'string', enum: ['accept', 'accept-with-fixes', 'reject'] },
    verdict_rationale: { type: 'string' },
    pinned_values_unchanged: { type: 'boolean' },
    zero_dependency_import_ok: { type: 'boolean' },
    all_green: { type: 'boolean' },
    total_tests: { type: 'number' },
    acceptance_check: { type: 'array', items: { type: 'object', additionalProperties: false, properties: { criterion: { type: 'string' }, met: { type: 'boolean' }, note: { type: 'string' } }, required: ['criterion', 'met', 'note'] } },
    findings: { type: 'array', items: { type: 'object', additionalProperties: false, properties: { severity: { type: 'string', enum: ['blocker', 'major', 'minor', 'nit'] }, issue: { type: 'string' }, fix: { type: 'string' } }, required: ['severity', 'issue', 'fix'] } },
    must_fix: { type: 'array', items: { type: 'string' } },
  },
  required: ['verdict', 'verdict_rationale', 'pinned_values_unchanged', 'zero_dependency_import_ok', 'all_green', 'total_tests', 'acceptance_check', 'findings', 'must_fix'],
}

phase('Adversarial')
const adv = await agent([
  'STAGE 3 — ADVERSARIAL REVIEW. Be rigorous, specific, and independent. You may read anything and RUN read-only commands',
  '(PowerShell), but do NOT edit source files (the Improve stage applies fixes) and do NOT use shell fs-mutation commands.',
  '', CONTEXT, '',
  'Independently verify the implementation of ' + MS + ':',
  ' - RUN the full suite yourself (Set-Location "' + REPO + '"; python -m pytest -q). Report all_green and total_tests. Confirm',
  '   it is strictly more than 42 and every prior pinned test still passes.',
  ' - Confirm NO COMPUTED VALUE CHANGED: this is an Epic-1 hardening milestone. Read each touched SOURCE file (not tests) and',
  '   verify no formula/constant/pinned output was altered — only citations, notes, docstrings, new optional fields (with',
  '   defaults), or new additive functions. Set pinned_values_unchanged accordingly.',
  ' - Confirm the zero-dependency import still holds: run  python -c "import bridgeland_stability"  and verify any new core',
  '   module is stdlib-only at import time (no top-level matplotlib/plotly/CAS import).',
  ' - Check EACH acceptance criterion from the roadmap and mark met/not-met with a note.',
  ' - Scrutinize the new tests: do their asserted values actually match the literature / an exact re-derivation? Recompute at',
  '   least the load-bearing ones by hand. Catch any test that would pass while asserting a WRONG value (a real risk in a',
  '   citation-accuracy project). Check frozen-dataclass defaults, convention correctness, and over-claims.',
  ' - Verify citations/allowlist content against primary sources where relevant (load web tools via ToolSearch',
  '   "select:WebSearch,WebFetch" if needed).',
  'Give a per-criterion acceptance_check, findings (severity + concrete fix), a must_fix list (blockers/majors the Improve stage',
  'MUST apply), and an overall verdict. Do not manufacture objections.',
  '',
  '===== IMPLEMENTATION SUMMARY =====',
  (impl && impl.summary) ? impl.summary : '[missing]',
  'Files changed: ' + ((impl && impl.files_changed) ? impl.files_changed.join(', ') : '?') + ' | tests added: ' + ((impl && impl.tests_added) ? impl.tests_added.join(', ') : '?'),
  'Implementer pytest tail: ' + ((impl && impl.pytest_tail) ? impl.pytest_tail : '?'),
  'ACCEPTANCE CRITERIA: ' + ((research && research.acceptance_criteria) ? research.acceptance_criteria.join(' | ') : '(see ROADMAP)'),
].join('\n'), { label: MS + ':adversarial', phase: 'Adversarial', schema: ADV_SCHEMA, agentType: 'general-purpose' })

// ---------------- Improve ----------------
const IMPROVE_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    summary: { type: 'string' },
    changes_applied: { type: 'array', items: { type: 'string' } },
    all_green: { type: 'boolean' },
    total_tests: { type: 'number' },
    pytest_tail: { type: 'string' },
    remaining_issues: { type: 'array', items: { type: 'string' } },
    ready: { type: 'boolean', description: 'True iff every must_fix was resolved, the suite is green, and acceptance is met.' },
  },
  required: ['summary', 'changes_applied', 'all_green', 'total_tests', 'pytest_tail', 'remaining_issues', 'ready'],
}

phase('Improve')
const mustFix = (adv && adv.must_fix && adv.must_fix.length) ? ('- ' + adv.must_fix.join('\n- ')) : '(none)'
const findings = (adv && adv.findings) ? adv.findings.filter(f => f.severity === 'blocker' || f.severity === 'major').map(f => '- [' + f.severity + '] ' + f.issue + ' -> ' + f.fix).join('\n') : '(none)'
const improve = (adv && adv.verdict === 'accept' && (!adv.must_fix || adv.must_fix.length === 0))
  ? { summary: 'No fixes required; adversarial verdict was accept with no must-fix items.', changes_applied: [], all_green: adv.all_green, total_tests: adv.total_tests, pytest_tail: '(unchanged since adversarial pass)', remaining_issues: [], ready: adv.all_green }
  : await agent([
      'STAGE 4 — IMPROVE. Apply the adversarial must-fix items to the REAL repo with Edit/Write, then re-run the full suite.',
      '', CONTEXT, '',
      'Apply every MUST-FIX and every blocker/major finding below for milestone ' + MS + '. Respect all invariants (no computed',
      'value change, zero-dependency import, frozen-dataclass defaults, pinned tests never weakened). Then RUN the full suite',
      '(Set-Location "' + REPO + '"; python -m pytest -q) and confirm green with > 42 tests. Report exactly what you changed, the',
      'pytest summary line, any remaining issues, and whether the milestone is ready (all must-fix resolved + green + acceptance met).',
      '',
      'MUST-FIX:\n' + mustFix,
      '\nBLOCKER/MAJOR FINDINGS:\n' + findings,
    ].join('\n'), { label: MS + ':improve', phase: 'Improve', schema: IMPROVE_SCHEMA, agentType: 'general-purpose' })

return {
  milestone: MS,
  research_plan_present: !!(research && research.plan_markdown),
  acceptance_criteria: (research && research.acceptance_criteria) ? research.acceptance_criteria : [],
  implement: impl ? { files_changed: impl.files_changed, tests_added: impl.tests_added, all_green: impl.all_green, total_tests: impl.total_tests } : null,
  adversarial: adv ? { verdict: adv.verdict, pinned_values_unchanged: adv.pinned_values_unchanged, zero_dependency_import_ok: adv.zero_dependency_import_ok, all_green: adv.all_green, total_tests: adv.total_tests, acceptance_check: adv.acceptance_check, must_fix: adv.must_fix, findings: adv.findings } : null,
  improve: improve ? { changes_applied: improve.changes_applied, all_green: improve.all_green, total_tests: improve.total_tests, remaining_issues: improve.remaining_issues, ready: improve.ready, pytest_tail: improve.pytest_tail } : null,
  gate_ready: !!(improve && improve.ready && improve.all_green),
}
