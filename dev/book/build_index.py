import collections
import html
import json
import pathlib
import re
import subprocess
import hashlib

repo = pathlib.Path(__file__).resolve().parents[2]
import os
os.chdir(repo)
book = repo / 'docs' / 'book'
out = book / 'downloads'
out.mkdir(exist_ok=True)
subprocess.run(['Rscript', '--vanilla', str(repo / 'dev/book/inventory.R'), str(out / 'r-functions.json')], check=True)
data = json.loads((out / 'r-functions.json').read_text())
rows = data['functions']
ns = pathlib.Path('NAMESPACE').read_text()
exports = set(re.findall(r'^export\((.*?)\)', ns, re.M))
methods = {a.strip('"') + '.' + b for a, b in re.findall(r'^S3method\((.*?),(.*?)\)', ns, re.M)}
for row in rows:
    row['language'] = 'R'
    if row['nested']:
        row['category'] = 'Local R helpers'
    elif row['name'] in exports:
        row['category'] = 'Public R functions'
    elif row['name'] in methods:
        row['category'] = 'Registered R methods'
    elif row['file'].endswith('extendr-wrappers.R'):
        row['category'] = 'R-to-Rust bindings'
    elif '/legacy_' in row['file']:
        row['category'] = 'Legacy internal R'
    else:
        row['category'] = 'Modern internal R'

# Rust declarations are grouped by their enclosing impl. The explicit tests
# modules and cfg(test) method are separated from production definitions.
for path in sorted(pathlib.Path('src/rust/src').glob('*.rs')):
    lines = path.read_text().splitlines()
    impl = ''
    in_tests = False
    test_next = False
    for i, line in enumerate(lines):
        if re.match(r'^mod tests\s*\{', line):
            in_tests = True
            impl = ''
        match = re.match(r'^impl(?:<[^>]*>)?\s+([A-Za-z_][A-Za-z_0-9]*)(?:<[^>]*>)?\s*\{', line)
        if match:
            impl = match.group(1)
        elif line == '}':
            impl = ''
        if line.strip() == '#[cfg(test)]':
            test_next = True
        match = re.match(r'^\s*(?:pub(?:\([^)]*\))?\s+)?fn\s+(\w+)\s*(?:<[^;]*?>)?\s*\(', line)
        if not match:
            continue
        name = match.group(1)
        if impl:
            name = impl + '::' + name
        is_test = in_tests or test_next
        rows.append(dict(name=name, file=str(path), line=i+1,
            language='Rust', nested=False, scope='test configuration' if is_test else impl,
            category='Rust tests and test helpers' if is_test else 'Rust implementation'))
        test_next = False

sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
counts = collections.Counter(r['category'] for r in rows)
assert {r['name'] for r in rows if r['category'] == 'Public R functions'} == exports
assert {r['name'] for r in rows if r['category'] == 'Registered R methods'} == methods
categories = ['Public R functions', 'Registered R methods', 'Modern internal R',
              'Legacy internal R', 'R-to-Rust bindings', 'Local R helpers',
              'Rust implementation', 'Rust tests and test helpers']
e = html.escape
descriptions = {
 'Public R functions': 'User-facing functions exported in NAMESPACE. Worked examples belong in the main book.',
 'Registered R methods': 'Demonstrate these through their generic operations, such as print(), plot(), predict(), subsetting, and concatenation.',
 'Modern internal R': 'Package-level implementation functions, including scientific algorithms, validation, formatting, I/O, and runtime support.',
 'Legacy internal R': 'Private helpers supporting the frozen create_phenotypes() engine.',
 'R-to-Rust bindings': 'Seven R wrappers for compiled Rust entry points; these are not seven additional scientific methods.',
 'Local R helpers': 'Named functions defined inside another function. Scope identifies the enclosing function; these are not independently callable package functions.',
 'Rust implementation': 'Functions and associated methods in the package-owned Rust source. Rust pub visibility does not make a function part of the exported R API.',
 'Rust tests and test helpers': 'Verification code, separated from implemented methods. Includes the test-specific Bits::tail_is_clean definition.'}
sections = []
for category in categories:
    byfile = collections.defaultdict(list)
    for row in rows:
        if row['category'] == category:
            byfile[row['file']].append(row)
    groups = []
    for file, entries in sorted(byfile.items()):
        trs = []
        for row in sorted(entries, key=lambda r: (r['name'], r['line'])):
            url = f'https://github.com/samuelbfernandes/simplePHENOTYPES/blob/{sha}/{file}#L{row["line"]}'
            trs.append(f'<tr><td><code>{e(row["name"])}</code></td><td>{e(row["scope"]) or "—"}</td><td><a href="{url}">{e(file)}:{row["line"]}</a></td></tr>')
        groups.append(f'<details open><summary>{e(file)} <span>{len(entries)}</span></summary><table><thead><tr><th>Function / method</th><th>Enclosing scope</th><th>Source</th></tr></thead><tbody>{"".join(trs)}</tbody></table></details>')
    sections.append(f'<section data-category="{e(category)}"><h2>{e(category)} <span>{counts[category]}</span></h2><p>{descriptions[category]}</p>{"".join(groups)}</section>')
summary = ''.join(f'<tr><td>{e(c)}</td><td>{counts[c]}</td></tr>' for c in categories)
page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>simplePHENOTYPES — function inventory</title><style>
:root{font-family:system-ui,sans-serif;color:#182c39;background:#f5f7fa}body{max-width:1200px;margin:auto;padding:32px}h1{font-size:34px;letter-spacing:-1px}h2{margin-top:42px}p{max-width:900px;line-height:1.6}code{font-size:13px;color:#11465a}a{color:#176282}table{border-collapse:collapse;width:100%;background:white}td,th{text-align:left;padding:9px 12px;border-bottom:1px solid #e2e8ed;overflow-wrap:anywhere}th{font-size:12px;text-transform:uppercase;color:#51616d}td:nth-child(1){width:35%}td:nth-child(2){width:25%;font-size:13px}td:nth-child(3){font-size:12px}summary{cursor:pointer;padding:14px 12px;background:#e9eff3;font-weight:600}details{margin:12px 0;border:1px solid #d4dfe5;border-radius:6px;overflow:hidden}span{font-size:13px;color:#56707f;margin-left:8px}.controls{position:sticky;top:0;background:#f5f7fa;padding:15px 0;display:flex;gap:10px;flex-wrap:wrap;border-bottom:1px solid #cbd9df}input,select,button{padding:10px;border:1px solid #a5bac5;border-radius:5px;font:inherit}input{flex:1;min-width:220px}.badge{color:#356179;font-weight:600;font-size:13px}.summary{max-width:650px}.note{border-left:4px solid #277f86;padding:12px 18px;background:#e9f4f3}footer{padding:28px 0;color:#51616d} @media(max-width:700px){body{padding:15px}td,th{padding:7px}h1{font-size:28px}}
</style><body><div class="badge">DOCUMENTATION SCOPE • SOURCE INVENTORY</div><h1>simplePHENOTYPES: all named functions</h1>
<p>Inventory of package-owned R and Rust source at revision <code>SHA</code>. Generated from the current working tree. R functions were extracted with R's parser; Rust declarations were scanned from <code>src/rust/src/</code>. Source links are pinned to this revision.</p>
<p class="note">This inventory supports the decision about documentation depth. It is not a theory audit. A function count is not a method count: one public function may expose several scientific methods, and one method may require many internal helpers.</p>
<table class="summary"><thead><tr><th>Category</th><th>Definitions</th></tr></thead><tbody>SUMMARY</tbody></table>
<p>Anonymous R callbacks (ANONYMOUS) are excluded from the named list. Vendored crates, generated native registration glue, tests outside Rust source, development scripts, and historical <code>context/</code> sources are outside scope. Rust test functions are included in their own section. R wrappers and their Rust kernels are separate definitions of the same cross-language interface.</p>
<div class="controls"><input id="search" type="search" aria-label="Search functions" placeholder="Search function, enclosing scope, or file…"><select id="category" aria-label="Filter category"><option value="">All categories</option>OPTIONS</select><button id="expand">Expand all</button><button id="collapse">Collapse all</button><span id="matches" aria-live="polite"></span></div>
SECTIONS
<footer><a href="function-inventory.json" download>Download the full machine-readable inventory</a> · Package-owned source definitions.</footer>
<script>
const search=document.querySelector('#search'),category=document.querySelector('#category');
function filter(){let n=0;for(const s of document.querySelectorAll('section')){let total=0;const chosen=!category.value||s.dataset.category===category.value;for(const d of s.querySelectorAll('details')){let visible=0;for(const r of d.querySelectorAll('tbody tr')){const ok=chosen&&r.textContent.toLowerCase().includes(search.value.toLowerCase());r.hidden=!ok;if(ok)visible++}d.hidden=!visible;if(search.value&&visible)d.open=true;total+=visible}s.hidden=!total;n+=total}document.querySelector('#matches').textContent=n+' definitions'}
search.addEventListener('input',filter);category.addEventListener('change',filter);document.querySelector('#expand').onclick=()=>document.querySelectorAll('details').forEach(d=>d.open=true);document.querySelector('#collapse').onclick=()=>document.querySelectorAll('details').forEach(d=>d.open=false);filter();
</script></body></html>'''
page = page.replace('SHA', sha).replace('SUMMARY', summary).replace('ANONYMOUS', str(data['anonymous_callbacks'])).replace('OPTIONS', ''.join(f'<option>{e(c)}</option>' for c in categories)).replace('SECTIONS', ''.join(sections))
(out / 'function-inventory.html').write_text(page)
(out / 'function-inventory.json').write_text(json.dumps({'revision':sha,'categories':dict(counts),'anonymous_R_callbacks_excluded':data['anonymous_callbacks'],'functions':rows},indent=2))
print(json.dumps(dict(counts), indent=2))
print('Total named definitions:', len(rows))

# Authored coverage and evidence stay separate from generated source inventories.
coverage = json.loads((book / 'coverage.json').read_text())
assignments = coverage['functions']
expected = exports | methods
if set(assignments) != expected:
    raise ValueError(f'Coverage differs from NAMESPACE: missing={expected-set(assignments)}, obsolete={set(assignments)-expected}')
chapters = {c['id']: c for c in coverage['chapters']}
labels = {}
for path in book.glob('*.qmd'):
    for label in re.findall(r'^#\| label: ([\w-]+)', path.read_text(), re.M):
        if label in labels:
            raise ValueError(f'Duplicate chunk label: {label}')
        labels[label] = path.name
for name, item in assignments.items():
    if item['chapter'] not in chapters:
        raise ValueError(f'Unknown chapter for {name}')
    for label in item['examples']:
        if label not in labels:
            raise ValueError(f'Missing example {label} for {name}')
for option in coverage['method_options']:
    if option['function'] not in expected:
        raise ValueError(f'Unknown public method: {option}')

generated = book / '_generated'
generated.mkdir(exist_ok=True)
def write_qmd(name, lines):
    (generated / name).write_text('\n'.join(lines) + '\n')

public = ['| Function or method | Planned chapter | Authored example |', '|---|---|---|']
for name, item in sorted(assignments.items()):
    chapter = chapters[item['chapter']]['title']
    examples = ', '.join(f'[{label}]({labels[label]}#{label})' for label in item['examples']) or 'Planned'
    public.append(f'| `{name}()` | {chapter} | {examples} |')
public += ['', '## Method option coverage', '',
           'This initial checklist is seeded from selection, prediction, mating, and architecture options. It will expand with each methods chapter; it is not yet an exhaustive option inventory.', '',
           '| Function | Method or configuration | Example status |', '|---|---|---|']
for option in coverage['method_options']:
    public.append(f'| `{option["function"]}()` | `{option["option"]}` | {option["status"]} |')
write_qmd('public-index.qmd', public)

developer = []
for category in categories[2:]:
    developer += [f'## {category}', '', descriptions[category], '',
                  '| Definition | Enclosing scope | Source |', '|---|---|---|']
    for row in sorted((r for r in rows if r['category'] == category), key=lambda r: (r['file'], r['line'])):
        link = f'https://github.com/samuelbfernandes/simplePHENOTYPES/blob/{sha}/{row["file"]}#L{row["line"]}'
        developer.append(f'| `{row["name"]}` | {row["scope"] or "—"} | [{row["file"]}:{row["line"]}]({link}) |')
    developer.append('')
write_qmd('developer-index.qmd', developer)

plan = ['| Chapter | Scope | Draft status |', '|---|---|---|']
for chapter in coverage['chapters']:
    plan.append(f'| {chapter["title"]} | {chapter["scope"]} | {chapter["status"]} |')
write_qmd('chapter-plan.qmd', plan)

equations = json.loads((book / 'evidence.json').read_text())
review_file = book / 'review-status.json'
review_current = False
if review_file.exists():
    review = json.loads(review_file.read_text())
    source_files = sorted(list((repo / 'R').glob('*.R')) + list((repo / 'src/rust/src').glob('*.rs')) + [repo / 'NAMESPACE', repo / 'DESCRIPTION'])
    source_hash = hashlib.sha256()
    for path in source_files:
        source_hash.update(str(path.relative_to(repo)).encode() + b'\0' + path.read_bytes() + b'\0')
    review_current = source_hash.hexdigest() == review.get('package_source_sha256') and all(
        (repo / path).exists() and hashlib.sha256((repo / path).read_bytes()).hexdigest() == expected_hash
        for path, expected_hash in review['files'].items())
if not review_current:
    print('Scientific review is pending or stale for the current source and chapter contents.')
evidence = ['| Claim | Evidence type and provenance | Executable check | Independent review |', '|---|---|---|']
for item in equations:
    refs = []
    for fn in item['functions']:
        matches = [r for r in rows if r['name'] == fn and not r['nested']]
        if len(matches) > 1:
            # An R binding and its Rust kernel share a name; cite the kernel.
            matches = [r for r in matches if r['category'] == 'Rust implementation'] or matches
        if len(matches) != 1:
            raise ValueError(f'Ambiguous evidence function: {fn}')
        r = matches[0]
        link = f'https://github.com/samuelbfernandes/simplePHENOTYPES/blob/{sha}/{r["file"]}#L{r["line"]}'
        refs.append(f'[`{fn}`]({link})')
    checks = []
    for label in item['checks']:
        if label not in labels:
            raise ValueError(f'Missing evidence check: {label}')
        checks.append(f'[{label}]({labels[label]}#{label})')
    provenance = ', '.join(refs) or item.get('provenance', 'derivation in Appendix A')
    review_label = item['review'] if review_current else 'Pending or stale: content differs from the recorded review'
    evidence.append(f'| {item["claim"]} | {item["kind"]}; {provenance} | {", ".join(checks) or item.get("check_note", "Not yet recorded")} | {review_label} |')
write_qmd('evidence.qmd', evidence)

historical = (repo / 'docs/simplePHENOTYPES_equation_code_map.md').read_text()
start = historical.index('# Appendix B')
end = historical.find('# Appendix C', start)
index = historical[start:end if end >= 0 else len(historical)]
table = [line for line in index.splitlines() if line.lstrip().startswith('|')]
if not table:
    raise ValueError('Could not extract historical equation index')
write_qmd('equation-map-index.qmd', table)
import shutil
shutil.copyfile(repo / 'docs/simplePHENOTYPES_equation_code_map.pdf', out / 'equation-code-map.pdf')
version = re.search(r'^Version:\s*(.*)$', (repo / 'DESCRIPTION').read_text(), re.M).group(1)
dirty = subprocess.check_output(['git', '-c', 'core.fileMode=false', 'status', '--porcelain', '--', 'R', 'src/rust/src', 'NAMESPACE'], text=True).strip()
record = [f'- Package version: `{version}`.', f'- Source revision: `{sha}`.',
          f'- Package source differs from that revision: **{"yes; source links describe HEAD, not local changes" if dirty else "no"}**.',
          f'- Public coverage assignments: **{len(assignments)}**.',
          f'- Public functions/methods with authored examples: **{sum(bool(v["examples"]) for v in assignments.values())}**.',
          f'- Pilot scientific review matches the current source and reviewed files: **{"yes" if review_current else "no"}**.',
          '- This count records authored chunks, not completed independent scientific review. A successful fresh render is required to establish execution.']
write_qmd('build-record.qmd', record)
print('Coverage and source-index checks passed.')
