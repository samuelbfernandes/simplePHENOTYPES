"""Check local files and fragment targets in a rendered book."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

root = Path(__file__).resolve().parents[2] / 'docs/book/_book'

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.hrefs.append(attrs['href'])

pages = {}
for path in root.rglob('*.html'):
    page = Links()
    page.feed(path.read_text())
    pages[path.resolve()] = page
if not pages:
    raise SystemExit('No rendered HTML pages found; run quarto render docs/book first.')

broken = []
for path, page in pages.items():
    for href in page.hrefs:
        url = urlsplit(href)
        if url.scheme or url.netloc or not (url.path or url.fragment):
            continue
        target = (path.parent / unquote(url.path)).resolve() if url.path else path
        if target.is_dir():
            target /= 'index.html'
        if not target.exists():
            broken.append(f'{path.name}: missing file {href}')
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            broken.append(f'{path.name}: missing anchor {href}')
if broken:
    raise SystemExit('\n'.join(broken))
print(f'Checked {len(pages)} HTML pages: no broken local links or fragment targets.')
