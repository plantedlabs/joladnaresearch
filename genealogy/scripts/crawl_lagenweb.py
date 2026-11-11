import re, html, subprocess, sys, time, os, urllib.parse, json
out = sys.argv[1]; start = 'https://lagenweb.org/catahoula/index.html'; root = 'https://lagenweb.org/catahoula/'
seen, queue, pages = set(), [start], {}
while queue and len(pages) < 400:
    u = queue.pop(0)
    if u in seen: continue
    seen.add(u)
    p = subprocess.run(['curl', '-sS', '-L', '--max-time', '30', u], capture_output=True)
    if p.returncode: continue
    t = p.stdout.decode('utf8', 'replace'); pages[u] = t
    fn = os.path.join(out, re.sub(r'[^A-Za-z0-9._-]', '_', u[len(root):]) or 'index')
    open(fn, 'w').write(t)
    if u.endswith(('.html', '.htm', '/')):
        for h in re.findall(r'href="([^"#]+)"', t, flags=re.I):
            v = urllib.parse.urljoin(u, h)
            if v.startswith(root) and v.endswith(('.html', '.htm', '.txt', '/')) and v not in seen:
                queue.append(v)
    time.sleep(0.4)
print(len(pages), 'pages fetched')
