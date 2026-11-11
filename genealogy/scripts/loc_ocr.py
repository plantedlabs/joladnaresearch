import json, subprocess, re, sys, time
files = sys.argv[1:]
pat = re.compile(r'(dr[ou]m-?\s?g[oa]{1,2}l[ae]?|drumgold|brungood|woodley)', re.I)
seen = set()
for f in files:
    try: d = json.load(open(f))
    except Exception: continue
    for r in d.get('results', [])[:25]:
        u = r.get('id')
        if not u or u in seen: continue
        seen.add(u)
        j = subprocess.run(['curl', '-sS', '--max-time', '60', u.replace('http://', 'https://') + '&fo=json'], capture_output=True).stdout
        try: jj = json.loads(j)
        except Exception: print('noparse', u); continue
        ft = None
        def walk(o):
            global ft
            if isinstance(o, dict):
                for k, v in o.items():
                    if k in ('fulltext_file',) and isinstance(v, str): ft = v
                    walk(v)
            elif isinstance(o, list):
                for v in o: walk(v)
        walk(jj)
        if not ft: print('no fulltext', u); continue
        txt = subprocess.run(['curl', '-sS', '--max-time', '60', ft], capture_output=True).stdout.decode('utf8', 'replace')
        words = re.findall(r'CONTENT="([^"]*)"', txt)
        txt = ' '.join(words) if words else re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', txt))
        if not txt.strip(): print('EMPTY OCR', u, ft); continue
        hits = [txt[max(0, m.start() - 300):m.end() + 300] for m in pat.finditer(txt)]
        t = r.get('partof_title'); t = t[0] if isinstance(t, list) and t else ''
        print(f"== {r.get('date')} | {t[:60]} | {u} | exact-variant hits: {len(hits)}")
        for h in hits[:2]: print('   >>', h)
        time.sleep(0.3)
