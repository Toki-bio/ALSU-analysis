import os,re,sys,collections,urllib.parse
root="C:/work/ALSU-analysis"
pages=[]
for dp,dn,fn in os.walk(root):
    if ".git" in dp.split(os.sep): continue
    for f in fn:
        if f.endswith((".html",".htm")): pages.append(os.path.relpath(os.path.join(dp,f),root).replace("\\","/"))
links=collections.defaultdict(set)
for p in pages:
    t=open(os.path.join(root,p),encoding="utf-8",errors="ignore").read()
    for m in re.finditer(r'(?:href|src)\s*=\s*["\']([^"\'#?]+)',t):
        h=m.group(1)
        if re.match(r"(https?:|mailto:|data:|javascript:)",h): continue
        tgt=os.path.normpath(os.path.join(os.path.dirname(p),urllib.parse.unquote(h))).replace("\\","/")
        links[p].add(tgt)
seen={"index.html"}; q=["index.html"]
while q:
    p=q.pop()
    for t in links.get(p,()):
        if t not in seen: seen.add(t); q.append(t)
orph=[p for p in pages if p not in seen]
inbound=collections.Counter(t for p in links for t in links[p] if t!=p)
print("pages",len(pages),"reachable from index",len([p for p in pages if p in seen]))
print("ORPHANS (not reachable from index.html):")
for p in sorted(orph): print(" ",p,"| inbound from any page:",inbound.get(p,0), "| size",os.path.getsize(os.path.join(root,p)))
broken=[(p,t) for p in pages for t in links[p] if t.endswith((".html",".md",".py",".json",".tsv",".png")) and not os.path.exists(os.path.join(root,t))]
print("BROKEN internal links:",len(broken))
for b in broken[:60]: print(" ",b)
print("\nBROKEN in live pages (not archive/old_logs/daily):")
live=[b for b in broken if not b[0].startswith(("steps/archive","old_logs"))]
print(len(live))
import collections as C
for p,t in live: print(" ",p,"->",t)
