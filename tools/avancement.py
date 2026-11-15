"""Progress of PLAN.md: recompute the table and the list of free items (A.4.2).

python -m tools.avancement            # print the checks and the counts
python -m tools.avancement --write    # rewrite the progress table in PLAN.md
"""

import io
import re
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
P = Path(__file__).resolve().parents[2] / "PLAN.md"
txt = open(P, encoding="utf-8").read()
body_start = txt.index("## A. Environnement")
head, body = txt[:body_start], txt[body_start:]

blocks = {}  # letter -> title
items = {}  # "A.1.1" -> dict
cur_block = cur_pt = None
for line in body.splitlines():
    m = re.match(r"## ([A-L])\. (.+)", line)
    if m:
        cur_block = m.group(1)
        blocks[cur_block] = m.group(2)
        continue
    m = re.match(r"### (\d+)\. ", line)
    if m:
        cur_pt = m.group(1)
        continue
    m = re.match(r"- \[( |x)\] \*\*(\d+)\.(\d+)\*\* (.*)", line)
    if m:
        key = f"{cur_block}.{m.group(2)}.{m.group(3)}"
        assert m.group(2) == cur_pt, (key, cur_pt)
        text = m.group(4)
        tags = re.findall(r"`([^`]*)`", text)
        dep = [t for t in tags if t == "libre" or t.startswith(("après", "progressif", "attend"))]
        assert len(dep) == 1, (key, dep)
        items[key] = dict(done=m.group(1) == "x", dep=dep[0], block=cur_block)

ref_re = re.compile(r"\b([A-L])(?:\.(\d+))?(?:\.(\d+))?\b")


def expand(ref):
    """return list of item keys covered by a reference like D, D.2, D.2.3"""
    parts = ref.split(".")
    return [k for k in items if k.split(".")[: len(parts)] == parts]


errors = []
graph = {}
for k, it in items.items():
    d = it["dep"]
    if d == "libre":
        graph[k] = []
        continue
    if d.startswith("attend"):
        graph[k] = []
        continue
    refs = re.findall(
        r"\b[A-L](?:\.\d+){0,2}\b",
        d.split(":", 1)[-1] if d.startswith("progressif") else d[len("après") :],
    )
    deps = []
    for r in refs:
        e = expand(r)
        if not e:
            errors.append(f"{k}: référence inconnue {r}")
        deps += e
    if k in deps:
        errors.append(f"{k}: dépend de lui-même")
    graph[k] = deps if d.startswith("après") else []  # progressif: pas de blocage strict

# cycles (après seulement)
state = {}


def dfs(n, path):
    state[n] = 1
    for m in graph[n]:
        if state.get(m) == 1:
            errors.append("cycle: " + " -> ".join(path + [n, m]))
        elif m not in state:
            dfs(m, path + [n])
    state[n] = 2


sys.setrecursionlimit(10000)
for k in items:
    if k not in state:
        dfs(k, [])

# références hors étiquettes (texte du plan)
all_refs = set(re.findall(r"\b([A-L]\.\d+(?:\.\d+)?)\b", txt))
for r in sorted(all_refs):
    if not expand(r):
        errors.append(f"référence au texte sans cible: {r}")

# dépendances inverses
rev = {}
for k, ds in graph.items():
    for d in ds:
        rev.setdefault(d, set()).add(k)
libre = [k for k, it in items.items() if it["dep"] == "libre" and not it["done"]]
attend = [k for k, it in items.items() if it["dep"].startswith("attend")]


def sortkey(k):
    b, p, s = k.split(".")
    return (b, int(p), int(s))


print("items:", len(items))
print("libre:", " · ".join(sorted(libre, key=sortkey)))
print("attend:", sorted(attend, key=sortkey))
top = sorted(rev.items(), key=lambda kv: -len(kv[1]))[:12]
print("top débloque:", [(k, len(v)) for k, v in top])
print("ERREURS:", errors if errors else "aucune")

# blocage transitif par attend
blocked = set()
changed = True
for a in attend:
    blocked.add(a)
while changed:
    changed = False
    for k, ds in graph.items():
        if k not in blocked and any(d in blocked for d in ds):
            blocked.add(k)
            changed = True
print("bloqués transitivement par un attend:", len(blocked))

# tableau d'avancement
bar = lambda r: "█" * round(r * 10) + "░" * (10 - round(r * 10))
rows = []
T = D = 0
W = 0
for b, title in blocks.items():
    ks = [k for k in items if k.startswith(b + ".")]
    t = len(ks)
    d = sum(items[k]["done"] for k in ks)
    w = sum(1 for k in ks if not items[k]["done"] and items[k]["dep"].startswith("attend"))
    T += t
    D += d
    W += w
    rows.append(f"| {b}. {title} | {d} | {t} | `{bar(d / t)}` {round(100 * d / t)} % | {w} |")
tbl = ["| Bloc | Fait | Total | Avancement | En attente |", "|---|---:|---:|---|---:|"] + rows
tbl.append(
    f"| **Ensemble du plan** | **{D}** | **{T}** | `{bar(D / T)}` **{round(100 * D / T)} %** | **{W}** |"
)
tbl.append(
    f"| Hors attentes extérieures | {D} | {T - W} | `{bar(D / (T - W))}` {round(100 * D / (T - W))} % | — |"
)
tbl = "\n".join(tbl)

if "--write" in sys.argv:
    new = re.sub(
        r"\| Bloc \| Fait.*?\| Hors attentes extérieures[^\n]*", lambda m: tbl, txt, flags=re.S
    )
    new = re.sub(
        r"(### Départs immédiats \(`libre`\)\n)[^\n]*",
        lambda m: m.group(1) + " · ".join(sorted(libre, key=sortkey)),
        new,
    )
    open(P, "w", encoding="utf-8", newline="\n").write(new)
    print("écrit")

if "--check" in sys.argv:
    stale = tbl not in txt
    if errors or stale:
        print("PLAN.md: références invalides" if errors else "PLAN.md: tableau d'avancement périmé")
        sys.exit(1)
    print("PLAN.md: références valides et tableau à jour")
