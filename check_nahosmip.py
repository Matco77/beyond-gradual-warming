#!/usr/bin/env python3
"""
Ispeziona i tar NAHosMIP (Zenodo 14884357) SENZA estrarli e riporta,
per ogni modello e variabile, quali esperimenti contengono.

Uso:
    python3 check_nahosmip.py /percorso/ai/tar
    python3 check_nahosmip.py .

Servono sia u03-hos (hosing) sia picon (controllo): il delta si calcola
tra i due, quindi un modello con solo uno dei due e' inutilizzabile.
"""

import sys, os, re, tarfile
from collections import defaultdict

RICHIESTI = ("u03-hos", "picon")

# Ordine importante: le stringhe piu' lunghe vanno confrontate per prime
ESPERIMENTI_NOTI = [
    "u03-r100", "u03-r70", "u03-r50", "u03-r20", "u03-hos",
    "g01-r100", "g01-r50", "g01-hos",
    "u05-r100", "u05-r50", "u05-hos",
    "g03-r100", "g03-r50", "g03-hos",
    "picon", "piControl", "control",
]

RE_ANNI = re.compile(r"(\d{4})\d{2}-(\d{4})\d{2}\.nc$")


def analizza(nome):
    """
    Estrae variabile / modello / esperimento da un nome CMIP tipo
        tas_Amon_IPSL-CM6A-LR_u03-hos_r1i1p1f1_gn_215001-215912.nc

    NON si fida della posizione dei campi: nel record alcuni tar
    omettono il variant (CESM2, MPI-*) e CanESM5 usa r1i1p2f1.
    L'esperimento viene trovato per corrispondenza con la lista nota.
    """
    base = os.path.basename(nome)
    if not base.endswith(".nc"):
        return None
    campi = base[:-3].split("_")
    if len(campi) < 3:
        return None

    esperimento, idx = None, None
    for i, c in enumerate(campi):
        if c in ESPERIMENTI_NOTI:
            esperimento, idx = c, i
            break

    modello = "_".join(campi[2:idx]) if (idx and idx > 2) else campi[2]

    anni = RE_ANNI.search(base)
    a_i, a_f = (int(anni.group(1)), int(anni.group(2))) if anni else (None, None)

    return dict(variabile=campi[0], modello=modello,
                esperimento=esperimento or "SCONOSCIUTO",
                a_i=a_i, a_f=a_f, nome=base)


def main():
    radice = sys.argv[1] if len(sys.argv) > 1 else "."
    tars = sorted(f for f in os.listdir(radice) if f.endswith((".tar", ".tar.gz")))
    if not tars:
        print(f"Nessun .tar in {os.path.abspath(radice)}")
        return 1

    inv = defaultdict(lambda: defaultdict(list))   # (modello, var) -> exp -> [(a_i,a_f)]
    ignoti = []

    for t in tars:
        p = os.path.join(radice, t)
        try:
            with tarfile.open(p, "r:*") as tf:
                membri = [m.name for m in tf.getmembers() if m.isfile()]
        except tarfile.TarError as e:
            print(f"!! {t}: illeggibile ({e})")
            continue

        info_utili = [i for i in (analizza(m) for m in membri) if i]
        print(f"{t}: {len(membri)} membri, {len(info_utili)} file .nc riconosciuti")
        for i in info_utili:
            inv[(i["modello"], i["variabile"])][i["esperimento"]].append((i["a_i"], i["a_f"]))
            if i["esperimento"] == "SCONOSCIUTO":
                ignoti.append(i["nome"])

    print("\n" + "=" * 76)
    print("DETTAGLIO PER MODELLO E VARIABILE")
    print("=" * 76)

    for (mod, var) in sorted(inv):
        exps = inv[(mod, var)]
        print(f"\n{mod}  [{var}]")
        for e in sorted(exps):
            ints = [x for x in exps[e] if x[0] is not None]
            if ints:
                p0, p1 = min(i[0] for i in ints), max(i[1] for i in ints)
                anni = f"{p0}-{p1} ({p1 - p0 + 1} anni)"
            else:
                anni = "anni non leggibili"
            print(f"    {e:<12} {len(exps[e]):>4} file   {anni}"
                  f"{'   <-- serve' if e in RICHIESTI else ''}")
        mancanti = [e for e in RICHIESTI if e not in exps]
        print(f"    >> MANCA: {', '.join(mancanti)}" if mancanti
              else "    >> OK, u03-hos e picon presenti")

    if ignoti:
        print("\nEsperimento non riconosciuto in questi file:")
        for n in ignoti[:10]:
            print(f"    {n}")
        if len(ignoti) > 10:
            print(f"    ... e altri {len(ignoti) - 10}")

    print("\n" + "=" * 76)
    print("RIEPILOGO")
    print("=" * 76)
    per_mod = defaultdict(dict)
    for (mod, var), exps in inv.items():
        per_mod[mod][var] = all(e in exps for e in RICHIESTI)
    for mod in sorted(per_mod):
        ok = sorted(v for v, b in per_mod[mod].items() if b)
        ko = sorted(v for v, b in per_mod[mod].items() if not b)
        print(f"  {mod:<22} usabile: {', '.join(ok) if ok else '(nessuna)'}"
              + (f"   |  incompleto: {', '.join(ko)}" if ko else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
