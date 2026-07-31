#!/usr/bin/env python3
"""
dataio.py -- chargement ROBUSTE et BRUYANT des fichiers .dat du solveur.

POURQUOI CE MODULE EXISTE (2026-07-22)
--------------------------------------
`postprocess.py` et `analyse_weber_ensemble.py` avaient chacun leur propre
chargeur, avec la meme faille :

    if ncol is None: ncol = len(v)      # largeur fixee par la 1re ligne
    if len(v) == ncol: rows.append(v)   # les autres largeurs sont JETEES

... silencieusement. Or `event bubble` a ete corrige le 2026-07-22 et ecrit
desormais **16** colonnes au lieu de 13 (ajout de n_regions, V_tag, straddle).
Donc des qu'on reprend un run existant avec le nouveau binaire, `bubble.dat`
contient des lignes a 13 puis a 16 colonnes, et l'ancien chargeur aurait ignore
**toute la partie posterieure a la reprise** sans le moindre message.

C'est exactement la classe de bug qui a coute deux semaines a ce projet
(cf. la "perte de volume" fantome et la "fausse fragmentation") : une mesure
fausse qui ne previent pas. La regle retenue :

    UN CHARGEUR QUI JETTE DES DONNEES DOIT LE DIRE, FORT.

Comportement
------------
- les lignes vides et les commentaires (#) sont ignores, normalement ;
- on recense les largeurs de ligne presentes. Une largeur vue >= 2 fois est un
  "format" ; une largeur vue une seule fois est presumee etre une ligne tronquee
  par un run tue en plein ecriture -> ecartee ;
- s'il y a plusieurs formats (reprise avec un binaire plus recent), on TRONQUE
  toutes les lignes a la largeur commune minimale. C'est sur ici parce que les
  colonnes 1-13 de bubble.dat ont garde leur ordre ET leur signification ;
- tout ecart (plusieurs formats, lignes ecartees) est signale sur stderr.
"""
import os
import sys
from collections import Counter

import numpy as np


def load_table(path, verbose=True):
    """Charge un .dat en tableau numpy. Renvoie None si absent/vide.

    Tronque a la largeur commune si le fichier melange plusieurs formats
    (typiquement : run repris avec un binaire ecrivant plus de colonnes),
    et previent sur stderr a chaque fois qu'il ecarte ou tronque quelque chose.
    """
    if not os.path.isfile(path):
        return None

    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            continue
    if not rows:
        return None

    widths = Counter(len(r) for r in rows)
    # une largeur vue une seule fois = ligne tronquee par un kill en plein write
    formats = sorted(w for w, c in widths.items() if c >= 2)
    if not formats:                      # fichier minuscule : on prend tel quel
        formats = sorted(widths)

    ncommon = min(formats)
    kept = [r[:ncommon] for r in rows if len(r) >= ncommon]
    dropped = len(rows) - len(kept)

    if verbose:
        name = os.path.relpath(path)
        if len(formats) > 1:
            detail = ", ".join(f"{w} col x{widths[w]}" for w in formats)
            print(f"[dataio] {name} : PLUSIEURS FORMATS ({detail}) -> tronque a "
                  f"{ncommon} colonnes. Reprise avec un binaire different ?",
                  file=sys.stderr)
        if dropped:
            print(f"[dataio] {name} : {dropped} ligne(s) ecartee(s) (trop "
                  f"courtes, ecriture interrompue ?)", file=sys.stderr)

    return np.array(kept) if kept else None


if __name__ == "__main__":
    for p in sys.argv[1:]:
        a = load_table(p)
        print(f"{p} -> {None if a is None else a.shape}")
