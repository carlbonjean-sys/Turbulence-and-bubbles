#!/usr/bin/env python3
"""Découpe sec9.tex (journal) en 3 sections du rapport final, par TITRE de
sous-section (robuste aux décalages de lignes). Retire les 2 sous-sections
redondantes avec le nouveau cadre/méthodes."""
import io, re
from pathlib import Path

SRC = Path("../sections/sec9.tex")
DST = Path(".")
txt = io.open(SRC, encoding="utf-8").read()

# découpe en blocs à chaque \subsection (le préambule de section est ignoré)
parts = re.split(r"(?=^\\subsection\{)", txt, flags=re.M)

def title(block):
    m = re.match(r"\\subsection\{(.*?)(\}|\n)", block, flags=re.S)
    return (m.group(1) if m else "").strip()

# clé = fragment distinctif du titre -> destination
DROP = ["Cadre adimensionnel", "Protocole numérique"]         # -> sec2/sec3
RESULT = ["Perspective", "Premiers résultats", "Campagne en cours",
          "Confrontation à la littérature"]
MECA = ["Identification du mécanisme"]

valid, result, meca = [], [], []
for b in parts:
    if not b.startswith("\\subsection"):
        continue                                   # entête de section, ignoré
    t = title(b)
    if any(k in t for k in DROP):
        continue
    if any(k in t for k in MECA):
        meca.append(b)
    elif any(k in t for k in RESULT):
        result.append(b)
    else:
        valid.append(b)

HDR = ("% genere par split_sec9.py -- contenu issu de sec9.tex (audite), "
       "reorganise pour le rapport final.\n")

def write(fname, section_title, blocks, intro=""):
    with io.open(DST / fname, "w", encoding="utf-8") as f:
        f.write(HDR)
        f.write("\\section{%s}\n\n" % section_title)
        if intro:
            f.write(intro + "\n\n")
        f.write("\n".join(blocks))
    n = sum(b.count("\\subsection{") for b in blocks)
    print(f"  {fname}: {n} sous-sections, {sum(len(b) for b in blocks)} car.")

write("sec4_validation.tex", "Validation du modèle numérique", valid,
      intro="Cette section établit que la chaîne de calcul est digne de confiance : "
            "résolution des échelles, découverte et correction d'un artefact de gravité, "
            "conservation du volume, et convergence en maillage en régime laminaire "
            "comme turbulent.")
write("sec5_resultats.tex", "Résultats : effet de la turbulence sur la vitesse d'ascension",
      result)
write("sec6_mecanisme.tex", "Mécanisme du ralentissement", meca)
print("OK")
