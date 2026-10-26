# מדריך צוות — בניית ייחוס לאורגניזם והרצת הפייפליין

> מדריך תפעולי בעברית: מה הריפו עושה, מה מריצים, ואיך מוסיפים אורגניזם חדש מקצה לקצה.
>
> **בנאי:** `multiomics-annotation-prep` · **פייפליין:** `multiomic-core` · **מקור נתונים:** KEGG REST

---

## איזה Action אני צריכה?

לא צריך להריץ קוד. הכל נעשה דרך **GitHub Actions**: ממלאים טופס קצר, מחשב של GitHub מוריד את נתוני הייחוס ובונה את הקבצים, והתוצאה מתפרסמת כ‑**Release** — דף הורדה קבוע ומתוארך — שאליו מפנים את הפייפליין.

| מה אני מנתחת | האורגניזם | ה‑Action |
|---|---|---|
| **חלבונים** | כל אורגניזם, כשצריך GO / KEGG לפי מזהי החלבונים שלי (חובה לאורגניזם לא‑מודלי) | [1. Build gene annotation](#1-build-gene-annotation-rna--proteomics) |
| **RNA‑seq** | לא‑מודלי (טרנסקריפטום de novo, או אורגניזם מאתר כמו GiardiaDB) | [1. Build gene annotation](#1-build-gene-annotation-rna--proteomics) |
| **RNA‑seq** | מודלי, ב‑Ensembl או ב‑UniProt | [2. Build model-organism annotation](#2-build-model-organism-annotation) |
| **מטבולומיקה** | כל אורגניזם ש‑KEGG מכיר (או קרוב שלו) | [3. Build & publish organism artifacts](#3-build--publish-organism-artifacts) |
| — | פרסמתי משהו בטעות | [4. Delete releases](#4-delete-releases) |

> **חלבונים של אורגניזם מודלי** (אדם, עכבר…) בדרך כלל לא צריכים כלום מהריפו הזה: `multiomics-core` טוען GO / KEGG לבד לפי `annotation.organism`.

**איך מריצים כל אחד מהם:** בריפו → לשונית **Actions** → בוחרים את ה‑Action משמאל → **Run workflow** (למעלה מימין) → ממלאים → **Run workflow**. כשהריצה ירוקה, הקבצים נמצאים בדף **Releases** של הריפו.

---

## 1. Build gene annotation (RNA / proteomics)

**מה הוא עושה:** הופך את האנוטציה של הפרויקט שלך לקבוצות GO ו‑KEGG **לפי המזהים שבנתונים שלך** (מזהי גנים ב‑RNA, מזהי חלבונים בפרוטאומיקה).

### שלב 1 — מעלים את הקלט לריפו

יוצרים תיקייה `annotation_inputs/<project>/`. ב‑GitHub: נכנסים ל‑`annotation_inputs` → **Add file → Upload files** → **גוררים תיקייה** בשם הפרויקט (גרירת תיקייה יוצרת אותה). שמים בה את מה שיש:

| קובץ | מה זה | נותן |
|---|---|---|
| `query.ko.txt` | תוצאת KAAS: `מזהה ⇥ KO` בכל שורה | KEGG |
| `trinotate_go.txt` | טבלת GO: `מזהה ⇥ מזהי GO` (מ‑Trinotate, eggNOG…) | GO |
| קובץ `.gaf` | קובץ GO מאתר אורגניזם (למשל GiardiaDB). **לא משנים את שם הקובץ** — יש בו את גרסת המאגר, והיא נרשמת ב‑Release. | GO |
| `id_map.txt` | *אופציונלי*, רק אם צריך לתרגם מזהים: `מזהה_בקובץ ⇥ מזהה_בנתונים` | — |

כללים: `trinotate_go.txt` **או** `.gaf`, לא שניהם; `query.ko.txt` **או** השדה `kegg_org`, לא שניהם. אם כל מה שצריך הוא KEGG מ‑`kegg_org`, אפשר בכלל בלי תיקייה.

### שלב 2 — ממלאים את הטופס

| שדה | מה ממלאים |
|---|---|
| `project` | שם התיקייה משלב 1 |
| `omics` | `rna` או `proteomics` |
| `kegg_org` | *אופציונלי.* קוד אורגניזם ב‑KEGG (למשל `gla` ל‑Giardia, `cre` ל‑Chlamydomonas). בונה KEGG ישירות מהאנוטציה של KEGG — **בלי KAAS**. לבדיקה שהקוד קיים: `https://rest.kegg.jp/list/<code>`. |
| `ncbi_assembly` | *אופציונלי.* מספר assembly ב‑NCBI (למשל `GCA_000002435.2`). משתמשים בו כשבנתונים יש מזהי חלבון של NCBI (כמו `KAE8301220.1`) וב‑GAF / KEGG יש שמות גנים (כמו `GL50803_0017558`): ה‑Action מוריד את הטבלה של NCBI שמקשרת ביניהם ומתרגם לבד. |
| `id_label` | *אופציונלי.* מילה קצרה לסוג המזהים, נכנסת לשם ה‑Release (למשל `kae`). מאפשרת שתי בניות של אותו פרויקט באותו יום. ברירת מחדל `ncbi` כשממלאים `ncbi_assembly`. |
| `date` | משאירים ריק (היום) |
| `strip_isoform` | רק ל‑RNA: מסומן אם ה‑counts לפי גן (`..._g1`), לא מסומן אם לפי טרנסקריפט (`..._g1_i1`) |

### שלב 3 — מה מקבלים

Release בשם `prot_annot_<project>[_<label>]_<date>` (חלבונים) או `rna_annot_<project>[_<label>]_<date>` (RNA), ובו:

- `GO_BP.gmt`, `GO_MF.gmt`, `GO_CC.gmt`, `KEGG_pathway.gmt` — קבוצות הגנים
- `GO2gene_*.tab`, `GO2name_*.tab`, `KEGG_pathway2gene.tab`, `KEGG_pathway2name.tab` — אותו מידע כטבלאות
- `GO_provenance.json` / `KEGG_provenance.json` — מאיפה הגיע כל קלט (שמות קבצים, גרסאות, checksums)
- `SHA256SUMS.txt` ו‑zip של הכל

בהערות ה‑Release רשום גם איזה commit ואיזו ריצה בנו אותו.

### שלב 4 — משתמשים ב‑`multiomics-core`

**חלבונים** קוראים את קבצי ה‑`.gmt`. מורידים אותם ל‑`annotation/` בתוך **תיקיית ה‑raw** של הפרויקט (`paths.raw`; הנתיבים היחסיים למטה יחסיים אליה):

```yaml
modes:
  proteomics:
    annotation:
      organism: "Giardia lamblia"
      skip_annotation: true        # ה-GMT לפי מזהי Protein.Group כמו שהם
    pathway:
      gmt_file:
        - annotation/GO_BP.gmt
        - annotation/GO_MF.gmt
        - annotation/GO_CC.gmt
        - annotation/KEGG_pathway.gmt
```

**RNA** קורא את קבצי ה‑`.tab`. פורסים את ה‑zip לתיקייה בפרויקט ומפנים אליה (נתיב יחסי כאן יחסי ל**תיקיית הפרויקט**, `project.dir`):

```yaml
modes:
  rna:
    enrichment:
      annotation_dir: "data/annotation"   # התיקייה עם קבצי ה-.tab
```

> ⚠️ **המזהים חייבים להתאים לנתונים.** אם ב‑GMT יש `KAE…` ובנתונים `GL50803_…`, שום דבר לא יתאים וה‑enrichment לא ימצא כלום. אם יש ספק — פותחים GMT, מסתכלים על כמה מזהים, ומשווים לעמודת Protein.Group / gene_id.

**דוגמה — חלבונים של Giardia:** ה‑GAF מ‑GiardiaDB הועלה ל‑`annotation_inputs/giardia/`, וה‑Action הורץ עם `omics: proteomics`, `kegg_org: gla`, `ncbi_assembly: GCA_000002435.2`. התוצאה: GO + KEGG לפי מזהי החלבונים `KAE…` שיש בנתוני הפרוטאומיקה.

---

## 2. Build model-organism annotation

**מה הוא עושה:** לאורגניזם עם גנום ב‑**Ensembl** או פרוטאום ב‑**UniProt** — מוריד GO (ואופציונלית KEGG) ישירות, בלי קבצי קלט.

| שדה | מה ממלאים |
|---|---|
| `source` | `ensembl` או `uniprot` |
| `label` | שם קצר ל‑Release (למשל `athaliana`) |
| `dataset`, `mart`, `host` | רק Ensembl. ה‑dataset ב‑BioMart (למשל `athaliana_eg_gene`); לצמחים/אצות/פרוטיסטים גם `mart` מתאים (למשל `plants_mart`) ו‑`host` (למשל `https://plants.ensembl.org`) |
| `taxon_id` | רק UniProt. מספר ה‑taxon ב‑NCBI (למשל `3702`) |
| `reviewed` | רק UniProt. מסומן = רק רשומות Swiss-Prot |
| `kegg_org` | *אופציונלי.* קוד KEGG כדי לבנות גם טבלאות KEGG |
| `date` | משאירים ריק |

**מקבלים:** Release בשם `<label>_<source>_annot_<date>` עם **טבלאות `.tab`** של `GO2gene_*` / `GO2name_*` (ו‑KEGG).

**משתמשים ב:** פייפליין ה‑**RNA** — פורסים ומגדירים `modes.rna.enrichment.annotation_dir` (כמו ב‑Action 1). המזהים הם Ensembl gene ids או UniProt accessions, ולכן הם חייבים להתאים ל‑counts.

---

## 3. Build & publish organism artifacts

**מה הוא עושה:** בונה קבצי מסלולים למטבולומיקה מ‑KEGG, לקוד אורגניזם אחד או יותר. לאורגניזם ש‑KEGG לא מכיר (למשל *Coelastrella*) משתמשים בקרובים ורושמים את האורגניזם האמיתי כ‑target. הפרטים המלאים — בסעיפים "התמונה הגדולה" ו"מה מריצים" למטה.

| שדה | מה ממלאים |
|---|---|
| `organisms` | קודי KEGG מופרדים בפסיק (ברירת מחדל `cre,cvr,mng`) |
| `target_organism` | האורגניזם שחוקרים בפועל, אם הקודים הם תחליף (ברירת מחדל `Coelastrella sp.`); ריק אם הם האורגניזם עצמו |
| `date` | משאירים ריק |

**מקבלים:** Release לכל אורגניזם, `<code>_kegg_<date>`, ובו מודל ל‑**mummichog** (`.json`), GMT של קבוצות תרכובות להעשרה **מבוססת‑ID** (`.compound_pathway.gmt`), ו‑manifest עם גרסאות ו‑checksums.

**משתמשים ב:** פייפליין ה‑**מטבולומיקה** — `modes.metabolomics.enrichment.gmt_file` (הקובץ שהורד, יחסי ל‑`paths.raw`) ו‑`modes.metabolomics.enrichment.mummichog.model_ref` (URL + sha256 מה‑Release; הפייפליין מוריד ומאמת לבד). דוגמה מלאה בסעיף "הרצת הפייפליין" למטה.

---

## 4. Delete releases

Release לא ניתן לדריסה, אז כדי לבנות מחדש באותו שם מוחקים קודם. ממלאים את שם ה‑Release (או כמה, מופרדים בפסיק) ב‑`tags`, וכותבים `DELETE` ב‑`confirm`. ה‑Release וה‑tag שלו נמחקים.

---

## טוב לדעת

- **Release הוא קבוע.** הרצה חוזרת של אותו Action באותו יום ובאותו שם נכשלת בכוונה. מוחקים קודם את הישן, או משתמשים ב‑`id_label` / `date` אחר.
- **רק דבר אחד ידני:** קבצים שדורשים התחברות כדי להוריד (למשל GiardiaDB) או שמגיעים מטופס באתר (KAAS). כל מה שציבורי — KEGG, Gene Ontology, NCBI, Ensembl, UniProt — ה‑Action מוריד לבד.
- **הכל מתועד.** כל Release מפרט את הקלטים עם checksums, את גרסאות המאגרים, ואת ה‑commit והריצה שבנו אותו — כך שאפשר לשחזר בדיוק ממה נבנתה כל תוצאה.

---

# פרטים נוספים — מטבולומיקה, sha256, ותחזוקה

## התמונה הגדולה

שני ריפוזיטוריז עובדים יחד. אחד **בונה** קבצי ייחוס לאורגניזם, והשני (הפייפליין) **צורך** אותם. הזרימה תמיד אותו דבר, והקבצים עצמם **לא** נשמרים בתוך ה‑git — הם artifacts שמתפרסמים כ‑GitHub Releases.

```
multiomics-annotation-prep  ──►  GitHub Releases  ──►  multiomic-core
        (הבנאי)                     (המדף)                (הפייפליין)
```

| שלב | מה קורה |
|-----|---------|
| **הבנאי** (`multiomics-annotation-prep`) | מוריד מ‑KEGG ובונה, לכל אורגניזם: מודל מטבולי (`.json`) עבור mummichog, GMT של קבוצות תרכובות להעשרה מבוססת‑ID, טבלה קריאה, ו‑`manifest` עם sha256. |
| **המדף** (GitHub Releases) | כל בנייה מתפרסמת כ‑Release מתוארך ובלתי‑משתנה. תג: `<code>_kegg_<date>` — 4 קבצים לכל אורגניזם. |
| **הפייפליין** (`multiomic-core`) | צורך את המודל דרך `model_ref` (URL+sha256 — נשלף ומאומת אוטומטית), ואת ה‑GMT דרך `gmt_file` (קובץ מקומי). |

**למה שני פורמטים לאותו אורגניזם?** המטבולומיקה רצה בשתי דרכים: מבוססת‑**m/z** (mummichog, עם המודל `.json`) ומבוססת‑**ID** (ORA/GSEA/QEA, עם ה‑GMT). שניהם נבנים מאותו מקור KEGG — אז אותם מסלולים ואותה ביולוגיה.

---

## מה זה sha256?

`sha256` הוא **"טביעת אצבע" דיגיטלית של קובץ** — מחרוזת של 64 תווי hex שמחושבת מהתוכן:

- **ייחודי לתוכן** — שינוי של בית אחד → sha256 שונה לגמרי.
- **חד‑כיווני** — אי אפשר לשחזר את הקובץ מה‑sha.
- **דטרמיניסטי** — אותו קובץ נותן תמיד אותו sha, בכל מחשב.

**למה זה נחוץ:** ב‑`model_ref` יש `url` (מאיפה להוריד) ו‑`sha256` (בדיוק איזה קובץ מצפים לקבל). כשהפייפליין מוריד את ה‑`.json`, הוא מחשב מחדש את ה‑sha256 ומשווה: תואם → זה הקובץ הנכון; לא תואם → נעצר בקול רם במקום להריץ מודל שגוי. זה מה שהופך את העבודה ל‑reproducible.

> GitHub מציג את זה כ‑`sha256:abc…`. התחילית `sha256:` היא רק שם האלגוריתם — בקונפיג מכניסים **רק את 64 ה‑hex** (ראו מלכודת 1).

---

## מה מריצים

### א. בניית קבצים לאורגניזם — דרך ה‑Workflow (מומלץ)

אין צורך במחשב מקומי או בהתקנות. הבנייה רצה על שרתי GitHub שיש להם גישת רשת ל‑KEGG והרשאה לפרסם Releases.

1. בריפו **multiomics-annotation-prep** → לשונית **Actions** → **Build & publish organism artifacts** → **Run workflow**.
2. ממלאים:
   - `organisms` — קודי KEGG מופרדים בפסיק (למשל `cre,cvr,mng`)
   - `target_organism` — ברירת מחדל `Coelastrella sp.`
   - `date` — ריק = היום (UTC)
3. התוצאה: Release לכל אורגניזם עם 4 הקבצים, אחרי אימות אוטומטי (mummichog smoke + בדיקת מסות).

### ב. בנייה מקומית (חלופה — לפיתוח/דיבוג)

דורש מחשב עם גישת רשת ל‑`rest.kegg.jp`.

```bash
# התקנה חד-פעמית
git clone https://github.com/multiomics-center-Israel/multiomics-annotation-prep
cd multiomics-annotation-prep
python -m venv .venv && source .venv/bin/activate
pip install requests pyyaml
pip install -r requirements-mummichog.txt

# בנייה: מודל + GMT יחד, לאורגניזם אחד
python scripts/run_mummichog_model.py --kegg-code cre \
  --model-organism "Chlamydomonas reinhardtii" \
  --target-organism "Coelastrella sp." \
  --emit-compound-sets --validate --out results --cache data
```

### ג. הרצת הפייפליין (multiomic-core)

בקונפיג המטבולומי, תחת `enrichment`, מצביעים על המודל וה‑GMT ומריצים כרגיל.

```yaml
modes:
  metabolomics:
    enrichment:
      run_enrichment: true
      gmt_file: "metabolomics/cre_kegg_20260723.compound_pathway.gmt"  # קובץ מקומי
      mummichog:
        enabled: true
        ionization_mode: positive        # pos_default | positive | negative
        model_ref:
          url:    "https://github.com/.../releases/download/cre_kegg_20260723/cre_kegg_20260723.json"
          sha256: "344851a4…ed68"        # 64 hex נטו, בלי sha256:
```

---

## איך מוסיפים אורגניזם חדש

התהליך המלא, מקוד KEGG ועד ריצה בפייפליין — חמישה שלבים:

1. **מצא את קוד ה‑KEGG.** ב‑[genome.jp/kegg](https://www.genome.jp/kegg/) — קוד בן 3–4 אותיות (למשל `cre` = Chlamydomonas reinhardtii). ודא שלאורגניזם יש גנום ב‑KEGG. אם אין — בוחרים אורגניזם קרוב כ‑surrogate.
2. **הרץ את ה‑Workflow.** Actions → **Run workflow**: `organisms=<code>`, `target_organism` לפי הצורך, `date` ריק. לשם מלא בתיעוד, אפשר להוסיף את הקוד למַפָּה שבתוך `publish-organism-artifacts.yml` (`NAMES`); אחרת ייבנה עם שם ריק — עדיין תקין.
3. **המתן לסיום (~5 דק').** נוצר Release `<code>_kegg_<date>` עם 4 קבצים. ריצה כושלת נעצרת עם הודעה ברורה; אותו תאריך שכבר קיים ייכשל בכוונה (artifacts הם immutable).
4. **אסוף sha256 והורד את ה‑GMT.** ב‑Release, ה‑`digest` של קובץ ה‑`.json` הוא ה‑sha256 של המודל. הורד את קובץ ה‑`.compound_pathway.gmt` למקום מקומי (הפייפליין לא מוריד GMT מ‑URL).
5. **חבר בקונפיג של הפייפליין והרץ.** `model_ref` = ה‑URL של ה‑`.json` + ה‑sha256 (נקי), ו‑`gmt_file` = הנתיב המקומי של ה‑GMT. מודל ו‑GMT — מאותו קוד אורגניזם. הרץ.

---

## מלכודות נפוצות

חמש טעויות שקל ליפול בהן — שווה לעבור עליהן לפני ריצה:

1. **ה‑sha256 — 64 hex נטו, בלי `sha256:`.** ה‑`digest` ב‑GitHub מוצג כ‑`sha256:abc…`. מעתיקים רק את ה‑hex. אחרת: `must be a 64-character hex sha256 digest`.
2. **`gmt_file` הוא נתיב מקומי — לא URL.** המודל נשלף אוטומטית מ‑URL; ה‑GMT לא. מורידים את הקובץ ומצביעים עליו בנתיב מקומי (מוחלט או יחסי לפרויקט).
3. **מודל ו‑GMT מאותו אורגניזם.** באותה ריצה, `model_ref` ו‑`gmt_file` צריכים להיות מאותו קוד — אחרת ההעשרה מבוססת‑m/z ומבוססת‑ID מדברות על ביולוגיות שונות.
4. **`ionization_mode`: `pos_default`, `positive` או `negative` בלבד.** ערך לא מוכר (כמו `mixed`) נבלע בשקט ל‑`pos_default` (ברירת המחדל של mummichog) — בלי שגיאה, אבל בפולריות שאולי לא התכוונת אליה. mummichog מריץ פולריות אחת בכל ריצה.
5. **עדכון = תאריך חדש, לא דריסה.** ה‑artifacts הם immutable. כדי לרענן (עדכון KEGG / תיקון) מריצים את ה‑workflow עם `date` חדש — נוצר Release חדש, הישן נשאר.

---

## הייחוסים הנוכחיים

שלוש אצות ירוקות שנבנו כ‑surrogates ל‑**Coelastrella sp.** (תג `<code>_kegg_20260723`). הספירות: compounds / reactions / pathways במודל, ו‑pathways / compounds ב‑GMT.

| אורגניזם | קוד | מודל (cpd/rxn/path) | GMT (path/cpd) | sha256 של המודל (`.json`) |
|----------|-----|---------------------|----------------|---------------------------|
| Chlamydomonas reinhardtii | `cre` | 1111 / 1322 / 80 | 81 / 1170 | `344851a45d310d4e6712bcea10a0a3fa5d7c31ed82ad19a8a98d1ac2e3ebed68` |
| Chlorella variabilis | `cvr` | 1150 / 1379 / 82 | 84 / 1229 | `7d298afa195b53a0ecce853c15401278b43ad4a4392e2dbb7d3fba13a882480f` |
| Monoraphidium neglectum | `mng` | 1057 / 1152 / 81 | 83 / 1139 | `182d524e38808ff44e327c5f4ebb2a15378471e09337fb29f0a99f21a2e96c8c` |

דפוס ה‑URL: `…/releases/download/<tag>/<tag>.json` (מודל) ו‑`…/<tag>.compound_pathway.gmt` (GMT). ה‑sha256 של ה‑GMT זמין ב‑`digest` של ה‑Release (לאימות ידני אחרי הורדה).

---

## כמה מקום ה‑Releases תופסים?

מעט מאוד. כל אורגניזם תופס **~0.7 MB**, ושלושת הייחוסים הנוכחיים יחד ≈ **2.1 MB**.

| קובץ | גודל טיפוסי |
|------|-------------|
| `.json` (מודל) | ~0.5 MB |
| `.pathway2compound.tab` | ~0.16 MB |
| `.compound_pathway.gmt` | ~21 KB |
| `.manifest.json` | ~3.4 KB |
| **סה"כ לאורגניזם** | **~0.7 MB** |

- הרוב הוא ה‑`.json` (הוא מכיל את כל ה‑compounds/reactions/pathways); ה‑GMT זעיר.
- כל **גרסה מתוארכת נוספת** מוסיפה עוד ~0.7 MB לאורגניזם (Release חדש, לא דריסה).
- **לא נספר בגודל ה‑git** — Release assets נשמרים באחסון נפרד של GitHub עם מכסה נדיבה. אפשר לפרסם מאות גרסאות בלי להתקרב לגבול.

---

## RNA / GO enrichment — טבלאות מבוססות‑גנים (נתיב נפרד)

ההעשרה של **RNA-seq** (GO + KEGG, ORA/GSEA) **לא** משתמשת ב‑Release עם URL+sha כמו המודל המטבולומי — היא קוראת **תיקייה מקומית** של טבלאות. הסיבה: הטבלאות ממופתחות לפי **מזהי הגנים של הטרנסקריפטום**, לא לפי compound IDs אוניברסליים.

**מה ה‑config מצפה** (תחת `modes.rna.enrichment`):

```yaml
enrichment:
  enabled: true
  annotation_dir: "data/Func_annot_data_expanded"   # תיקייה מקומית עם הקבצים למטה
  databases: ["KEGG", "GO_BP", "GO_MF", "GO_CC"]     # מה שחסר — מדולג עם warning
```

עד 8 קבצי `.tab` (שתי עמודות; ה‑loader לוקח את שתי העמודות הראשונות, שורת ה‑header לא משנה):

| database | TERM2GENE | TERM2NAME |
|----------|-----------|-----------|
| KEGG | `KEGG_pathway2gene.tab` | `KEGG_pathway2name.tab` |
| GO_BP | `GO2gene_BP.tab` | `GO2name_BP.tab` |
| GO_MF | `GO2gene_MF.tab` | `GO2name_MF.tab` |
| GO_CC | `GO2gene_CC.tab` | `GO2name_CC.tab` |

**`annotation-prep` כבר מייצר בדיוק את זה** — אין צורך בקוד חדש:

```bash
# שתי טבלאות ה-KEGG (מ-KAAS)
python scripts/run_kegg_nonmodel.py --kaas <query.ko.txt> --out annot_dir --cache data
# שש טבלאות ה-GO (מ-Trinotate), עם הרחבת היררכיה — "pre-expanded"
python scripts/run_go.py --go-table <trinotate_go.txt> --out annot_dir --cache data
```

מריצים את שניהם עם אותו `--out`, מקבלים את 8 הקבצים בתיקייה אחת, ומכוונים אליה את `annotation_dir`.

**אוטומטי (workflow):** אפשר גם דרך GitHub Actions — **Build gene annotation (RNA / proteomics)**. מעלים את `query.ko.txt` ו/או `trinotate_go.txt` (או קובץ `.gaf`, ראו למטה) ל‑`annotation_inputs/<project>/` (Add file → Upload files), מריצים את ה‑workflow עם שם הפרויקט ו‑`omics: rna`, ומקבלים Release `rna_annot_<project>_<date>` עם קבצי ה‑annotation (+ zip). מורידים, פורסים לתיקייה, ומכוונים אליה את `annotation_dir`. (הדגל `strip_isoform` שם — כמו הגוצ'ה למטה.)

> ⚠️ **מזהי הגנים חייבים להתאים ל‑counts.** ה‑loader (`R/core/09_enrichment.R`) בודק חפיפה בין הגנים בטבלאות למזהי הגנים במטריצת ה‑counts, ומזהיר אם החפיפה <5% (*"Check that gene ID types match"*). לכן הטבלאות נבנות מ‑KAAS + Trinotate של **אותו** פרויקט — כדי שמרחב ה‑IDs יהיה זהה.

**למה זה נתיב נפרד מהמטבולומיקה:**

| | מטבולומיקה (model / GMT) | RNA (GO / KEGG) |
|---|---|---|
| מפתח | compound IDs (KEGG) | gene IDs (טרנסקריפטום) |
| מקור | קוד אורגניזם KEGG | KAAS + Trinotate של הפרויקט |
| הפצה | Release מתוארך (URL+sha) | תיקייה מקומית לפרויקט |

*(מבוסס על PR #128 ב‑`multiomic-core`, `feature/enrichment-migration-v2`.)*

---

## נספח — הכנת הקלט לאנוטציה (KAAS / Trinotate / Ensembl / UniProt)

השלב היחיד שהוא ידני, *לפני* הריפו: להשיג את קבצי הקלט. הבחירה תלויה אם יש גנום מוער או רק טרנסקריפטום. **החוט המקשר בכל המסלולים:** מזהי הגנים בקבצי האנוטציה חייבים להיות אותו מרחב IDs כמו עמודת `gene_id` ב‑counts.

### א. לא‑מודל — טרנסקריפטום de novo

**KAAS → `query.ko.txt` (KEGG).** [KEGG KAAS](https://www.genome.jp/kegg/kaas/) מקצה KO לרצפים.
- קלט: FASTA, מומלץ **חלבונים** (ORFs מ‑TransDecoder).
- בחירות: `GHOSTX` (מהיר); שיטת **SBH** לטרנסקריפטום (BBH לגנום שלם); רשימת אורגניזמים ייחוס קרובים (אצות/צמחים).
- פלט: `gene ⇥ KO` — משנים שם ל‑`query.ko.txt`. לסטים גדולים: BlastKOALA/GhostKOALA.

**Trinotate → `trinotate_go.txt` (GO).** [Trinotate](https://github.com/Trinotate/Trinotate) מאחד blastx/blastp מול Swiss-Prot, hmmscan מול Pfam ועוד.
- פייפליין: Trinity → TransDecoder → חיפושי הומולוגיה → טעינה ל‑SQLite → `Trinotate --report`.
- חילוץ GO: הסקריפט `extract_GO_assignments_from_Trinotate_xls.pl` נותן `gene ⇥ GO`. ה‑parser של `run_go.py` עמיד לפורמט Trinotate המלא (מחלץ `GO:\d{7}`).

**סדר הפעולות:** Trinity → TransDecoder → (KAAS + Trinotate במקביל) → שני הקבצים → העלאה ל‑`annotation_inputs/<project>/` → הרצת ה‑workflow.

**GAF מאתר אורגניזם (למשל GiardiaDB / VEuPathDB) → קובץ `.gaf` (GO).** אתרי VEuPathDB מפרסמים קובץ GO לכל אורגניזם, אבל ההורדה דורשת התחברות לאתר, ולכן מורידים ידנית ומעלים:
- מתחברים לאתר, מורידים את קובץ ה‑GO GAF של האורגניזם, ו**לא משנים את שמו** (למשל `GiardiaDB-71_GintestinalisAssemblageAWB_GO.gaf`) — מספר הגרסה בשם מתועד ב‑Release. קובץ GAF אחד לתיקייה.
- אם המזהים בנתונים שונים ממזהי האתר (למשל חלבונים `KAE…` מ‑NCBI מול locus tags `GL50803_…`), ממלאים ב‑workflow את `ncbi_assembly` (למשל `GCA_000002435.2`): ה‑runner מוריד בעצמו את ה‑feature table הציבורי של ה‑assembly (ממפה `locus_tag` → `product_accession`) ומתעד את הכתובת וה‑sha256 — אין צורך להעלות אותו ידנית. למיפוי אחר — מעלים `id_map.txt` (או `.gz`): טבלה בשתי עמודות `מזהה_ב‑GAF ⇥ מזהה_בנתונים`.
- **תיעוד:** ה‑Release כולל את ה‑commit, קישור לריצה, sha256 של כל קובץ קלט ואת שם וכותרת קובץ ה‑GAF (מקור + גרסה), ו‑`GO_provenance.json` עם הגרסאות וה‑checksums של ה‑GAF, המיפוי וה‑GO ontology.
- קובץ ה‑GAF ו‑`trinotate_go.txt` הם חלופות — לא שניהם באותה תיקייה.

**KEGG בלי KAAS:** לאורגניזם ש‑KEGG כבר מאנטט (למשל Giardia — `gla`) ממלאים ב‑workflow את `kegg_org`. ה‑runner מושך את הקישור גן→pathway מ‑KEGG REST, ממפה אותו עם אותו `ncbi_assembly`/`id_map` של ה‑GO, וכותב `KEGG_pathway.gmt` + `KEGG_provenance.json` (כולל גרסת KEGG). לא יחד עם `query.ko.txt`.

**לחלבונים:** מריצים את אותו workflow עם `omics: proteomics`. מקבלים Release `prot_annot_<project>[_<id_label>]_<date>` (השדה `id_label`, למשל `kae`, מציין את סוג המזהים; ריק עם `ncbi_assembly` → `ncbi`), ומכוונים את `modes.proteomics.pathway.gmt_file` לקבצי `GO_*.gmt` (ו‑`KEGG_pathway.gmt`). אם ה‑GMT ממופתח ל‑Protein.Group כמו שהוא — `skip_annotation: true`.

### ב. מודל — גנום ב‑Ensembl או פרוטאום ב‑UniProt

לא צריך KAAS/Trinotate — מורידים GO ישירות. שני מודולים, אותם קבצי פלט (`GO2gene_{BP,MF,CC}.tab` + `GO2name_…`, מורחבים היררכית → נכנסים לאותו `annotation_dir`):

- **Ensembl / BioMart** (`run_ensembl.py`) — לפי Ensembl gene ID (שאילתה ישירה ל‑BioMart REST, בלי חבילות נוספות). ⚠️ בוחרים את חלוקת BioMart הנכונה (`mart: plants_mart` + host מתאים לצמח/אצה, לא ה‑mart הדיפולטי).
- **UniProt** (`run_uniprot.py`) — לפי `taxon_id` (Swiss-Prot/reviewed), ממופתח ל‑UniProt accessions.

**workflow ייעודי:** GitHub Actions → **Build model-organism annotation** → בוחרים `source: ensembl|uniprot`, ממלאים `dataset` (Ensembl) או `taxon_id` (UniProt), ואופציונלית `kegg_org` (גם טבלאות KEGG). מקבלים Release `<label>_<source>_annot_<date>` עם קבצי ה‑annotation. *(או מקומית: `run_ensembl.py` / `run_uniprot.py`, או בלוק config + `run_all.py`.)*

### התאמת רמת ה‑ID (חוזר בכל מסלול)

מזהי Trinity: `TRINITY_DN1000_c0_g1` (גֵן) מול `..._g1_i1` (טרנסקריפט). הדגל `strip_isoform` (ברירת מחדל: מסיר `_iN`) מיישר לרמת גֵן. תמיד מרימים את קבצי האנוטציה לאותה רמה כמו ה‑counts — אחרת אזהרת ה‑<5% חפיפה.

---

**קישורים:** [multiomics-annotation-prep](https://github.com/multiomics-center-Israel/multiomics-annotation-prep) (בנאי) · [multiomic-core](https://github.com/multiomics-center-israel/multiomics-core) (פייפליין) · חוזה הפורמט: `MODEL_CONTRACT.md`
