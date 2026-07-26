# מדריך צוות — בניית ייחוס לאורגניזם והרצת הפייפליין

> מדריך תפעולי בעברית: מה הריפו עושה, מה מריצים, ואיך מוסיפים אורגניזם חדש מקצה לקצה.
>
> **בנאי:** `multiomics-annotation-prep` · **פייפליין:** `multiomic-core` · **מקור נתונים:** KEGG REST

---

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
        ionization_mode: positive        # positive | negative
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
4. **`ionization_mode`: `positive` או `negative` בלבד.** ערך לא מוכר (כמו `mixed`) נבלע בשקט ל‑`positive` — בלי שגיאה, אבל בפולריות שאולי לא התכוונת אליה. mummichog מריץ פולריות אחת בכל ריצה.
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

**קישורים:** [multiomics-annotation-prep](https://github.com/multiomics-center-Israel/multiomics-annotation-prep) (בנאי) · [multiomic-core](https://github.com/multiomics-center-israel/multiomics-core) (פייפליין) · חוזה הפורמט: `MODEL_CONTRACT.md`
