#RawDogging variations :P :P :P
#PS1
'''
KAINOURGIO VARIANT!! -- DHLADH c.HGVS poy prokalei gnwsto p.HGVS
(8a doyme parakatw ti ennoyme gnwsto p.hgvs)
TI EINAI??

TO BAROS APODEIKSHS GIA MISSENSE PEFTEI SE AYTON POY TO BRHKE KAI 8ELEI NA TO TESTAREI GIA PS1
OPWS KAI TO BAROS GIA TO GEGONOS OTI H NOUKLEOTIDIKH ALLAGH DEN FEREI PI8ANOTHTA GIA:
Cryptic splice site,    Splice enhancer/silencer disruption,   mRNA structure/stability impact

EXONTAS AYTO YPOPSHN ELEGXOYME AN YPARXEI ALLO
1 (ENA) (einai arketo) 
KALA SUPPORTED PATHOGENIC VARIANT (DHLADH 2 ASTERIA KAI ANW), (ASTERIA >= 2)

EPISHS AFOY DE 8A BREI TO C.HGVS 8A ELEKSOYME ME TO P.HGVS
'''
def stars_from_review_status(review_status: str) -> int:
    """Metatrepei to review_status se arithmo asteriwn ClinVar"""
    review_status = review_status.strip().lower()
    if review_status in {"no assertion provided", "no assertion criteria provided"}:
        return 0
    elif review_status in {"criteria provided, single submitter", "criteria provided, conflicting interpretations"}:
        return 1
    elif review_status in {"criteria provided, multiple submitters, no conflicts", "reviewed by expert panel"}:
        return 2
    elif review_status == "reviewed by practice guideline":
        return 3
    elif review_status == "reviewed by professional society":
        return 4
    return 0
