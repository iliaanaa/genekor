# Σημαντικές Σημειώσεις:
#Δομή Βάσης Δεδομένων:

Ο πίνακας brca1_variants χρησιμοποιεί το variation_id ως primary key

Τα κριτήρια ACMG αποθηκεύονται ως JSONB

Οι διενέξεις (conflicting interpretations) αποθηκεύονται ως JSONB

Τα RCV accession IDs αποθηκεύονται ως πίνακας κειμένων (TEXT[])

# Φιλτράρισμα:

Κρατάμε μόνο μεταλλάξεις του BRCA1

Φιλτράρουμε βάσει assembly (π.χ. GRCh38)

Διατηρούμε όλες τις υποβολές για conflicting interpretations

# Αποθήκευση HGVS:

HGVS_c για DNA changes (π.χ. c.1510A>G)

HGVS_p για protein changes (π.χ. p.Arg504Gly)

# Επεκτάσεις:

Προσθέστε περισσότερα ACMG criteria στη συνάρτηση apply_acmg_criteria

Προσαρμόστε το significance_map για επιπλέον κλινικές ταξινομήσεις

# Ανάλυση Διπλότυπων:

Το VariationID εγγυάται μοναδικότητα

Ομαδοποιούμε υποβολές με ίδιο VariationID για conflicting interpretations

