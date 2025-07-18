# Οδηγίες: Πώς να κάνετε Push στα Branch σας

## Αρχικό Checkout του Δικού σας Branch

Αν δεν έχετε ακόμα το τοπικό branch:
```
git fetch origin
git checkout branch-onoma
```
Παράδειγμα
```
git checkout branch2-xim
```

## Προσθήκη και Αποστολή Αλλαγών (Push)

Μόλις κάνετε αλλαγές:
```
git add .
git commit -m "Σχόλιο για τις αλλαγές"
git push origin branch-onoma
```

Παράδειγμα:
```
git add .
git commit -m "Διόρθωση στο dataparse"
git push origin branch2-xim
```

##

* Μην κάνετε push στο `main` 
* Αν δείξει conflict, ενημερώνουμε και κάνουμε rebase ή merge.

## Branches Ομάδας

* branch1-iliana → Iliana
* branch2-xim → Xim
* branch3-aggeliki → Υaggeliki


