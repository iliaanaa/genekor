# Οδηγίες: Πώς να κάνετε Push στα Branch σας

## 1)Αρχικό Checkout του Δικού σας Branch

Αν δεν έχετε ακόμα το τοπικό branch:
```
git fetch origin
git checkout branch-onoma
```
Παράδειγμα
```
git checkout branch2-xim
```

Ανάλυση: Πώς δουλεύει το git checkout branch-onoma και το git fetch origin && git checkout branch-onoma
Περίπτωση 1: Το branch υπάρχει ήδη τοπικά
Εντολή:
```
git checkout branch-onoma
```
Τι κάνει:
Μεταφέρεται τοπικά στο υπάρχον branch, αν αυτό υπάρχει ήδη στο σύστημά σας.



Περίπτωση 2: Το branch υπάρχει στο remote (GitHub), αλλά όχι τοπικά
Εντολή (2 βήματα ή σε ένα):
```
git fetch origin
git checkout branch-onoma
```
Πότε χρησιμοποιείται:
Όταν κάποιος έχει ήδη κάνει fetch/pull στο παρελθόν ή έχει δημιουργήσει το branch ο ίδιος.

Τι κάνει το git fetch origin:
Φέρνει όλες τις τελευταίες αλλαγές και branches από το απομακρυσμένο repository (GitHub), χωρίς να αλλάξει κάτι στο working directory σας ακόμα.

Τι κάνει το git checkout branch-onoma μετά το fetch:
Αν το branch τώρα υπάρχει μετά το fetch, μεταφέρεστε σε αυτό.

Γιατί χρειάζεται fetch;
Εάν κάποιος άλλος δημιούργησε ένα νέο branch και εσείς δεν το βλέπετε ακόμα τοπικά. Το fetch ενημερώνει τα refs/heads του τοπικού repository


Σημείωση για πρώτη φορά κλωνοποίηση:
Αν κάποιος κάνει:
```
git clone https://github.com/iliaanaa/genekor.git
```
Τότε θα πάρει μόνο το default branch (main) εκτός αν:

Είτε κάνει fetch μετά.

Είτε χρησιμοποιεί:
```
git clone --branch branch-onoma https://github.com/iliaanaa/genekor.git
```

!!!
Αν κάποιος δεν έχει ξαναδουλέψει με το συγκεκριμένο repository (repo) και δεν το έχει καθόλου τοπικά στον υπολογιστή του, πρέπει πρώτα να κάνει git clone.
Τι σημαίνει αυτό πρακτικά:
Το git clone κατεβάζει όλα τα αρχεία, το ιστορικό των commits, και τα branches από το GitHub στον υπολογιστή του χρήστη.
Έτσι αποκτά έναν τοπικό φάκελο που είναι συνδεδεμένος με το απομακρυσμένο repository (remote origin).

## 2)Προσθήκη και Αποστολή Αλλαγών (Push)

Μόλις κάνετε αλλαγές:
```
git add .
git commit -m "Σχόλιο για τις αλλαγές"
git push origin branch-onoma
```
Όταν κάνεις git push origin με HTTP, αντί για τον κωδικό σου GitHub, χρησιμοποιείς το PAT (Personal Access Token),

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
* branch3-aggeliki →aggeliki


