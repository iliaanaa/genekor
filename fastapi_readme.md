# Πώς λειτουργεί το FastAPI (σε βασικό επίπεδο):
Ορίζεις endpoints με decorators όπως @app.get() ή @app.post(), αντίστοιχα για GET και POST αιτήματα.

Ορίζεις τα δεδομένα εισόδου/εξόδου με Python τύπους ή Pydantic models (π.χ. για body, query parameters, path variables).

Το FastAPI δημιουργεί αυτόματα τεκμηρίωση (Swagger UI και ReDoc) για το API σου.



# Ροή εκτέλεσης:
Έρχεται ένα HTTP αίτημα (π.χ. από browser ή frontend app).

Το FastAPI το δρομολογεί στο κατάλληλο route με βάση το URL και το HTTP method.

Τα δεδομένα εισόδου αναλύονται και επικυρώνονται με Pydantic.

Η συνάρτηση (route handler) εκτελείται, και επιστρέφει αποτέλεσμα (συνήθως σε JSON).

Η απάντηση αποστέλλεται πίσω στον client.


# Παράδειγμα:
python
from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI()

class User(BaseModel):
    username: str
    email: str

@app.get("/")
def read_root():
    return {"message": "Καλώς ήρθες στο FastAPI!"}

@app.post("/users/")
def create_user(user: User):
    return {"user_created": user}


Εάν επισκεφτείς http://localhost:8000/docs αφού τρέξεις την εφαρμογή, θα δεις διαδραστική τεκμηρίωση (Swagger UI) που δημιουργείται αυτόματα!

# Εκκίνηση εφαρμογής:
bash
uvicorn main:app --reload
(όπου main.py είναι το όνομα του αρχείου)


# Πλεονεκτήματα FastAPI
Πολύ γρήγορο (χτισμένο πάνω σε uvicorn, Starlette).

Υποστηρίζει πλήρως async/await.

Αυτόματη τεκμηρίωση (Swagger, ReDoc).

Εύκολη επικύρωση εισόδου με Pydantic.

Ιδανικό για REST APIs, microservices, και ML εφαρμογές.


# Τι είναι τα endpoints και τα @app.get(), @app.post() στο FastAPI
🔹 Τι είναι endpoint;
Ένα endpoint είναι μια "διεύθυνση" στο API σου — δηλαδή ένα URL που συνδέεται με μια συγκεκριμένη λειτουργία του backend σου.



# Παράδειγμα:
Αν έχεις αυτό το endpoint:

python
@app.get("/users")
def get_users():
    return {"message": "Όλοι οι χρήστες"}
Τότε όταν κάποιος στείλει αίτημα GET στο:

bash
http://localhost:8000/users
θα πάρει πίσω αυτό το μήνυμα.

# Τι είναι το @app.get() και το @app.post();
Αυτά είναι decorators (ειδικές εντολές που "τυλίγουν" συναρτήσεις) και δηλώνουν:

@app.get("/κάποιο_path") → Αντιστοιχεί σε HTTP GET αίτημα. Δηλαδή, για λήψη δεδομένων.

@app.post("/κάποιο_path") → Αντιστοιχεί σε HTTP POST αίτημα. Δηλαδή, για αποστολή δεδομένων.




Παράδειγμα:
python
from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "Αρχική σελίδα"}

@app.post("/submit")
def submit_data():
    return {"message": "Δεδομένα καταχωρήθηκαν"}
Αν μπεις στο /, θα πάρεις το μήνυμα "Αρχική σελίδα" (με GET).

Αν στείλεις POST στο /submit, θα πάρεις "Δεδομένα καταχωρήθηκαν".

# GET vs POST
Μέθοδος	Πότε τη χρησιμοποιούμε
GET	Για να πάρουμε/διαβάσουμε δεδομένα
POST	Για να στείλουμε/δημιουργήσουμε δεδομένα\


# τι είναι η "διαδραστική τεκμηρίωση (Swagger UI)
Τι είναι η Swagger UI στο FastAPI;
Η Swagger UI είναι ένα γραφικό περιβάλλον στο browser, που δημιουργείται αυτόματα από το FastAPI και σου επιτρέπει:

✅ Να βλέπεις όλα τα endpoints του API σου

✅ Να δοκιμάζεις εύκολα GET/POST/PUT/DELETE αιτήματα χωρίς να γράφεις κώδικα

✅ Να βλέπεις ποιες παραμέτρους δέχεται κάθε endpoint

✅ Να βλέπεις παραδείγματα απαντήσεων

✅ Να κάνεις debug ή έλεγχο του API σου

🔗 Πού τη βρίσκω;
Όταν τρέχεις την FastAPI εφαρμογή σου με uvicorn, άνοιξε τον browser σου και πήγαινε στο:

bash
Copy
Edit
http://localhost:8000/docs
Εκεί θα εμφανιστεί η Swagger UI ✨

📸 Παράδειγμα εμφάνισης:
Φαντάσου ένα πίνακα με όλα τα endpoints, π.χ.:

bash
Copy
Edit
GET     /users          ➜ Επιστρέφει όλους τους χρήστες
POST    /users          ➜ Προσθέτει έναν νέο χρήστη
GET     /users/{id}     ➜ Παίρνει έναν χρήστη με id
Κάθε ένα έχει κουμπί "Try it out" για να το δοκιμάσεις live, στέλνοντας δεδομένα και βλέποντας την απάντηση.

✅ Πλεονεκτήματα:
Δεν χρειάζεσαι Postman ή άλλο εργαλείο — έχεις δικό σου interface για δοκιμές.

Ιδανικό για ομάδες developers που θέλουν να τεκμηριώσουν γρήγορα ένα API.

Αν συνδυαστεί με Pydantic, εμφανίζονται αυτόματα τύποι δεδομένων, validation errors κλπ.

# Αν θες, δοκίμασέ το:
Φτιάξε ένα αρχείο main.py με αυτό:

python
from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def root():
    return {"message": "Γεια σου από FastAPI!"}


Τρέξε:

bash
uvicorn main:app --reload
Πήγαινε στον browser:

bash
http://localhost:8000/docs

# Υπάρχει και εναλλακτική τεκμηρίωση στο:

bash
http://localhost:8000/redoc
(Είναι πιο "αναγνωστική", όχι διαδραστική).