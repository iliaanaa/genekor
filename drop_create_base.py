import psycopg2
from psycopg2 import sql

def drop_and_create_database(dbname="clinvar_db"):
    '''Ορίζουμε τη συνάρτηση drop_and_create_database.

Δέχεται όρισμα dbname με προεπιλογή "clinvar_db", δηλαδή μπορείς να την καλέσεις και για άλλη βάση αν θέλεις.'''

    print(f"Dropping και re-creating τη βάση '{dbname}'...")


'''Συνδέεται στη βάση postgres (όχι στη clinvar_db!) γιατί δε μπορείς να κάνεις DROP μια βάση στην οποία είσαι συνδεδεμένος.

Τα στοιχεία σύνδεσης έρχονται από ένα dict DB_CONFIG που περιέχει user, password, host, port.'''
    try:
        # Σύνδεση στην default βάση 'postgres'
        conn = psycopg2.connect(
        dbname="postgres",
        user=DB_CONFIG["user"],
        password=DB_CONFIG["password"],
        host=DB_CONFIG["host"],
        port=DB_CONFIG["port"]
        )
        '''Ενεργοποιούμε autocommit mode, ώστε οι εντολές DROP και CREATE DATABASE να εκτελούνται αμέσως χωρίς conn.commit().'''
    conn.autocommit = True


'''Δημιουργούμε cursor (αντικείμενο εκτέλεσης εντολών SQL).

Η with φροντίζει να κλείσει σωστά ο cursor μόλις τελειώσουμε.'''
with conn.cursor() as cur:
            # Τερματισμός ενεργών συνδέσεων
        '''Εντολή PostgreSQL που τερματίζει όλες τις ενεργές συνδέσεις στην βάση clinvar_db, εκτός από την τρέχουσα.

pg_stat_activity: προβολή με όλες τις συνδέσεις στη βάση.

pg_terminate_backend(pid): σκοτώνει τη σύνδεση.

%s: placeholder (προστατεύει από SQL injection).'''
cur.execute(sql.SQL("""
                SELECT pg_terminate_backend(pid)
                FROM pg_stat_activity
                WHERE datname = %s AND pid <> pg_backend_pid();
            """), [dbname])

            # Drop και Create της βάσης
            '''Εντολή που κάνει DROP της βάσης clinvar_db αν υπάρχει.

Το sql.Identifier(dbname) βάζει το όνομα της βάσης με ασφάλεια 
(αντί να το "κολλήσουμε" απευθείας στο string, που είναι επικίνδυνο).'''
            cur.execute(sql.SQL("DROP DATABASE IF EXISTS {}").format(sql.Identifier(dbname)))
        '''Δημιουργεί ξανά τη βάση clinvar_db.'''
            cur.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier(dbname)))

        print("Database reset completed.")

    except Exception as e:
        print("Σφάλμα κατά τη δημιουργία βάσης:", e)

    finally:
'''Εξασφαλίζει ότι η σύνδεση θα κλείσει οπωσδήποτε, είτε προέκυψε σφάλμα είτε όχι.'''
        if conn:
            conn.close()

def download_file(
    url: str,
    output_path: str,
    timeout: tuple = (10, 30),  # (connect, read)
    max_retries: int = 3,
    backoff_factor: float = 1.5,
    chunk_size: int = 8192
) -> None:
    """Κατέβασμα αρχείου με timeout, επανάληψη και progress tracking.
    
    Args:
        url: URL προς κατέβασμα
        output_path: Τοπική διαδρομή αποθήκευσης
        timeout: (connect_timeout, read_timeout) σε δευτερόλεπτα
        max_retries: Μέγιστες επαναλήψεις σε αποτυχία
        backoff_factor: Πολλαπλασιαστής καθυστέρησης μεταξύ επαναλήψεων
        chunk_size: Μέγεθος chunk για streaming
        
    Raises:
        RuntimeError: Αν όλες οι επαναλήψεις αποτύχουν
        ValueError: Αν το αρχείο δεν είναι έγκυρο gzip
    """
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    
    for attempt in range(max_retries):
        try:
            print(f"Προσπάθεια {attempt + 1}/{max_retries}...")
            
            with requests.get(url, headers=headers, timeout=timeout, stream=True) as r:
                r.raise_for_status()
                
                # Έλεγχος status code
                if r.status_code == 404:
                    raise ValueError("404 - Το αρχείο δεν βρέθηκε")
                
                total_size = int(r.headers.get('content-length', 0))
                downloaded = 0
                
                with open(output_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            print(f"\rΚατέβασμα: {downloaded}/{total_size} bytes ({downloaded/total_size*100:.1f}%)", end='')
                
                print("\nΕπιτυχές κατέβασμα!")
                
                # Έλεγχος και αποσυμπίεση gzip
                if output_path.endswith('.gz'):
                    print("Έλεγχος gzip αρχείου...")
                    with open(output_path, 'rb') as f:
                        if f.read(2) != b'\x1f\x8b':
                            raise ValueError("Μη έγκυρο gzip αρχείο")
                    
                    print("Αποσυμπίεση...")
                    unzipped_path = output_path[:-3]
                    with gzip.open(output_path, 'rb') as f_in:
                        with open(unzipped_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    print(f"Αποσυμπιέστηκε στο: {unzipped_path}")
                
                return  # Επιτυχία
            
        except (requests.exceptions.RequestException, ValueError) as e:
            print(f"Σφάλμα: {str(e)}")
            if attempt == max_retries - 1:
                raise RuntimeError(f"Αποτυχία μετά από {max_retries} προσπάθειες")
            
            wait_time = backoff_factor * (attempt + 1)
            print(f"Αναμονή {wait_time:.1f} δευτερόλεπτα...")
            time.sleep(wait_time)





#Χρήση pretty print (indent=2)
'''Δέχεται δεδομένα είτε ως pandas.DataFrame, dict ή list.

Το filename είναι το path του αρχείου που θα αποθηκευτούν τα δεδομένα σε JSON μορφή.'''
def save_to_json(data: Union[pd.DataFrame, dict, list], filename: str) -> None:
    """
    Αποθήκευση δεδομένων σε JSON αρχείο με έλεγχο ορθότητας.
    Υποστηρίζει:
    - DataFrames
    - Λεξικά
    - Λίστες
    """
    try:
        print(f"Αποθήκευση δεδομένων στο {filename}...")
        
        # Δημιουργία φακέλου μόνο αν υπάρχει path
        folder = os.path.dirname(filename)
        if folder:
            os.makedirs(folder, exist_ok=True)
            '''Αν ο φάκελος δεν υπάρχει, τον δημιουργεί.
            Το exist_ok=True αποτρέπει σφάλμα αν ήδη υπάρχει.'''
        
        if isinstance(data, pd.DataFrame):
            #Αν είναι DataFrame:
            data.to_json(filename, orient="records", indent=2, date_format="iso")
            '''orient="records" → λίστα από λεξικά.
            indent=2 → μορφοποιημένο JSON για αναγνωσιμότητα.
            date_format="iso" → ημερομηνίες σε ISO 8601 (YYYY-MM-DDTHH:MM:SSZ).'''
        else:
            '''indent=2 → για ωραία μορφοποίηση.
            ensure_ascii=False → επιτρέπει ελληνικούς χαρακτήρες (χωρίς unicode escapes όπως \u03b1).'''
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        
        # Έλεγχος επιτυχίας
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
            print(f"Επιτυχής αποθήκευση ({os.path.getsize(filename)/1024:.1f} KB)")
        else:
            raise IOError("Το αρχείο δεν δημιουργήθηκε σωστά")
            
    except Exception as e:
        print(f"Σφάλμα κατά την αποθήκευση: {e}")
        raise
