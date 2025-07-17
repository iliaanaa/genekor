CLINVAR_README_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/README.txt"
METADATA_DIR = "metadata"
METADATA_FILE = os.path.join(METADATA_DIR, "clinvar_metadata.json")



def get_clinvar_release_date() -> str:
    """Επιστρέφει την ημερομηνία release του ClinVar σε μορφή YYYYMMDD"""
    try:
        response = requests.get(CLINVAR_README_URL, timeout=10)
        response.raise_for_status()
        
        for line in response.text.splitlines():
            if "Release date:" in line:
                parts = line.strip().split(":")
                if len(parts) == 2:
                    release_date_str = parts[1].strip()
                    # Εξασφάλιση ότι η ημερομηνία είναι σε σωστή μορφή
                    datetime.strptime(release_date_str, "%Y%m%d")
                    return release_date_str
        raise ValueError("Δεν βρέθηκε ημερομηνία release στο README")
    except Exception as e:
        raise RuntimeError(f"Σφάλμα ανάγνωσης README: {str(e)}")
    


def load_local_metadata() -> Optional[Dict]:
    """Φορτώνει τα τοπικά μεταδεδομένα αν υπάρχουν"""
    if os.path.exists(METADATA_FILE):
        try:
            with open(METADATA_FILE, "r") as f:
                data = json.load(f)
                # Επαλήθευση της δομής των μεταδεδομένων
                if "release_date" in data:
                    datetime.strptime(data["release_date"], "%Y%m%d")
                    return data
        except (json.JSONDecodeError, ValueError) as e:
            print(f"Προειδοποίηση: Άκυρα μεταδεδομένα - {str(e)}")
    return None


def save_local_metadata(release_date: str) -> None:
    """Αποθηκεύει τα μεταδεδομένα με έλεγχο ορθότητας"""
    try:
        # Επαλήθευση της ημερομηνίας
        datetime.strptime(release_date, "%Y%m%d")
        os.makedirs(METADATA_DIR, exist_ok=True)
        
        metadata = {
            "release_date": release_date,
            "last_updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        with open(METADATA_FILE, "w") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
    except ValueError as e:
        raise ValueError(f"Άκυρη μορφή ημερομηνίας: {str(e)}")
    


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
            print(f"✅ Επιτυχής αποθήκευση ({os.path.getsize(filename)/1024:.1f} KB)")
        else:
            raise IOError("❌ Το αρχείο δεν δημιουργήθηκε σωστά")
            
    except Exception as e:
        print(f"Σφάλμα κατά την αποθήκευση: {e}")
        raise




def is_after_first_thursday() -> bool:
    """Ελέγχει αν η τρέχουσα ημερομηνία είναι μετά την 1η Πέμπτη του μήνα"""
    today = datetime.now()
    first_day = today.replace(day=1)
    
    # Βρίσκουμε την 1η Πέμπτη (weekday=3)
    first_thursday = first_day
    while first_thursday.weekday() != 3:
        first_thursday += timedelta(days=1)
    
    # Επιστρέφει True μόνο αν σήμερα >= Παρασκευή μετά την 1η Πέμπτη
    return today >= (first_thursday + timedelta(days=1))




def needs_update(remote_date: str, local_date: Optional[str], force: bool = False) -> bool:
    """Κρίνει αν χρειάζεται ενημέρωση"""
    if force:
        return True
    
    if not local_date:  # Πρώτη ενημέρωση
        return True
    
    try:
        remote_dt = datetime.strptime(remote_date, "%Y%m%d")
        local_dt = datetime.strptime(local_date, "%Y%m%d")
        
        if remote_dt > local_dt:
            return True
        elif remote_dt < local_dt:
            print(f"Προσοχή: Τοπική έκδοση ({local_date}) είναι νεότερη από την remote ({remote_date})!")
            return False
        else:
            print("Η βάση είναι ήδη ενημερωμένη με την τελευταία έκδοση.")
            return False
    except ValueError as e:
        print(f"Σφάλμα ανάλυσης ημερομηνιών: {str(e)}")
        return True
    


def update_database(data_file: str) -> None:
    """Ενημερώνει τη βάση δεδομένων με τα νέα δεδομένα"""
    try:
        # Φόρτωση και επεξεργασία δεδομένων
        print("Φόρτωση και επεξεργασία δεδομένων...")
        df = pd.read_csv(data_file, sep='\t', low_memory=False)
        
        # Εδώ μπορείτε να προσθέσετε τη λογική επεξεργασίας σας
        # π.χ. φιλτράρισμα για BRCA1, μετατροπές κλπ.
        
        # Σύνδεση στη βάση δεδομένων
        conn = psycopg2.connect(**DB_CONFIG)
        
        try:
            with conn.cursor() as cur:
                # Εδώ μπορείτε να προσθέσετε τη λογική ενημέρωσης της βάσης
                # π.χ. TRUNCATE πίνακα και νέα εισαγωγή ή UPSERT
                print("Ενημέρωση βάσης δεδομένων...")
                
            conn.commit()
            print("Ενημέρωση βάσης ολοκληρώθηκε επιτυχώς!")
        except Exception as e:
            conn.rollback()
            raise e
        finally:
            conn.close()
    except Exception as e:
        raise RuntimeError(f"Σφάλμα ενημέρωσης βάσης: {str(e)}")
    finally:
        if os.path.exists(data_file):
            os.remove(data_file)

# --- Κύρια Λειτουργία ---
def main(force_update: bool = False) -> None:
    print("=== Έναρξη διαδικασίας ενημέρωσης ClinVar ===")
    
    # Έλεγχος ημερομηνίας (εκτός αν είναι forced update)
    if not force_update and not is_after_first_thursday():
        print("Η ενημέρωση επιτρέπεται μόνο μετά την 1η Πέμπτη του μήνα.")
        print("Χρησιμοποιήστε --force για να παρακάμψετε αυτόν τον έλεγχο.")
        return
    
    try:
        # Λήψη ημερομηνίας release
        remote_date = get_clinvar_release_date()
        print(f"Ημερομηνία τελευταίας έκδοσης ClinVar: {remote_date}")
        
        # Φόρτωση τοπικών μεταδεδομένων
        local_metadata = load_local_metadata()
        local_date = local_metadata.get("release_date") if local_metadata else None
        
        if local_date:
            print(f"Ημερομηνία τοπικής έκδοσης: {local_date}")
        
        # Απόφαση ενημέρωσης
        if needs_update(remote_date, local_date, force_update):
            print("Απαιτείται ενημέρωση...")
            
            # Κατέβασμα και ενημέρωση
            data_file = download_clinvar_data()
            update_database(data_file)
            
            # Αποθήκευση νέων μεταδεδομένων
            save_local_metadata(remote_date)
            print("Η ενημέρωση ολοκληρώθηκε επιτυχώς!")
        else:
            print("Δεν απαιτείται ενημέρωση.")
            
    except Exception as e:
        print(f"Κρίσιμο σφάλμα: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ενημέρωση βάσης δεδομένων ClinVar")
    parser.add_argument(
        "--force", 
        action="store_true",
        help="Εξαναγκασμός ενημέρωσης ακόμα και αν η βάση είναι ενημερωμένη"
    )
    args = parser.parse_args()
    
    try:
        main(force_update=args.force)
    except Exception as e:
        print(f"Η διαδικασία ενημέρωσης απέτυχε: {str(e)}")
        exit(1)