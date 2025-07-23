from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "Γεια σου από το FastAPI στο project genekor!"}
