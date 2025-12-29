# pipeline.py
import os
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from logic import run_pipeline

app = FastAPI(title="Ubiquitin Linkage Predictor")

# ---------- Request schema ----------
class PredictRequest(BaseModel):
    substrate_uniprot: str
    e3_input: str
    cell_state: str
    cell_type: str
    top_n_e2: int = 3

# ---------- Health check ----------
@app.get("/")
def health():
    return {"status": "ok"}

# ---------- Prediction endpoint ----------
@app.post("/predict")
def predict(req: PredictRequest):
    try:
        result = run_pipeline(
            substrate_uniprot=req.substrate_uniprot,
            e3_input=req.e3_input,
            cell_state=req.cell_state,
            cell_type=req.cell_type,
            top_n_e2=req.top_n_e2
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
