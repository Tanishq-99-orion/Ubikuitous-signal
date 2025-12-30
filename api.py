from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import traceback

# import your pipeline function
from pipeline import run_pipeline

app = FastAPI(
    title="Ubiquitin Linkage & Process Inference API",
    version="1.0"
)

# -------------------------
# Request schema
# -------------------------
class PipelineRequest(BaseModel):
    substrate_uniprot: str
    e3_input: str
    cell_state: str              # DNA_damage | cycling | immune_active | quiescent
    cell_type: str               # one of 15 types
    top_n_e2: int = 3


# -------------------------
# API endpoint
# -------------------------
@app.post("/predict")
def predict(req: PipelineRequest):
    try:
        result = run_pipeline(
            substrate_uniprot=req.substrate_uniprot,
            e3_input=req.e3_input,
            cell_state=req.cell_state,
            cell_type=req.cell_type,
            top_n_e2=req.top_n_e2
        )

        # convert pandas DF safely
        process_df = result["process_df"].to_dict(orient="records")

        return {
            "dominant_linkage": result["dominant_linkage"],
            "linkage_probabilities": result["linkage_probs"],
            "process_inference": process_df
        }

    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
