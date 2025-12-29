import os
import json
import redis
import asyncio
from fastapi import FastAPI, HTTPException, BackgroundTasks
from pydantic import BaseModel
from typing import Optional

# Import your existing pipeline logic here
# (Assuming the _ns and run_pipeline logic is available)

app = FastAPI(title="Ubiquitin Linkage API")

# Initialize Redis (Assumes default port 6379)
# In production, use environment variables for host/port
r_cache = redis.Redis(host='localhost', port=6379, db=0, decode_responses=True)

# --- Request Models ---
class PipelineRequest(BaseModel):
    substrate_uniprot: str
    e3_input: str
    cell_state: str
    cell_type: str
    top_n_e2: Optional[int] = 3

# --- Helper for Redis Caching ---
def get_cached_result(key: str):
    data = r_cache.get(key)
    return json.loads(data) if data else None

def set_cached_result(key: str, value: dict, expire: int = 86400):
    r_cache.setex(key, expire, json.dumps(value))

# --- API Endpoints ---
@app.post("/predict")
async def predict_linkage(request: PipelineRequest):
    # Create a unique cache key based on inputs
    cache_key = f"{request.substrate_uniprot}:{request.e3_input}:{request.cell_state}:{request.cell_type}"
    
    # Check Redis Cache first
    cached = get_cached_result(cache_key)
    if cached:
        return {"status": "success", "source": "cache", "data": cached}

    try:
        # Run the pipeline in a separate thread to prevent blocking the Event Loop
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, 
            run_pipeline, 
            request.substrate_uniprot,
            request.e3_input,
            request.cell_state,
            request.cell_type,
            request.top_n_e2
        )
        
        # Convert DataFrame to JSON for storage/response
        result["process_df"] = result["process_df"].to_dict(orient="records")
        
        # Store in Redis
        set_cached_result(cache_key, result)
        
        return {"status": "success", "source": "live", "data": result}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/health")
def health_check():
    return {"status": "online", "redis": r_cache.ping()}