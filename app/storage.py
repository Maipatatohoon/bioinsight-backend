import os
import shutil
import uuid
from typing import Tuple
from fastapi import UploadFile

STORAGE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "uploads")

def ensure_storage_dir():
    os.makedirs(STORAGE_DIR, exist_ok=True)

def save_upload(file: UploadFile, prefix: str = "upload") -> str:
    """Saves an UploadFile to disk and returns the absolute file path."""
    ensure_storage_dir()
    file_id = f"{prefix}_{uuid.uuid4().hex}_{file.filename}"
    file_path = os.path.join(STORAGE_DIR, file_id)
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    return file_path

def get_file_path(filename: str) -> str:
    return os.path.join(STORAGE_DIR, filename)

def delete_file(filepath: str):
    if os.path.exists(filepath):
        os.remove(filepath)
