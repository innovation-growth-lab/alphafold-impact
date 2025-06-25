"""
Pipeline for PubMed and iCite data collection.

This pipeline handles data collection from:
- ClinicalTrials.gov (API v2)
- NIH iCite database (from figshare repository)
"""

from .pipeline import create_pipeline

__all__ = ["create_pipeline"]

__version__ = "0.1"
