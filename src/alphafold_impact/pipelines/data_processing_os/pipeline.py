"""
This is a boilerplate pipeline 'data_processing_os'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import get_regional_biology_syllabii, get_topic_biology_syllabii, label_biology_courses


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
        node(
            get_regional_biology_syllabii,
            inputs={"data_loaders": "os.data_collection.regions.raw"},
            outputs="regional_biology_syllabii",
            tags=["data_processing_os"]
        ),
        node(
            get_topic_biology_syllabii,
            inputs={"data_loaders": "os.data_collection.keywords.raw"},
            outputs="topic_biology_syllabii",
            tags=["data_processing_os"]
        ),
        node(
            label_biology_courses,
            inputs={"regional_data": "regional_biology_syllabii", "keyword_data": "topic_biology_syllabii"},
            outputs="os.labelled_biology_courses.primary",
            tags=["data_processing_os"]
        ),
    ])
